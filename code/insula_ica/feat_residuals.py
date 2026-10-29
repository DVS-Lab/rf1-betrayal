"""Non-destructive FEAT inventory, reconstruction and residual-scale verification.

This module never fits a new nuisance model and never edits historical FEATs.
The original task COPEs supply the separate task-preserved analysis branch.
"""
from pathlib import Path
import html
import json
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter

import nibabel as nib
import numpy as np
import insula_ica as core


def files_identity(paths):
    return [core.identity(p, content=not str(p).endswith(('.nii', '.nii.gz'))) for p in paths if Path(p).is_file()]


def freeze_code():
    return {p.name: core.sha(p) for p in (Path(__file__), Path(core.__file__))}


def output_only_fsf(source, destination):
    text = Path(source).read_text()
    updated, n = re.subn(r'^set\s+fmri\(outputdir\)\s+.*$',
                         lambda _: f'set fmri(outputdir) "{destination}"', text, flags=re.M)
    core.require(n == 1, 'Expected exactly one original FSF output destination')
    # Every byte outside this one setting is retained.
    return updated


def film_command(feat):
    """Read the historical FILM invocation; do not invent its fit options."""
    found = []
    for path in (feat / 'report_log.html', feat / 'report.log', feat / 'logs/feat2_pre'):
        if not path.is_file():
            continue
        text = html.unescape(re.sub('<[^>]*>', '\n', path.read_text(errors='replace')))
        for line in text.splitlines():
            tokens = shlex.split(line) if 'film_gls ' in line else []
            positions = [i for i, t in enumerate(tokens) if Path(t).name == 'film_gls']
            if positions:
                cmd = tokens[positions[0]:]
                core.require(not any(t in (';', '&&', '|', '>', '2>&1') or '$' in t or '`' in t for t in cmd), 'Unsafe/ambiguous historical FILM command')
                found.append((cmd, path))
    unique = {tuple(cmd) for cmd, _ in found}
    core.require(len(unique) == 1, 'No unique stored FILM command for stats-only recovery')
    return found[0]


def option(cmd, name):
    for i, token in enumerate(cmd):
        if token.startswith(name + '='):
            return token.split('=', 1)[1]
        if token == name:
            return cmd[i+1]
    return None


def replay_command(original, target, prewhiten):
    cmd, log = film_command(original)
    core.require(option(cmd, '--in') in ('filtered_func_data', 'filtered_func_data.nii.gz', str(original/'filtered_func_data.nii.gz')), 'Historical FILM input differs from retained FEAT data')
    core.require(option(cmd, '--pd') in ('design.mat', str(original/'design.mat')), 'Historical FILM design differs')
    core.require(option(cmd, '--rn') == 'stats', 'Unexpected original FILM results directory')
    core.require(option(cmd, '--con') in ('design.con', str(original/'design.con')), 'Stored FILM command lacks verified contrast input')
    core.require(('--noest' not in cmd) == prewhiten, 'Stored FILM command and FSF disagree about prewhitening')
    core.require(not any(t.startswith(('--vef', '--ven', '--mf', '--mft', '--mode=surface')) for t in cmd), 'Unsupported voxelwise/perfusion/surface model')
    result = [shutil.which('film_gls') or 'film_gls']
    skip = False
    for i, token in enumerate(cmd[1:], 1):
        if skip:
            skip = False
            continue
        if token == '--rn':
            result.extend(['--rn', str(target/'stats')]); skip = True
        elif token.startswith('--rn='):
            result.append('--rn=' + str(target/'stats'))
        else:
            result.append(token)
    return result, log


def validate_model(settings, design, contrast, task):
    _, index, label = core.TASKS[task]
    core.require(settings.get('fmri(level)') == '1', 'Not a first-level FEAT model')
    core.require(settings.get('fmri(motionevs)', '0') == '0', 'Unexpected additional motion-EV generation')
    core.require(not settings.get('fmri(scriptevsbeta)', ''), 'Custom EV script requires separate review')
    core.require(settings[f'fmri(conname_real.{index})'] == label, 'Unexpected stored task contrast')
    n = int(settings['fmri(evs_real)'])
    expected = np.array([0,0,0,-1,1,-1,1,-1,1,0] if task == 'trust' else [0,1,0,1,0,1,0,1,0,0,0])
    weights = np.array([float(settings.get(f'fmri(con_real{index}.{j})', 0)) for j in range(1,n+1)])
    core.require(np.array_equal(weights, expected), 'Unexpected stored task contrast weights/order')
    core.require(contrast.shape[0] >= index and contrast.shape[1] == design.shape[1], 'Original contrast/design shape mismatch')
    core.require(np.array_equal(contrast[index-1,:n], expected) and np.all(contrast[index-1,n:] == 0), 'FEAT contrast matrix differs from saved FSF')
    core.require(int(settings['fmri(npts)']) - int(settings['fmri(ndelete)']) == len(design), 'Saved design/time sample count mismatch')
    core.require(int(settings['fmri(ndelete)']) == 0, 'Nonzero deletion needs an explicit timing review')
    core.require(int(settings['fmri(prewhiten_yn)']) in (0,1), 'Unknown prewhitening setting')
    return index


def safe_roots(c):
    historical=(core.ROOT/'derivatives/fsl').resolve()
    for name in ('output_root','derivatives_root'):
        location=Path(c[name]).resolve()
        core.require(location != historical and historical not in location.parents,
                     'Analysis/reconstruction root cannot be inside historical FEAT storage')
        core.require(not any(x.name.endswith('.feat') for x in (location,*location.parents)),
                     'Analysis/reconstruction root cannot be inside any FEAT directory')


def recovery_route(residual,raw,conf,events,filtered,replay):
    if residual.is_file():return 'reuse'
    if raw.is_file() and conf.is_file() and all(p.is_file() for p in events):return 'full_feat'
    if filtered.is_file() and replay is not None:return 'replay_film'
    return 'blocked'


def audit(c):
    safe_roots(c)
    out = Path(c['output_root']) / 'audit'
    core.require(not out.exists(), f'Audit exists: {out}; use a fresh output_root')
    out.mkdir(parents=True)
    sw = core.software()
    for exe in ('melodic','film_gls'):
        p = subprocess.run([exe,'--help'], text=True, capture_output=True)
        (out/f'{exe}_help.txt').write_text(p.stdout+p.stderr)
    core.require('--update_mask' in (out/'melodic_help.txt').read_text(), 'Unsupported MELODIC version')
    cohort_file = core.ROOT/'code/sublist_n132.txt'
    cohort = cohort_file.read_text().split()
    core.require(len(cohort) == len(set(cohort)) == c['expected_cohort_size'], 'Manuscript cohort mismatch')
    qc_file = Path(c['master_qc'])
    core.require(qc_file.is_file(), f'Missing saved QC table: {qc_file}')
    qc_rows = core.read_table(qc_file, ',')
    qc = {str(int(float(r['sub_id']))):r for r in qc_rows}
    core.require(len(qc)==len(qc_rows), 'Duplicate QC IDs')
    final = {s for s,r in qc.items() if r.get('include_final','').upper() in ('TRUE','1')}
    core.require(final == set(cohort), 'Corrected QC include_final must match manuscript N exactly')
    registry = files_identity([cohort_file,qc_file])
    rows, blockers = [], []
    reference_anchor=None
    for sub in cohort:
        for task,(model,_,_) in core.TASKS.items():
            for rn in (1,2):
                key = f'sub-{sub}_ses-01_task-{task}_run-{rn}'
                old = core.ROOT/f'derivatives/fsl/sub-{sub}/ses-01/L1_task-{task}_model-{model}_type-act_run-{rn}_sm-5.feat'
                target = Path(c['derivatives_root'])/'feat_reconstruction'/f'{key}.feat'
                row = {'key':key,'subject':sub,'task':task,'run':rn,'original_feat':str(old),'reconstructed_feat':str(target)}
                rows.append(row)
                flag = qc[sub].get(f'{task}_run{rn}_full_ok','').upper()
                for kind in ('beh','ev','mri'):
                    row[f'original_{kind}_ok'] = qc[sub].get(f'{task}_{kind}_run{rn}_ok','')
                if flag in ('FALSE','0'):
                    row.update(status='excluded_original_qc',reason='Saved corrected full_ok is false'); continue
                if flag not in ('TRUE','1'):
                    row.update(status='blocked',reason='Unknown original QC flag'); continue
                try:
                    fsf = old/'design.fsf'
                    if not fsf.exists():
                        fsf = old.parent/f'L1_sub-{sub}_task-{task}_model-{model}_type-act_run-{rn}.fsf'
                    core.require(fsf.is_file(), f'Missing original generated run FSF: {fsf}')
                    settings = core.fsl_settings(fsf)
                    core.require(Path(settings['fmri(outputdir)'] if settings['fmri(outputdir)'].endswith('.feat') else settings['fmri(outputdir)']+'.feat').resolve() == old.resolve(), 'Saved FSF points to a different historical output')
                    required = [old/'design.mat',old/'design.con',old/'mask.nii.gz',old/'mean_func.nii.gz',old/'stats/dof']
                    core.require(all(p.is_file() for p in required), 'Missing historical model/reference files: ' + ', '.join(str(p) for p in required if not p.is_file()))
                    dm, con = core.fsl_matrix(old/'design.mat'),core.fsl_matrix(old/'design.con')
                    index = validate_model(settings,dm,con,task)
                    stats_files = []
                    for j in range(1,con.shape[0]+1):
                        stats_files += [old/f'stats/cope{j}.nii.gz',old/f'stats/zstat{j}.nii.gz']
                    core.require(all(p.is_file() for p in stats_files), 'Historical COPE/ZSTATs incomplete; cannot prove reconstruction equivalence')
                    raw = Path(settings['feat_files(1)'])
                    conf = Path(settings.get('confoundev_files(1)',''))
                    core.require(raw.is_absolute() and conf.is_absolute(), 'Unresolved relative original BOLD/confound path')
                    internal = old/'confoundevs.txt'
                    matrix = conf if conf.is_file() else internal
                    core.require(matrix.is_file(), f'Neither original nor FEAT-preserved confound matrix exists: {conf}')
                    nuisance = core.load_numeric(matrix)
                    core.require(nuisance.shape[0]==len(dm), 'Original confound/design length mismatch')
                    if conf.is_file() and internal.is_file():
                        core.require(np.allclose(nuisance,core.load_numeric(internal),rtol=1e-6,atol=1e-6), 'External and FEAT-preserved confounds differ')
                    events = []
                    for j in range(1,int(settings['fmri(evs_orig)'])+1):
                        if settings.get(f'fmri(shape{j})') == '3':
                            path = Path(settings[f'fmri(custom{j})'])
                            core.require(path.is_absolute(), 'Relative EV path in original FSF')
                            events.append(path)
                    original_res = old/'stats/res4d.nii.gz'
                    filtered = old/'filtered_func_data.nii.gz'
                    prewhiten = settings['fmri(prewhiten_yn)'] == '1'
                    replay, replay_error = None, ''
                    try:
                        replay, log = replay_command(old,target,prewhiten)
                        registry.extend(files_identity([log]))
                    except ValueError as exc:
                        replay_error = str(exc)
                    mode=recovery_route(original_res,raw,conf,events,filtered,replay)
                    active=old if mode=='reuse' else target
                    if mode=='blocked':
                        raise ValueError('Cannot reconstruct: original res4d absent; BOLD=' + str(raw.is_file()) +
                            ', filtered_func_data=' + str(filtered.is_file()) + ', original confounds=' + str(conf.is_file()) +
                            ', EVs present=' + str(sum(p.is_file() for p in events)) + '/' + str(len(events)) + '; ' + replay_error +
                            '. Restore exact historical inputs from archive; no new-version substitute.')
                    tr = float(settings['fmri(tr)'])
                    core.require(np.isclose(tr,c['expected_tr'],atol=1e-5), 'Original model TR mismatch')
                    reference = nib.load(old/'mean_func.nii.gz')
                    if reference_anchor is None:reference_anchor=reference
                    core.require(core.same_grid(reference_anchor,reference),'Incompatible historical grids across accepted runs')
                    if original_res.is_file():
                        im=nib.load(original_res)
                        core.require(im.shape==(*reference.shape,len(dm)) and core.same_grid(im,reference), 'Existing residual grid/time mismatch')
                    if raw.is_file():
                        im=nib.load(raw)
                        core.require(im.shape==(*reference.shape,len(dm)) and core.same_grid(im,reference), 'Original input/model grid mismatch')
                    actual_dof=float((old/'stats/dof').read_text().strip())
                    core.require(actual_dof >= c['min_residual_dof']+7, 'Insufficient original FILM residual dof')
                    row.update(status='included',reason='',source_fsf=str(fsf),source_bold=str(raw),bold_exists=raw.is_file(),
                        source_confounds=str(conf),confound_matrix=str(matrix),filtered_data=str(filtered),filtered_exists=filtered.is_file(),
                        recovery_mode=mode,active_feat=str(active),residual=str(active/'stats/res4d.nii.gz'),
                        tr=tr,nvol=len(dm),nss=0,grid_shape=json.dumps(reference.shape),voxel_mm=json.dumps([float(x) for x in reference.header.get_zooms()]),
                        orientation=''.join(nib.aff2axcodes(reference.affine)),affine=json.dumps(reference.affine.tolist()),confound_columns=nuisance.shape[1],nuisance_rank=core.basis(dm).shape[1],original_film_dof=actual_dof,
                        smoothing=settings.get('fmri(smooth)'),highpass_enabled=settings.get('fmri(temphp_yn)'),
                        highpass_seconds=settings.get('fmri(paradigm_hp)'),prewhiten=int(prewhiten),contrast_index=index,
                        replay_argv=json.dumps(replay) if replay else '', duration=len(dm)*tr)
                    paths = [fsf,*required,matrix,*stats_files,*events,raw,original_res,filtered,old/'design.fts',old/'stats/threshac1.nii.gz']
                    registry.extend(files_identity(paths))
                    planned=out/key
                    planned.mkdir()
                    (planned/'reconstruction.fsf').write_text(output_only_fsf(fsf,target))
                    core.write_json(planned/'settings.json',settings)
                    core.write_json(planned/'task_contrast.json',{'index':index,'weights':con[index-1].tolist(),
                        'original_name':settings[f'fmri(conname_real.{index})']})
                except (ValueError,OSError,KeyError) as exc:
                    row.update(status='blocked',reason=str(exc))
    included=[r for r in rows if r['status']=='included']
    blocked=[r for r in rows if r['status']=='blocked']
    if blocked:blockers.append(f'{len(blocked)} accepted runs lack a safe verified recovery plan; see manifest.tsv')
    if not included:blockers.append('No recoverable accepted runs in this environment')
    pilot=[]
    for task in core.TASKS:
        rr=sorted([r for r in included if r['task']==task],key=lambda r:(r['recovery_mode'],r['nuisance_rank'],r['key']))
        if rr:
            pilot.extend([rr[0]['key'],rr[-1]['key']])
    pilot=list(dict.fromkeys(pilot))
    for r in sorted(included,key=lambda r:r['nuisance_rank'],reverse=True):
        if len(pilot)>=5:break
        if r['key'] not in pilot:pilot.append(r['key'])
    subs=sorted({r['subject'] for r in included},key=lambda s:max(r['nuisance_rank'] for r in included if r['subject']==s))
    pilot_subs=[subs[i] for i in sorted(set(np.linspace(0,len(subs)-1,min(5,len(subs))).astype(int)))] if subs else []
    core.write_table(out/'manifest.tsv',rows)
    a={'passed':not blockers,'blockers':blockers,'config':c,'software':sw,'code_hashes':freeze_code(),'inputs':registry,
        'cohort_n':len(cohort),'pilot_keys':pilot,'pilot_subjects':pilot_subs,'run_counts':dict(Counter(r['task'] for r in included)),
        'recovery_counts':dict(Counter(r['recovery_mode'] for r in included)),
        'analysis_commit':subprocess.check_output(['git','-C',str(core.ROOT),'rev-parse','HEAD'],text=True).strip()}
    core.write_json(out/'audit.json',a)
    report=['# FEAT residual input inventory',f'Cohort: {len(cohort)}. Corrected include_final identities match exactly.',
        f'Recovery plans: {a["recovery_counts"]}; task counts: {a["run_counts"]}. These are not completed residual QC.',
        f'Separate reconstruction root: {Path(c["derivatives_root"])/"feat_reconstruction"}',
        'Neither the original wrappers nor original output directories will be executed/overwritten.',
        'Existing res4d can be reused. Missing res4d can use an output-only FSF rerun if all historical sources remain, or replay the exact logged FILM command against retained FEAT filtered data/design.',
        'Prewhitened residuals require a separate scale review before ICA/background connectivity. No further nuisance regression or temporal filter is added.',
        '## Blockers\n'+'\n'.join('- '+b for b in blockers)]
    (out/'report.md').write_text('\n\n'.join(report)+'\n')
    print(out/'report.md',flush=True)
    return 0 if not blockers else 2


def inventory(c):
    out=Path(c['output_root']);a=json.loads((out/'audit/audit.json').read_text())
    core.require(a['code_hashes']==freeze_code(),'Pipeline changed since inventory; use a fresh output root')
    ignored={'network_templates','network_template_notes','residual_scale_review','ugr_pmods_increase_with_fairness','ugr_coding_evidence'}
    core.require({k:v for k,v in a['config'].items() if k not in ignored}=={k:v for k,v in c.items() if k not in ignored},'Input/reconstruction config changed since inventory')
    core.require(a['software']==core.software(),'Software changed since inventory')
    for item in a['inputs']:
        core.require(core.identity(item['path'],item['sha256'] is not None)==item,'Historical input changed: '+item['path'])
    return out,a,[r for r in core.read_table(out/'audit/manifest.tsv') if r['status']=='included']


def scope_rows(rows,a,scope):
    return [r for r in rows if r['key'] in a['pilot_keys']] if scope=='pilot' else rows


def reconstruct(c,scope,execute=False):
    out,a,rows=inventory(c);rows=scope_rows(rows,a,scope)
    safe_roots(c)
    if execute and scope=='pilot':
        core.require(3<=len(rows)<=5 and {r['task'] for r in rows}==set(core.TASKS),'Execution pilot requires 3-5 runs covering both tasks')
    if execute and scope=='full':
        check=json.loads((out/'checks/pilot.json').read_text())
        core.require(check['comparison_passed'],'Pilot historical reproduction has not passed')
    plan=[]
    for r in rows:
        if r['recovery_mode']=='reuse':
            plan.append({'key':r['key'],'mode':'reuse','argv':[],'output':r['residual']});continue
        target=Path(r['reconstructed_feat'])
        core.require(not target.exists() or (target/'stats/res4d.nii.gz').is_file(),'Partial reconstruction already exists; preserve it and inspect instead of rerunning')
        cmd=(['feat',str(out/'audit'/r['key']/'reconstruction.fsf')] if r['recovery_mode']=='full_feat' else json.loads(r['replay_argv']))
        plan.append({'key':r['key'],'mode':r['recovery_mode'],'argv':cmd,'cwd':r['original_feat'] if r['recovery_mode']=='replay_film' else str(core.ROOT),'output':r['residual']})
        if not execute or (target/'stats/res4d.nii.gz').is_file():continue
        target.parent.mkdir(parents=True,exist_ok=True)
        if r['recovery_mode']=='replay_film':
            target.mkdir()
            for name in ('design.mat','design.con','design.fts','mask.nii.gz','mean_func.nii.gz'):
                source=Path(r['original_feat'])/name
                if source.is_file():shutil.copy2(source,target/name)
            shutil.copy2(out/'audit'/r['key']/'reconstruction.fsf',target/'design.fsf')
        logdir=out/'reconstruction_logs'/r['key'];logdir.mkdir(parents=True,exist_ok=True)
        core.write_json(logdir/'command.json',plan[-1])
        with (logdir/'command.log').open('a') as log:
            result=subprocess.run(cmd,cwd=plan[-1]['cwd'],stdout=log,stderr=subprocess.STDOUT)
        core.require(result.returncode==0,'Reconstruction command failed; inspect '+str(logdir))
        core.require((target/'stats/res4d.nii.gz').is_file(),'FEAT has not completed or failed; inspect logs and do not rerun in place')
    core.write_json(out/f'reconstruction_plan_{scope}.json',plan)
    print(('Executed' if execute else 'DRY RUN:')+' '+str(out/f'reconstruction_plan_{scope}.json'),flush=True)


def whiten_design(design, autocorr):
    """Diagnostic reproduction of FILM's voxelwise FFT design transform only.

    Does not transform or re-clean any analysis residual image. Source: installed
    FSL AutoCorrEstimator::setDesignMatrix/preWhiten, tested against film_gls.
    """
    t=len(design);n=1 << (t-1).bit_length();half=t//2
    ac=np.zeros(half);ac[:min(half,len(autocorr))]=autocorr[:half]
    mirror=np.zeros(n);mirror[:half]=ac
    mirror[-(half-1):]=ac[1:][::-1]
    spectrum=np.fft.fft(mirror).real
    weights=np.zeros(n);weights[1:]=1/np.sqrt(np.abs(spectrum[1:]))
    core.require(np.isfinite(weights).all(),'Nonfinite FILM diagnostic whitening filter')
    weights/=np.sqrt(np.mean(weights**2))
    mean=design.mean(axis=0)
    transformed=np.fft.ifft(np.fft.fft(design-mean,n=n,axis=0)*weights[:,None],axis=0).real[:t]+mean
    return transformed,weights


def relative_error(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-12))


def residual_qc(feat,original,row,c):
    im=nib.load(feat/'stats/res4d.nii.gz');brain=nib.load(original/'mask.nii.gz')
    dm=core.fsl_matrix(original/'design.mat');mask=brain.get_fdata()>0
    core.require(im.shape==(*mask.shape,len(dm)) and core.same_grid(im,brain),'Residual geometry/time sample mismatch')
    core.require(np.isclose(im.header.get_zooms()[3],float(row['tr']),atol=1e-5),'Residual TR mismatch')
    core.require(im.header.get_xyzt_units()[1]=='sec','Residual temporal units are not verified seconds')
    qform,qcode=im.get_qform(coded=True);sform,scode=im.get_sform(coded=True)
    core.require(not(qcode and scode) or np.allclose(qform,sform,atol=1e-3),'Conflicting residual qform/sform')
    data=im.get_fdata(dtype=np.float32);values=data[mask].T
    core.require(np.isfinite(values).all(),'Nonfinite residuals')
    valid=values.std(axis=0)>1e-8
    core.require(valid.sum()>100,'Insufficient nonconstant residual coverage')
    valid_indices=np.flatnonzero(valid)
    indices=valid_indices[np.linspace(0,len(valid_indices)-1,min(32,len(valid_indices))).astype(int)]
    leakage=[];ordinary=[]
    ac_path=feat/'stats/threshac1.nii.gz'
    prewhiten=int(row['prewhiten'])==1
    ac=nib.load(ac_path).get_fdata()[mask] if prewhiten and ac_path.exists() else None
    if ac is not None and ac.ndim==1:ac=ac[:,None]
    core.require(not prewhiten or ac is not None,'Cannot verify prewhitened residuals without retained threshac1')
    for j in indices:
        y=values[:,j].astype(float)
        q=core.basis(dm)
        ordinary.append(float(np.linalg.norm(q.T@y)/max(np.linalg.norm(y),1e-12)))
        effective=whiten_design(dm,ac[j])[0] if prewhiten else dm
        q=core.basis(effective)
        leakage.append(float(np.linalg.norm(q.T@y)/max(np.linalg.norm(y),1e-12)))
    core.require(max(leakage)<=c['residual_projection_tolerance'],'Strong fitted-design signal remains under matching FILM convention')
    scale='voxelwise_prewhitened' if prewhiten else 'unwhitened'
    record={'key':row['key'],'scale':scale,'sampled_voxels':len(indices),'max_effective_design_projection':max(leakage),
        'max_unwhitened_design_projection':max(ordinary),'nvol':im.shape[3],'voxels':int(mask.sum())}
    # If fitted input is retained, additionally verify Y-X*beta on original scale.
    filtered=feat/'filtered_func_data.nii.gz'
    if not filtered.exists():filtered=Path(row['filtered_data'])
    if filtered.exists():
        y=nib.load(filtered).get_fdata(dtype=np.float32)[mask].T[:,indices].astype(float)
        y-=y.mean(axis=0)
        beta=np.vstack([nib.load(feat/f'stats/pe{k+1}.nii.gz').get_fdata()[mask][indices] for k in range(dm.shape[1])])
        ordinary_res=y-dm@beta
        record['unwhitened_joint_fit_vs_res4d_relative_error']=relative_error(ordinary_res,values[:,indices])
        if not prewhiten:
            core.require(record['unwhitened_joint_fit_vs_res4d_relative_error']<c['residual_projection_tolerance'],'Unwhitened joint decomposition does not reconstruct residuals')
        else:
            errors=[]
            for v,j in enumerate(indices):
                effective,w=whiten_design(dm,ac[j]);n=len(w)
                yw=np.fft.ifft(np.fft.fft(y[:,v],n=n)*w).real[:len(y)]
                errors.append(relative_error(yw-effective@beta[:,v],values[:,j]))
            record['max_whitened_decomposition_relative_error']=max(errors)
            core.require(max(errors)<c['residual_projection_tolerance'],'Native FILM residual scale could not be reconstructed numerically')
    valid_mask=np.zeros(mask.shape,dtype='uint8');valid_mask[mask]=valid
    valid_path=Path(c['output_root'])/'checks/valid_masks'/f'{row["key"]}.nii.gz'
    valid_path.parent.mkdir(parents=True,exist_ok=True)
    core.save_image(valid_path,valid_mask,brain)
    record['valid_mask']=str(valid_path)
    record['zero_variance_voxels']=int((~valid).sum())
    return record


def verify(c,scope):
    out,a,all_rows=inventory(c);rows=scope_rows(all_rows,a,scope)
    checks=[];failures=[];metrics=[]
    for r in rows:
        try:
            original=Path(r['original_feat']);active=Path(r['active_feat'])
            core.require((active/'stats/res4d.nii.gz').is_file(),'Residual not available; execute reconstruction or wait for FEAT')
            for name in ('design.mat','design.con','design.fts'):
                if (original/name).exists():
                    core.require((active/name).exists() and np.allclose(core.fsl_matrix(original/name),core.fsl_matrix(active/name),rtol=1e-8,atol=1e-8),'Rerun model differs: '+name)
            reference=nib.load(original/'mean_func.nii.gz')
            for source in sorted((original/'stats').glob('*.nii.gz')):
                if not re.fullmatch(r'(cope|zstat)\d+\.nii\.gz',source.name):continue
                target=active/'stats'/source.name
                core.require(target.exists(),'Missing rerun statistic: '+source.name)
                old,new=nib.load(source),nib.load(target)
                core.require(core.same_grid(old,new) and core.same_grid(old,reference),'Rerun/historical statistic grid differs')
                x,y=old.get_fdata(),new.get_fdata()
                close=np.allclose(x,y,atol=c['comparison_atol'],rtol=c['comparison_rtol'])
                metrics.append({'key':r['key'],'statistic':source.name,'max_absolute_difference':float(np.max(np.abs(x-y))),'relative_l2':relative_error(x,y),'passed':close})
                core.require(close,'Material historical COPE/ZSTAT difference: '+source.name)
            for name in ('mask.nii.gz','mean_func.nii.gz'):
                source,target=nib.load(original/name),nib.load(active/name)
                core.require(core.same_grid(source,target) and np.allclose(source.get_fdata(),target.get_fdata(),atol=c['comparison_atol'],rtol=c['comparison_rtol']),'Rerun reference/mask differs: '+name)
            active_settings=core.fsl_settings(active/'design.fsf')
            original_settings=core.fsl_settings(r['source_fsf'])
            for key in ('tr','npts','ndelete','smooth','temphp_yn','paradigm_hp','prewhiten_yn'):
                core.require(float(active_settings[f'fmri({key})'])==float(original_settings[f'fmri({key})']),'Rerun metadata differs: '+key)
            checks.append(residual_qc(active,original,r,c))
        except (ValueError,OSError,KeyError) as exc:
            failures.append({'key':r['key'],'reason':str(exc)})
    core.write_table(out/f'checks/{scope}_statistics.tsv',metrics)
    core.write_table(out/f'checks/{scope}_residual_qc.tsv',checks)
    core.write_table(out/f'checks/{scope}_failures.tsv',failures)
    reviewed=c.get('residual_scale_review',{})
    scale_ok=all(r['scale']=='unwhitened' or (reviewed.get('accept_voxelwise_prewhitened') is True and reviewed.get('reason')) for r in checks)
    record={'comparison_passed':not failures and bool(rows),'scale_approved':bool(checks) and scale_ok,
        'verified_inputs':files_identity([path for r in rows for path in (Path(r['active_feat'])/'stats').glob('*.nii.gz')]+[Path(x['valid_mask']) for x in checks]),
        'rows':len(rows),'checks':checks,'failures':failures,'input_inventory_sha256':core.sha(out/'audit/audit.json')}
    core.write_json(out/f'checks/{scope}.json',record)
    print(f'{scope}: historical comparison {"PASS" if record["comparison_passed"] else "BLOCKED"}; residual scale {"APPROVED" if record["scale_approved"] else "REVIEW REQUIRED"}',flush=True)
    return 0 if record['comparison_passed'] else 2


def checked_audit(c):
    out,a,rows=inventory(c)
    core.require(a['passed'],'Some original accepted runs have no recoverable source; inspect inventory')
    check=json.loads((out/'checks/full.json').read_text())
    core.require(check['comparison_passed'],'Full historical reproduction/residual QC has not passed')
    review=c.get('residual_scale_review',{})
    core.require(all(r['scale']=='unwhitened' or (review.get('accept_voxelwise_prewhitened') is True and review.get('reason')) for r in check['checks']),
                 'Prewhitened res4d requires explicit scientific scale review before connectivity')
    for item in check['verified_inputs']:
        core.require(core.identity(item['path'],False)==item,'Verified residual changed; rerun QC')
    for r in rows:
        active=Path(r['active_feat'])
        core.require((active/'stats/res4d.nii.gz').exists(),'Verified residual disappeared')
    return out,a,rows


def prepare(c):
    out,a,rows=checked_audit(c)
    destination=out/'prepared';core.require(not destination.exists(),'Prepared data already exist')
    destination.mkdir()
    ref=nib.load(Path(rows[0]['original_feat'])/'mean_func.nii.gz')
    common=np.ones(ref.shape,dtype=bool)
    for r in rows:
        original=Path(r['original_feat']);mask=nib.load(original/'mask.nii.gz')
        res=nib.load(r['residual'])
        core.require(core.same_grid(ref,mask) and core.same_grid(ref,res),'Incompatible run geometry/coverage')
        valid=nib.load(out/'checks/valid_masks'/f'{r["key"]}.nii.gz').get_fdata()>0
        common &= (mask.get_fdata()>0)&valid
    fsldir=Path(__import__('os').environ['FSLDIR'])
    atlas=fsldir/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
    xml=fsldir/'data/atlases/HarvardOxford-Cortical.xml'
    masks,counts,label=core.anatomical_masks(atlas,xml,ref)
    coverage={}
    for roi in core.ROIS:
        final=masks[roi]&common;lost=masks[roi]&~common
        xyz=nib.affines.apply_affine(ref.affine,np.argwhere(lost))
        coverage[roi]=dict(counts[roi],final_voxels=int(final.sum()),retained_fraction=float(final.sum()/masks[roi].sum()),label=label,
                          lost_mni_bounds=[xyz.min(axis=0).tolist(),xyz.max(axis=0).tolist()] if len(xyz) else [])
        core.require(coverage[roi]['retained_fraction']>=c['minimum_roi_coverage'],'Incompatible '+roi+' coverage')
        core.save_image(out/f'audit/{roi}_mask.nii.gz',final.astype('uint8'),ref)
        core.save_image(out/f'audit/{roi}_lost_coverage.nii.gz',lost.astype('uint8'),ref)
    core.save_image(out/'audit/wholebrain_coverage.nii.gz',common.astype('uint8'),ref)
    shutil.copy2(Path(rows[0]['original_feat'])/'mean_func.nii.gz',destination/'reference_mean.nii.gz')
    prepared=[]
    for r in rows:
        # Reference the verified native residual; do not rescale/filter/re-clean it.
        prepared.append({'key':r['key'],'residual':r['residual'],**core.identity(r['residual'],content=False)})
    core.write_table(destination/'residual_inputs.tsv',prepared)
    core.write_json(destination/'coverage.json',coverage)
    core.write_json(destination/'complete.json',{'inventory_sha256':core.sha(out/'audit/audit.json'),'native_residuals_unchanged':True,
        'source_atlas':core.identity(atlas),'source_xml':core.identity(xml)})


def regress(c,scope):
    out,a,rows=checked_audit(c);rows=core.scoped_rows(rows,a,scope);dest=out/scope
    core.require((dest/'ica_complete.json').exists(),'ICA not complete')
    wb_image=nib.load(out/'audit/wholebrain_coverage.nii.gz');wb=wb_image.get_fdata()>0
    diagnostics=[]
    for roi in core.ROIS:
        mask=nib.load(out/f'audit/{roi}_mask.nii.gz').get_fdata()>0
        maps=nib.load(dest/f'{roi}/dim-5/oriented_IC.nii.gz').get_fdata()[mask]
        spatial=np.column_stack([np.ones(mask.sum()),maps-maps.mean(axis=0)])
        core.require(np.linalg.matrix_rank(spatial)==6,'Nonestimable spatial ICA basis')
        weights=np.linalg.pinv(spatial)[1:]
        np.savetxt(dest/roi/'spatial_projection_weights.tsv',weights,delimiter='\t')
        for row in rows:
            rd=dest/roi/'runs'/row['key'];rd.mkdir(parents=True,exist_ok=True)
            residual=nib.load(row['residual']).get_fdata(dtype=np.float32)
            tc=(weights@residual[mask]).T
            scale=tc.std(axis=0,ddof=1);core.require((scale>1e-8).all(),'Zero residual component variance')
            standardized=tc/scale
            # Both stages use exactly the same residuals. No task/nuisance EVs or
            # temporal filter are appended to the five-component temporal model.
            core.require(np.linalg.matrix_rank(standardized)==5,'Residual component temporal basis rank deficient')
            condition=float(np.linalg.cond(standardized))
            core.require(condition<c['max_condition_number'],'Ill-conditioned residual temporal basis')
            coefficients=np.linalg.pinv(standardized)@residual[wb].T
            image=np.zeros((*wb.shape,5),dtype='float32');image[wb]=coefficients.T
            core.save_image(rd/'partial_connectivity.nii.gz',image,wb_image)
            np.savetxt(rd/'timecourses.tsv',standardized,delimiter='\t')
            core.write_json(rd/'scales.json',{'residual_tc_sd':scale.tolist(),'unit':'native FILM residual unit per TC SD'})
            diagnostics.append({'roi':roi,'key':row['key'],'residual_scale': 'prewhitened' if row['prewhiten']=='1' else 'unwhitened',
                'temporal_condition':float(np.linalg.cond(standardized)),'no_additional_temporal_regression':True})
    core.write_table(dest/'regression_qc.tsv',diagnostics)
    core.write_json(dest/'regression_complete.json',{'scope':scope,
        'source':'native FEAT residual-to-residual background connectivity only','n_runs':len(rows)})

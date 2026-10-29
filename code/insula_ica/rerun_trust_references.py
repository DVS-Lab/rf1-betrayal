#!/usr/bin/env python3
"""Redo affected Trust runs with the unchanged stock filter and cached anatomy."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import threading
import traceback
import uuid

import insula_ica as core
import rebuild_trust as rebuild
from diagnose_trust_registration import sbref_selection, SOURCE
from export_trust_qc import montage
from regenerate_trust import qc_run
from run_bold_controls import save, publish


STOCK = core.ROOT / 'code/insula_ica/records/2026-10-28-stock-filter-rerun/sources/fmriprep_config.json'
RAW = Path('/ZPOOL/data/projects/rf1-sra-linux2/bids')
WORK = Path('/ZPOOL/data/scratch') / os.environ.get('USER','tug87422') / 'rf1-trust-stock-filter'
IMAGE = Path('/ZPOOL/data/tools/fmriprep-24.1.1.simg')
LOCK = threading.Lock()

# Uses the installed sMRIPrep API and builds its anatomical graph without running
# it. Fail closed if native anatomy/TPMs/transforms are absent or a heavy anatomy
# interface remains. No subject starts processing without this check.
ANATOMY_PREFLIGHT = r'''
import inspect, json, sys
from pathlib import Path
import importlib.metadata as md
import nibabel as nib
import numpy as np
from smriprep.utils.bids import collect_derivatives
from smriprep.workflows.anatomical import init_anat_fit_wf
from niworkflows.utils.spaces import SpatialReferences, Reference
subject, target = sys.argv[1:]
assert md.version('fmriprep') == '24.1.1', 'Wrong fMRIPrep version'
spaces = ['MNI152NLin6Asym', 'MNI152NLin2009cAsym']
cache = collect_derivatives('/rebuild/anatomy-cache', subject.removeprefix('sub-'), spaces)
for key in ('t1w_preproc','t1w_mask','t1w_dseg','t1w_tpms'):
    assert key in cache, 'Missing cached anatomy: '+key
for key in ('t1w_preproc','t1w_mask','t1w_dseg'):
    assert isinstance(cache[key], str) and Path(cache[key]).is_file(), 'Ambiguous/missing '+key
assert isinstance(cache['t1w_tpms'], list) and len(cache['t1w_tpms']) == 3, 'Need three tissue probability maps'
native = nib.load(cache['t1w_preproc'])
for path in [cache['t1w_mask'],cache['t1w_dseg'],*cache['t1w_tpms']]:
    image = nib.load(path)
    assert image.shape == native.shape and np.allclose(image.affine,native.affine), 'Cached anatomy grids disagree'
for space in spaces:
    transforms = cache.get('transforms',{}).get(space,{})
    assert set(transforms)=={'forward','reverse'}, 'Missing cached transforms: '+space
    for path in transforms.values():
        assert isinstance(path,str) and Path(path).is_file(), 'Ambiguous/missing cached transform'
raw_t1w = sorted(str(p) for p in Path('/bids',subject,'ses-01/anat').glob('*_T1w.nii*'))
assert raw_t1w, 'Missing raw T1w in selected-session view'
kwargs = dict(bids_root='/bids',output_dir='/rebuild/fmriprep-24',freesurfer=False,
    hires=True,longitudinal=False,msm_sulc=False,t1w=raw_t1w,t2w=[],flair=[],
    skull_strip_mode='auto',skull_strip_template=Reference.from_string('OASIS30ANTs')[0],
    spaces=SpatialReferences(spaces),precomputed=cache,omp_nthreads=8)
signature = inspect.signature(init_anat_fit_wf)
required = {k for k,v in signature.parameters.items() if v.default is inspect.Parameter.empty
            and v.kind not in (v.VAR_POSITIONAL,v.VAR_KEYWORD)}
assert required <= set(kwargs), 'Unexpected installed anatomical API: '+str(required-set(kwargs))
wf = init_anat_fit_wf(**{k:v for k,v in kwargs.items() if k in signature.parameters})
nodes = sorted((n.fullname,type(n.interface).__name__) for n in wf._get_all_nodes())
heavy = {'Registration','RegistrationSynQuick','SpatialNormalization','RobustMNINormalization',
         'BrainExtraction','N4BiasFieldCorrection','StructuralReference','FAST'}
for name, interface in nodes:
    assert interface not in heavy, 'Anatomical processing would rerun: '+name+' '+interface
Path(target).write_text(json.dumps(dict(fmriprep=md.version('fmriprep'),smriprep=md.version('smriprep'),
    cache=cache,anatomical_graph=nodes,anatomical_registration_nodes_absent=True),indent=2)+'\n')
'''


def choose_runs(rows, audit, repaired_keys):
    by_key = {r['key']:r for r in rows}
    unknown = [a['key'] for a in audit if a['status']=='unknown_no_retained_log_evidence' and a['key'] not in repaired_keys]
    core.require(not unknown, 'Reference selection remains unknown: '+', '.join(unknown))
    selected = [dict(by_key[a['key']]) for a in audit if a['status']=='different_run_reference' and a['key'] not in repaired_keys]
    core.require(len(selected)==108 and all(int(r['run'])==2 for r in selected), 'Expected 108 remaining affected Trust run-2 inputs')
    core.require(len({r['subject'] for r in selected})==108, 'Expected one affected run per participant')
    return selected


def link_file(source, destination):
    core.require(source.is_file(), 'Missing source: '+str(source))
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.symlink_to(source.resolve())


def stage_inputs(row, case, raw=None, stock=None, task='trust'):
    """Limit the input view, not the JSON queries; expose only anatomy for reuse."""
    raw=RAW if raw is None else raw
    stock=STOCK if stock is None else stock
    case.mkdir(parents=True,exist_ok=False)
    bids = case/'bids-view'; bids.mkdir()
    for source in raw.iterdir():
        if source.is_file(): link_file(source,bids/source.name)
    session = raw/row['subject']/'ses-01'
    view = bids/row['subject']/'ses-01'; view.mkdir(parents=True)
    for name in ('anat','fmap'):
        core.require((session/name).is_dir(),'Missing raw '+name)
        (view/name).symlink_to((session/name).resolve(),target_is_directory=True)
    for source in (raw/row['subject']).iterdir():
        if source.is_file():link_file(source,view.parent/source.name)
    for source in session.iterdir():
        if source.is_file():link_file(source,view/source.name)
    func=view/'func';func.mkdir()
    for name in row['raw_images']:
        source=Path(name)
        core.require(source.parent.resolve()==(session/'func').resolve(),'Unexpected raw functional root')
        link_file(source,func/source.name)
    # Keep all task SBRefs. The stock query inherits the BOLD run and selects
    # the matching one. Include JSON inheritance without adding other BOLD data.
    for source in (session/'func').iterdir():
        if source.is_file() and (source.suffix=='.json' or (f'_task-{task}_' in source.name and source.name.endswith(('_sbref.nii.gz','_sbref.nii')))):
            if not (func/source.name).exists():link_file(source,func/source.name)
    config=case/'fmriprep_config.json'; shutil.copyfile(stock,config)
    core.require(config.read_bytes()==stock.read_bytes(),'Stock filter copy changed')
    source_anat=(Path(row['new_bold']).parent.parent/'anat').resolve()
    cache=case/'anatomy-cache'; (cache/row['subject']/'ses-01/anat').mkdir(parents=True)
    registry=[]
    for source in sorted(source_anat.iterdir()):
        if source.is_file():
            link_file(source,cache/row['subject']/'ses-01/anat'/source.name)
            registry.append(core.identity(source,content=False))
    core.require(registry,'No cached anatomy')
    (cache/'dataset_description.json').write_text(json.dumps(dict(Name='Existing fMRIPrep 24.1.1 anatomy only',
        BIDSVersion='1.9.0',DatasetType='derivative',GeneratedBy=[dict(Name='fMRIPrep',Version='24.1.1')]))+'\n')
    return source_anat,registry


def command(runtime, case, work, subject, source_anat):
    return [runtime,'run','--cleanenv','-B',f'{case}/bids-view:/bids:ro',
        '-B',f'{RAW}:{RAW}:ro','-B',f'{source_anat}:{source_anat}:ro',
        '-B',f'{case}:/rebuild','-B',f'{work}:/scratch',
        '-B','/ZPOOL/data/tools/templateflow:/opt/templateflow',
        '-B','/ZPOOL/data/tools/mplconfigdir:/opt/mplconfigdir',
        '-B','/ZPOOL/data/tools/licenses:/opts:ro',str(IMAGE),'/bids','/rebuild/fmriprep-24',
        'participant','--participant-label',subject.removeprefix('sub-'),'--task-id','trust',
        '--stop-on-first-crash','--skip-bids-validation','--nthreads','14','--omp-nthreads','8','--mem-mb','28000',
        '--me-output-echos','--output-spaces','T1w','MNI152NLin6Asym',
        '--bids-filter-file','/rebuild/fmriprep_config.json','--derivatives','anatomy=/rebuild/anatomy-cache',
        '--fs-no-reconall','--fs-license-file','/opts/fs_license.txt','--notrack','-w','/scratch']


def execute(argv, env, log):
    """Stop before execution if actual workflow construction asks for new ANTs normalization."""
    with log.open('w') as stream:
        process=subprocess.Popen(argv,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                 text=True,start_new_session=True)
        try:
            for line in process.stdout:
                stream.write(line);stream.flush()
                if 'ANAT Stage 4: Preparing normalization workflow' in line:
                    raise ValueError('Actual workflow attempted new anatomical normalization; job stopped')
            core.require(process.wait()==0,'fMRIPrep failed; see '+str(log))
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL);process.wait()
            raise
        finally:
            process.stdout.close()


def run(record,state,stamp,only_subjects=None):
    source=core.ROOT/'code/insula_ica/records'/SOURCE
    rows=json.loads((source/'input-registry.json').read_text())
    repair=core.ROOT/'code/insula_ica/records/linux2-registration-fix-sub-10834-20261028T151154Z-acdf91'
    decision=json.loads((repair/'review-decision.json').read_text())
    core.require(decision['assessment']=='accept_repaired_pair_for_sub_10834_analysis_use','Missing reviewed sub-10834 repair')
    extra=core.ROOT/'code/insula_ica/records/linux2-bold-controls-20261027T113051Z-b1efa8/sub-10478-v24-a.fmriprep-log.txt'
    audit=json.loads((repair/'cohort-sbref-selection-audit.json').read_text())
    unknown_keys={a['key'] for a in audit if a['status']=='unknown_no_retained_log_evidence'}
    resolved={a['key']:a for a in sbref_selection([r for r in rows if r['key'] in unknown_keys],[extra])}
    audit=[resolved.get(a['key'],a) if a['key'] in unknown_keys else a for a in audit]
    save(record/'reference-selection-audit.json',audit)
    selected=choose_runs(rows,audit,{r['key'] for r in decision['selected_inputs']})
    if only_subjects:
        wanted={s if s.startswith('sub-') else 'sub-'+s for s in only_subjects}
        core.require(wanted <= {r['subject'] for r in selected},'Requested subject is outside the affected set')
        selected=[r for r in selected if r['subject'] in wanted]
    # Read the standard filter on linux2 and preserve its exact bytes, unchanged.
    actual=RAW.parent/'code/fmriprep_config.json'
    core.require(json.loads(actual.read_text())==json.loads(STOCK.read_text()),'Linux2 stock filter differs from the verified wrapper configuration')
    shutil.copyfile(actual,record/'fmriprep_config.json')
    original_script=RAW.parent/'code/fmriprep.sh'
    shutil.copyfile(original_script,record/'original-fmriprep.sh')
    save(record/'policy.json',dict(stock_filter_source=str(actual),stock_filter_sha256=core.sha(actual),
        original_script_sha256=core.sha(original_script),
        anatomical_reuse_required=True,concurrent_subjects=8,accepted_repair=str(repair),
        notes='24.1.1 volume settings retained. Original fmriprep.sh/JSON are not edited. Input/run scope is supplied by a read-only BIDS view, not altered queries.'))
    runtime=shutil.which('apptainer') or shutil.which('singularity')
    core.require(runtime and shutil.which('feat') and os.environ.get('FSLDIR'),'Normal container/FSL setup required')
    expected=json.loads((source/'container.json').read_text())
    core.require(rebuild.sha(IMAGE)==expected['sha256'],'24.1.1 container changed')
    root=core.ROOT/'derivatives/insula_ica'/('trust-stock-filter-'+stamp);root.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ)
    for prefix in ('APPTAINERENV_','SINGULARITYENV_'):
        env.update({prefix+'TEMPLATEFLOW_HOME':'/opt/templateflow',prefix+'MPLCONFIGDIR':'/opt/mplconfigdir'})
    state.update(state='running',selected_runs=len(selected),output_root=str(root));save(record/'state.json',state)
    for name in ('qc','cases','errors'):(record/name).mkdir()

    def launch(original):
        key=original['key'];subject=original['subject'];case=root/subject;work=WORK/root.name/subject
        result=dict(key=key,subject=subject,state='failed')
        try:
            core.require(core.sha(original['confounds'])==original['confounds_sha256'],'Preserved confounds changed')
            warpkit=rebuild.warpkit_inputs(RAW.parent,subject,int(original['run']),[Path(p) for p in original['raw_images']])
            core.require(warpkit['ready'],'Existing Warpkit inputs are incomplete: '+subject)
            source_anat,identities=stage_inputs(original,case,stock=actual)
            work.mkdir(parents=True,exist_ok=False)
            argv=command(runtime,case,work,subject,source_anat)
            save(record/'cases'/(subject+'.json'),dict(state='preflight',argv=argv,original_row=original,
                anatomical_inputs=identities,existing_warpkit=warpkit))
            # Same bind table as the run, but inspect the installed cache/graph.
            executable=argv[:argv.index(str(IMAGE))]+[str(IMAGE),'python','-c',ANATOMY_PREFLIGHT,subject,'/rebuild/anatomy-preflight.json']
            executable[1]='exec'
            with (case/'anatomy-preflight.log').open('w') as stream:
                subprocess.run(executable,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
            preflight=json.loads((case/'anatomy-preflight.json').read_text())
            core.require(preflight['anatomical_registration_nodes_absent'],'Anatomical reuse not established')
            print('Running functional rerun:',subject,'| cached anatomy verified',flush=True)
            log=case/(subject+'.fmriprep.log');execute(argv,env,log)
            text=log.read_text(errors='replace')
            core.require('ANAT Stage 4: Found pre-computed registrations' in text,'Missing actual normalization-reuse log evidence')
            row=dict(original,new_bold=str(case/'fmriprep-24'/subject/'ses-01/func'/Path(original['new_bold']).name))
            rebuild.verify_outputs([row]);rebuild.fit_run(row,case,strict_historical=False)
            qc=qc_run(row,case);save(record/'qc'/(key+'.json'),qc)
            montage(qc,record/'qc'/(key+'.png'),Path(os.environ['FSLDIR']))
            used=sbref_selection([row])[0]
            core.require(used['status'] in ('same_run','bold_derived_reference_no_sbref'),
                         'Stock reference selection not established from actual run log: '+used['status'])
            core.require(identities==[core.identity(x['path'],content=False) for x in identities],'Cached anatomical sources changed')
            result.update(state='complete_qc_review_pending',row=row,reference_selection=used,
                          anatomical_reuse_verified=True,review_flags=qc['review_flags'])
            save(record/'cases'/(subject+'.json'),dict(argv=argv,anatomy_preflight=preflight,
                anatomical_inputs=identities,existing_warpkit=warpkit,**result))
        except Exception as exc:
            result['error']=str(exc);(record/'errors'/(subject+'.txt')).write_text(traceback.format_exc())
            for name in ('anatomy-preflight.log',subject+'.fmriprep.log'):
                path=case/name
                if path.is_file(): (record/'errors'/(subject+'-'+name+'.tail.txt')).write_text('\n'.join(path.read_text(errors='replace').splitlines()[-35:])+'\n')
        with LOCK:
            state['subjects'][subject]=result;save(record/'state.json',state)
        print(subject+': '+result['state']+((' — '+result['error']) if 'error' in result else ''),flush=True)

    with ThreadPoolExecutor(max_workers=8) as pool:
        for future in as_completed([pool.submit(launch,r) for r in selected]):future.result()
    failed=[x for x in state['subjects'].values() if x['state']=='failed']
    state.update(state='complete_with_errors' if failed else 'functional_reruns_complete_qc_pending',
                 completed=len(selected)-len(failed),failed=len(failed),finished_utc=datetime.now(timezone.utc).isoformat())
    (record/'README.md').write_text('# Trust reference rerun with unchanged stock filter\n\n'+
        f'Completed: {len(selected)-len(failed)}/{len(selected)}; failed: {len(failed)}.\n\n'+
        'Native anatomy and both MNI transform pairs were mandatory precomputed inputs. The installed anatomical graph and actual run logs were checked. '
        'Functional reference selection, registration, preprocessing and FEAT fitting were redone in new directories. '
        'Original confounds/models preserved; stock fmriprep_config.json unchanged. QC images are in qc/. No ICA or exclusions were performed.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish',action='store_true')
    parser.add_argument('--subjects',nargs='+',help='Optional affected participant subset for a failed-job retry')
    args=parser.parse_args()
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    record=core.ROOT/'code/insula_ica/records'/('linux2-trust-stock-filter-'+stamp);record.mkdir(parents=True)
    state=dict(state='preflight',subjects={},started_utc=datetime.now(timezone.utc).isoformat(),
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=core.ROOT,text=True).strip())
    try:run(record,state,stamp,args.subjects)
    except Exception as exc:
        state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
    save(record/'state.json',state)
    if args.publish:
        try:publish(record,'Record Trust functional rerun with stock filter and verified anatomical reuse')
        except Exception as exc:
            (record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
    print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='functional_reruns_complete_qc_pending' else 1)


if __name__=='__main__':main()

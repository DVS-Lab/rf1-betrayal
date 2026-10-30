#!/usr/bin/env python3
"""Read-only original PPI seed, physiological-EV and analysis-chain audit."""
import argparse
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile
import traceback
import uuid

import nibabel as nib
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/insula_ica'))
import insula_ica as core
from run_bold_controls import save,publish

RECORDS=ROOT/'code/ppi_seed_audit/records'
SEEDLIKE=re.compile(r'ains|insula|(?:^|[_-])ai(?:[_-]|$)',re.I)
EVIDENCE=re.compile(r'fslmeants|seed[-_](?:AIns|aIns|AI)|\bMASK=|\bPHYS=',re.I)
REGISTRY=ROOT/'code/insula_ica/records/linux2-network-timecourses-20261030T012523Z-613cec/accepted-inputs.json'


def sha(data):return hashlib.sha256(data).hexdigest()


def table(path,rows):
    core.write_table(path,rows)
    path.write_text(path.read_text(),newline='\n')


class Archive:
    """Exact small source bytes; split archives stay below GitHub file limits."""
    def __init__(self,record):
        self.record=record;self.files=[];self.by_path={};self.blobs={};self.number=0;self.size=0;self.tar=None

    def close(self):
        if self.tar:self.tar.close();self.tar=None

    def add_bytes(self,data,name,role,**meta):
        digest=sha(data);entry=dict(source=name,role=role,sha256=digest,size=len(data),**meta)
        if digest not in self.blobs:
            if self.tar is None or self.size+len(data)>32*1024**2:
                self.close();self.number+=1;self.size=0
                self.tar=tarfile.open(self.record/f'source-files-{self.number:03d}.tar.gz','w:gz')
            member='files/'+digest;info=tarfile.TarInfo(member);info.size=len(data);info.mode=0o600;info.mtime=0
            self.tar.addfile(info,io.BytesIO(data));self.size+=len(data)
            self.blobs[digest]=dict(archive=f'source-files-{self.number:03d}.tar.gz',member=member)
        entry.update(self.blobs[digest]);self.files.append(entry);return entry

    def capture(self,path,role):
        path=Path(path);key=str(path.resolve())
        if key in self.by_path:return self.by_path[key]
        if not path.is_file():
            entry=dict(source=str(path),role=role,exists=False);self.files.append(entry);self.by_path[key]=entry;return entry
        stat=path.stat()
        if stat.st_size>8*1024**2:
            entry=dict(source=str(path),role=role,exists=True,archived=False,size=stat.st_size,reason='Source exceeds 8 MiB small-file limit')
            self.files.append(entry);self.by_path[key]=entry;return entry
        data=path.read_bytes();after=path.stat()
        core.require((stat.st_size,stat.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'Source changed while reading: '+str(path))
        entry=self.add_bytes(data,str(path),role,exists=True,mtime_ns=stat.st_mtime_ns)
        self.by_path[key]=entry;return entry


def locate(value):
    """Keep original spelling; allow only explicit same-project path relocation."""
    original=Path(value)
    if original.is_file() or original.is_dir():return original,'as_saved'
    if not original.is_absolute():
        candidate=ROOT/original
        if candidate.exists():return candidate,'project_relative'
    for anchor in ('/derivatives/fsl/','/masks/'):
        if anchor in value:
            candidate=ROOT/anchor.lstrip('/')/value.split(anchor,1)[1]
            if candidate.exists():return candidate,'same_project_suffix_relocation'
    return original,'missing_as_saved'


def affine_check(x,y):
    x=np.asarray(x,dtype=float).ravel();y=np.asarray(y,dtype=float).ravel()
    if len(x)!=len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():return dict(status='length_or_nonfinite_mismatch')
    xc=x-x.mean();yc=y-y.mean();den=np.dot(xc,xc)
    if den<=1e-20 or np.linalg.norm(yc)<=1e-10:return dict(status='constant_series')
    slope=float(np.dot(xc,yc)/den);residual=float(np.linalg.norm(yc-slope*xc)/np.linalg.norm(yc))
    return dict(status='matches_up_to_positive_centering_scaling' if slope>0 and residual<1e-6 else 'does_not_match_supported_model',
        spatial_or_temporal_r=float(np.corrcoef(x,y)[0,1]),positive_scaling=slope,relative_residual=residual,
        note='Agreement with saved model column supports use of this current external series; does not identify a mask by itself.')



def model_column_agreement(series,column):
    """Account for feat_model shape-2 input parsing via fscanf("%f", float*).

    Keep the strict residual threshold and record both representations. This is
    only a model-input check; raw mask-mean identity remains a separate test.
    """
    direct=affine_check(series,column)
    if direct['status']!='does_not_match_supported_model':return direct
    rounded=affine_check(np.asarray(series,dtype=np.float32).astype(float),column)
    if rounded['status']!='matches_up_to_positive_centering_scaling':return direct
    return dict(rounded,status='matches_after_float32_EV_input_rounding',
        raw_text_comparison=direct,input_representation='IEEE float32 custom shape-2 EV input',
        note='Saved column agrees after FSL custom-EV input rounding and positive centering/scaling; spatial mask identity is not established.')


def supported_model_match(result):
    return result.get('status') in ('matches_up_to_positive_centering_scaling','matches_after_float32_EV_input_rounding')


def check_model(settings,j,series,design):
    unsupported=[]
    if settings.get(f'fmri(shape{j})')!='2':unsupported.append('phys EV is not one-value-per-volume custom input')
    if settings.get(f'fmri(convolve{j})','0')!='0':unsupported.append('phys EV convolution enabled')
    if any(settings.get(f'fmri({k})','0')!='0' for k in ('temphp_yn','templp_yn','ndelete')):unsupported.append('temporal filtering or deleted volumes')
    if any(int(settings.get(f'fmri(convolve{k})','0')) not in (0,1,2,3,8) for k in range(1,j)):unsupported.append('basis expansion before physiological EV')
    if any(float(v)!=0 for k,v in settings.items() if re.fullmatch(rf'fmri\(ortho{j}\.\d+\)',k)):unsupported.append('phys EV orthogonalization')
    if unsupported:return dict(status='not_checked_unsupported_processing',reasons=unsupported)
    column=j-1+sum(int(settings.get(f'fmri(deriv_yn{k})','0')) for k in range(1,j))
    if design.ndim!=2 or column>=design.shape[1]:return dict(status='model_column_unavailable')
    return dict(model_column_agreement(series,design[:,column]),real_column_zero_based=column)


def source_excerpt(path,results,limit=30):
    try:
        with path.open('rb') as stream:
            size=path.stat().st_size;start=max(0,size-4*1024**2);stream.seek(start);data=stream.read().decode(errors='replace')
        for number,line in enumerate(data.splitlines(),1):
            if EVIDENCE.search(line):
                results.append(dict(source=str(path),line_number=number if start==0 else None,tail_byte_offset=start,text=line[:1200]))
                limit-=1
                if limit==0:break
    except (OSError,PermissionError) as exc:results.append(dict(source=str(path),error=str(exc)))


def script_history(archive,record):
    relevant=[];excerpts=[]
    for folder in (ROOT/'code',ROOT/'code/archive',ROOT/'code/ppi_audit',ROOT/'templates',ROOT/'templates/archive'):
        if not folder.is_dir():continue
        for path in sorted(folder.iterdir()):
            if path.is_file() and path.suffix in ('.sh','.py','.fsf','.md') and re.search(r'L[123]stats|ppi|AIns|aIns|intersection',path.name,re.I):
                relevant.append(path);archive.capture(path,'current_script_or_template');source_excerpt(path,excerpts)
    paths=[str(p.relative_to(ROOT)) for p in relevant if p.parent in (ROOT/'code',ROOT/'code/archive') and p.suffix in ('.sh','.py')]
    history=[]
    if paths:
        log=subprocess.check_output(['git','log','--all','--max-count=100','--format=COMMIT %H %aI %s','--name-only','--',*paths],cwd=ROOT,text=True)
        commit=None;seen=set()
        for line in log.splitlines():
            if line.startswith('COMMIT '):
                parts=line.split(' ',3);commit=dict(commit=parts[1],date=parts[2],subject=parts[3] if len(parts)>3 else '')
            elif line and commit:
                proc=subprocess.run(['git','show',commit['commit']+':'+line],cwd=ROOT,capture_output=True)
                if proc.returncode:continue
                digest=sha(proc.stdout);history.append(dict(commit,path=line,sha256=digest))
                if (line,digest) in seen:continue
                seen.add((line,digest));archive.add_bytes(proc.stdout,line,'historical_script',**commit)
                for number,text in enumerate(proc.stdout.decode(errors='replace').splitlines(),1):
                    if EVIDENCE.search(text):excerpts.append(dict(commit,source=line,line_number=number,text=text[:1200]))
    save(record/'script_history.json',history);save(record/'script_evidence.json',excerpts)


def group_fsfs(fsl):
    result=set()
    if not fsl.is_dir():return []
    for top in fsl.glob('L3*'):
        if not top.is_dir():continue
        dirs=[top]+[p for p in top.iterdir() if p.is_dir() and not p.name.endswith('.gfeat')]
        for parent in dirs:
            for group in parent.glob('*.gfeat'):
                text=str(group)
                if re.search(r'(?:[-_])n132(?:[-_.]|$)',text) and (group/'design.fsf').is_file():result.add(group/'design.fsf')
    return sorted(result)


def extraction_commands(evidence):
    result=[]
    for row in evidence:
        line=row.get('text','')
        if 'fslmeants' not in line:continue
        try:
            tokens=shlex.split(line[line.index('fslmeants'):]);options={}
            for k,value in enumerate(tokens):
                if value in ('-i','-m','-o') and k+1<len(tokens):options[value]=tokens[k+1]
            if '-m' not in options:continue
            result.append(dict(row,BOLD=options.get('-i'),mask=options['-m'],timecourse=options.get('-o'),unresolved_variables=any('$' in v or v in ('MASK','DATA','PHYS') for v in options.values()),
                evidence_scope='Recorded text; execution and original file content require corroboration.'))
        except ValueError:continue
    return result


def mask_candidates(archive,record,extra_paths=()):
    paths=[]
    for folder in (ROOT/'masks',ROOT/'masks/archive'):
        if folder.is_dir():paths.extend(p for p in folder.iterdir() if p.is_file() and p.name.endswith(('.nii','.nii.gz')) and (SEEDLIKE.search(p.name) or 'intersection' in p.name.lower()))
    paths.extend(Path(p) for p in extra_paths)
    historical=ROOT/'code/insula_ica/records/2026-10-30-seed-provenance'
    if historical.is_dir():
        archive.capture(historical/'historical-mask-inventory.json','historical_mask_provenance')
        paths.extend((historical/'primary_maps').glob('*.nii.gz'))
    candidates={};errors=[]
    for path in sorted(set(paths)):
        try:
            im=nib.load(path);core.require(len(im.shape)==3,'Candidate is not 3D')
            data=im.get_fdata();core.require(np.isfinite(data).all(),'Nonfinite mask');support=data>.5
            core.require(support.any(),'Empty mask');fingerprint=sha(np.array(im.shape,dtype='<i8').tobytes()+np.packbits(support).tobytes()+np.round(im.affine,5).astype('<f8').tobytes())
            raw=archive.capture(path,'candidate_mask_3D')
            if fingerprint not in candidates:candidates[fingerprint]=dict(id=fingerprint,aliases=[],image=im,indices=np.argwhere(support),nvox=int(support.sum()),binary_or_integer_labels=bool(np.all(data==np.round(data))))
            candidates[fingerprint]['aliases'].append(dict(path=str(path),sha256=raw.get('sha256')))
        except Exception as exc:errors.append(dict(path=str(path),error=str(exc)))
    save(record/'candidate_masks.json',[{k:v for k,v in c.items() if k not in ('image','indices')} for c in candidates.values()]);save(record/'candidate_mask_errors.json',errors)
    return list(candidates.values())


def exact_indices(mask,reference):
    transform=np.linalg.inv(reference.affine)@mask['image'].affine;rotation=transform[:3,:3];rounded=np.round(rotation)
    core.require(np.allclose(rotation,rounded,atol=1e-5) and np.all(np.abs(rounded).sum(0)==1) and np.all(np.abs(rounded).sum(1)==1),'Different voxel lattices; no interpolation permitted')
    ijk=nib.affines.apply_affine(transform,mask['indices']);integer=np.round(ijk).astype(int)
    core.require(np.allclose(ijk,integer,atol=1e-4) and np.all(integer>=0) and np.all(integer<np.array(reference.shape[:3])),'Mask not exactly embedded in original BOLD grid')
    return tuple(integer.T)


def raw_match(observed,predicted):
    result=affine_check(observed,predicted)
    if len(observed)!=len(predicted) or not np.isfinite(predicted).all():return dict(result,raw_match=False)
    delta=np.asarray(observed)-predicted;sd=float(np.std(observed,ddof=1));centered=delta-delta.mean()
    error=float(np.sqrt(np.mean(centered**2))/max(sd,1e-12))
    return dict(result,raw_max_absolute_difference=float(np.max(np.abs(delta))),centered_relative_RMS=error,
        raw_match=bool(np.allclose(observed,predicted,rtol=1e-7,atol=1e-5) and error<1e-4))


def reproduce(phys,candidates,series,record,archive,max_bolds):
    # Each selected original 4D image is read once for all candidate masks and
    # every saved physiological series referencing it. No new BOLD substitute.
    grouped=defaultdict(list)
    for row in phys:
        if row.get('series_key') in series and row.get('original_BOLD_exists'):
            grouped[row['original_BOLD_as_saved']].append(row)
    selected=[];families=set()
    for path,rows in sorted(grouped.items()):
        new={(r['task'],r['folder_seed_label']) for r in rows}-families
        if new:selected.append((path,rows));families.update(new)
        if len(selected)>=max_bolds:break
    results=[]
    for number,(path,rows) in enumerate(selected,1):
        try:
            image=nib.load(path);core.require(len(image.shape)==4,'Original input is not 4D')
            required=int(np.prod(image.shape))*4*2+4*1024**3
            if Path('/proc/meminfo').is_file():
                available=int(re.search(r'MemAvailable:\s+(\d+)',Path('/proc/meminfo').read_text())[1])*1024
                core.require(available>required,'Insufficient available RAM for serial original-BOLD check')
            before=core.identity(path);data=np.asarray(image.dataobj,dtype='float32');predictions={}
            for candidate in candidates:
                try:
                    core.require(candidate['binary_or_integer_labels'],'Noninteger mask: extraction semantics require separate review')
                    idx=exact_indices(candidate,image);predicted=data[idx].mean(axis=0,dtype=np.float64)
                    destination=record/'numerical_checks';destination.mkdir(exist_ok=True)
                    np.savetxt(destination/f'original-BOLD-{number}-{candidate["id"][:12]}.tsv',predicted,fmt='%.12g')
                    predictions[candidate['id']]=predicted
                except Exception as exc:results.append(dict(BOLD=path,candidate_id=candidate['id'],status='candidate_not_compared',error=str(exc)))
            del data
            core.require(core.identity(path,False)==dict(before,sha256=None),'Original BOLD changed during read')
            for row in rows:
                for candidate_id,predicted in predictions.items():
                    value=raw_match(series[row['series_key']],predicted)
                    result=dict(case=row['case'],task=row['task'],folder_seed_label=row['folder_seed_label'],physiological_file=row['series_file_used_for_check'],external_file_as_saved=row['physiological_file_as_saved'],series_source=row['series_source'],candidate_id=candidate_id,BOLD_identity=before,**value)
                    result['supports_saved_model_and_current_original_BOLD']=bool(value['raw_match'] and supported_model_match(row['model_check']))
                    results.append(result)
            print(f'Original BOLD comparison {number}/{len(selected)} completed',flush=True)
        except Exception as exc:results.append(dict(BOLD=path,status='original_BOLD_not_compared',error=str(exc)))
    save(record/'candidate_timecourse_comparisons.json',results)
    return dict(original_BOLD_images_present=len(grouped),original_BOLD_images_selected=len(selected),comparisons=len(results),raw_candidate_matches=sum(r.get('raw_match',False) for r in results),model_supported_matches=sum(r.get('supports_saved_model_and_current_original_BOLD',False) for r in results))


def run(record,args,state):
    registry=json.loads(REGISTRY.read_text());subjects={str(r['subject']).removeprefix('sub-') for r in registry};accepted={r['key'] for r in registry}
    core.require(len(subjects)==132 and len(accepted)==490,'Unexpected manuscript cohort')
    archive=Archive(record);phys=[];fsfs=[];links=[];missing=[];evidence=[];series={};queue=deque();seen=set();chain=set()
    try:
        archive.capture(REGISTRY,'accepted_run_registry');script_history(archive,record)
        fsl=ROOT/'derivatives/fsl'
        for source in group_fsfs(fsl):
            settings=core.fsl_settings(source);refs=[v for k,v in settings.items() if k.startswith('feat_files(')]
            relevant=bool(SEEDLIKE.search(str(source)) or any(SEEDLIKE.search(v) for v in refs))
            if relevant:queue.append((source,'L3_group',True))
            fsfs.append(dict(path=str(source),level='L3_inventory',relevant=relevant,lower_inputs=len(refs)))
        for subject in sorted(subjects):
            session=fsl/f'sub-{subject}/ses-01'
            if not session.is_dir():continue
            for directory in sorted(session.iterdir()):
                if directory.is_dir() and SEEDLIKE.search(directory.name) and directory.name.endswith(('.feat','.gfeat')):
                    queue.append((directory/'design.fsf','L1' if directory.name.startswith('L1') else 'L2',False))
            for path in session.iterdir():
                if path.is_file() and path.suffix=='.txt' and SEEDLIKE.search(path.name):archive.capture(path,'additional_saved_physiological_candidate')
        state.update(state='tracing_original_PPI_designs');save(record/'state.json',state)
        while queue:
            source,level,linked=queue.popleft();source,relocation=locate(str(source));key=str(source.resolve())
            if linked:chain.add(key)
            if key in seen:continue
            seen.add(key)
            if not source.is_file():missing.append(dict(path=str(source),level=level,reason='Saved design absent'));continue
            try:
                archive.capture(source,'saved_PPI_design_FSF');settings=core.fsl_settings(source);refs={k:v for k,v in settings.items() if k.startswith('feat_files(')}
                row=dict(path=str(source),level=level,group_linked=linked,path_resolution=relocation,lower_inputs=len(refs),npts=settings.get('fmri(npts)'),TR=settings.get('fmri(tr)'),contrast_names={k:v for k,v in settings.items() if k.startswith('fmri(conname_')})
                fsfs.append(row)
                for name in ('design.mat','design.con','design.grp'):archive.capture(source.parent/name,'saved_PPI_model')
                if level!='L1':
                    for name,value in refs.items():
                        located,how=locate(value);directory=located if located.name.endswith(('.feat','.gfeat')) else located.parent.parent
                        if directory.parent.name.endswith('.gfeat'):directory=directory.parent
                        child=directory/'design.fsf';child_level='L1' if directory.name.startswith('L1') else 'L2'
                        links.append(dict(source_FSF=str(source),setting=name,input_as_saved=value,located_input=str(located),path_resolution=how,child_FSF=str(child)))
                        if '/derivatives/fsl/' in str(directory) and re.search(r'/sub-(\d+)/',str(directory)) and re.search(r'/sub-(\d+)/',str(directory))[1] in subjects:
                            if linked:queue.appendleft((child,child_level,linked))
                            else:queue.append((child,child_level,linked))
                    continue
                match=re.search(r'L1_task-(trust|ugr)_model-(.+?)_type-(.+)_run-(\d+)_sm-',source.parent.name)
                subject=re.search(r'/sub-(\d+)/',str(source))
                if not match or not subject:missing.append(dict(path=str(source),reason='Unrecognized original PPI directory naming'));continue
                task,model,directory_type,run=match.groups();label=directory_type.split('_seed-',1)[1] if '_seed-' in directory_type else 'unlabeled_'+directory_type
                run_key=f'sub-{subject[1]}_ses-01_task-{task}_run-{run}'
                design=core.fsl_matrix(source.parent/'design.mat') if (source.parent/'design.mat').is_file() else None
                for k,title in settings.items():
                    ev=re.fullmatch(r'fmri\(evtitle(\d+)\)',k)
                    if not ev:continue
                    j=int(ev[1]);value=settings.get(f'fmri(custom{j})','')
                    named=bool(re.match(r'phys|seed.*(?:time|course)',title,re.I))
                    heuristic=settings.get(f'fmri(shape{j})')=='2' and bool(SEEDLIKE.search(Path(value).name))
                    if not named and not heuristic:continue
                    tc,how=locate(value)
                    item=dict(case=str(source.parent),run_key=run_key,subject=subject[1],task=task,run=int(run),model=model,folder_seed_label=label,accepted_manuscript_run=run_key in accepted,phys_EV=j,phys_title=title,physiological_file_as_saved=value,physiological_file=str(tc),path_resolution=how,original_BOLD_as_saved=settings.get('feat_files(1)',''),original_BOLD_exists=Path(settings.get('feat_files(1)','')).is_file(),model_check=dict(status='unavailable'))
                    captured=archive.capture(tc,'saved_physiological_timecourse');item['timecourse_source']=captured
                    cached=source.parent/'custom_timing_files'/f'ev{j}.txt'
                    item['FEAT_cached_EV_path']=str(cached);item['cached_timecourse_source']=archive.capture(cached,'FEAT_cached_physiological_EV')
                    available={}
                    for kind,path in (('external',tc),('FEAT_cached_EV',cached)):
                        try:
                            x=np.loadtxt(path,ndmin=2);core.require(x.shape[0]>2 and x.shape[1]==1 and np.isfinite(x).all(),'Phys input is not a finite per-volume single column')
                            available[kind]=x[:,0]
                            item[kind+'_model_check']=check_model(settings,j,x[:,0],design) if design is not None else dict(status='model_unavailable')
                        except Exception as exc:item[kind+'_timecourse_error']=str(exc)
                    try:
                        kind='FEAT_cached_EV' if 'FEAT_cached_EV' in available else 'external';x=available[kind]
                        item.update(series_source=kind,series_file_used_for_check=str(cached if kind=='FEAT_cached_EV' else tc),nvalues=len(x),mean=float(x.mean()),SD=float(x.std(ddof=1)),series_key=sha(x.astype('<f8').tobytes()))
                        series[item['series_key']]=x
                        item['model_check']=item[kind+'_model_check']
                        if set(available)=={'external','FEAT_cached_EV'}:item['external_vs_FEAT_cache']=raw_match(available['external'],available['FEAT_cached_EV'])
                    except Exception as exc:item['timecourse_error']=str(exc)
                    phys.append(item)
                logs=[source.parent/'report.log',source.parent/'report_log.html']
                logs.extend(sorted((source.parent/'logs').glob('feat*'))[:8])
                for path in logs:
                    if path.is_file():source_excerpt(path,evidence)
                if len(phys)%50==0:
                    state.update(physiological_EVs=len(phys),saved_designs_read=len(seen));save(record/'state.json',state)
            except Exception as exc:
                missing.append(dict(path=str(source),level=level,reason=str(exc)))
                (record/'errors').mkdir(exist_ok=True)
                (record/'errors'/(sha(str(source).encode())[:16]+'.txt')).write_text(traceback.format_exc())
        # Wrapper logs are bounded to the project's own known log directories.
        for folder in (ROOT/'logs',ROOT/'code'):
            if folder.is_dir():
                eligible=[p for p in folder.iterdir() if p.is_file() and p.suffix in ('.log','.out','.txt') and re.search(r'ppi|AIns|L1stats|L1.*(?:ugr|trust)',p.name,re.I)]
                for path in sorted(eligible)[:512]:source_excerpt(path,evidence)
        for row in fsfs:
            if row.get('level')!='L3_inventory':row['referenced_by_discovered_group_chain']=str(Path(row['path']).resolve()) in chain
        for row in phys:row['referenced_by_discovered_group_chain']=str((Path(row['case'])/'design.fsf').resolve()) in chain
        commands=extraction_commands(evidence+json.loads((record/'script_evidence.json').read_text()));save(record/'extraction_commands.json',commands)
        candidates=mask_candidates(archive,record,[locate(c['mask'])[0] for c in commands if not c['unresolved_variables']])
        save(record/'saved_designs.json',fsfs);save(record/'analysis_chain.json',links);save(record/'physiological_EVs.json',phys);save(record/'missing_sources.json',missing);save(record/'log_evidence.json',evidence)
        fields=('run_key','task','model','folder_seed_label','accepted_manuscript_run','referenced_by_discovered_group_chain','phys_EV','physiological_file_as_saved','physiological_file','FEAT_cached_EV_path','series_source','series_file_used_for_check','path_resolution','nvalues','mean','SD','original_BOLD_as_saved','original_BOLD_exists')
        table(record/'physiological_index.tsv',[dict({k:r.get(k,'') for k in fields},model_column_status=r['model_check']['status']) for r in phys])
        state.update(state='checking_available_original_BOLD_against_candidates');save(record/'state.json',state)
        reproduction=reproduce([r for r in phys if r['accepted_manuscript_run']],candidates,series,record,archive,args.max_bolds) if args.max_bolds else dict(skipped=True)
        duplicate=defaultdict(list)
        for r in phys:
            if 'series_key' in r:duplicate[r['run_key'],r['series_key']].append(dict(case=r['case'],folder_seed_label=r['folder_seed_label'],file=r['physiological_file']))
        save(record/'identical_series_across_seed_labels.json',[dict(run_key=k[0],series_key=k[1],uses=v) for k,v in duplicate.items() if len({r['folder_seed_label'] for r in v})>1])
        state.update(state='PPI_seed_provenance_audit_complete_review_required',subjects=132,accepted_runs=490,physiological_EVs=len(phys),accepted_run_physiological_EVs=sum(r['accepted_manuscript_run'] for r in phys),FEAT_cached_physiological_EVs=sum(r.get('series_source')=='FEAT_cached_EV' for r in phys),external_vs_cache_mismatches=sum(not r['external_vs_FEAT_cache']['raw_match'] for r in phys if 'external_vs_FEAT_cache' in r),model_column_status=dict(Counter(r['model_check']['status'] for r in phys)),case_counts=dict(Counter(r['task']+'/'+r['folder_seed_label'] for r in phys)),group_linked_seed_counts=dict(Counter(r['task']+'/'+r['folder_seed_label'] for r in phys if r['referenced_by_discovered_group_chain'])),distinct_candidate_supports=len(candidates),missing_sources=len(missing),numerical_reproduction=reproduction,production_modified=False)
        save(record/'state.json',state)
    finally:
        archive.close();save(record/'source_manifest.json',archive.files)
    (record/'report.md').write_text('# Original PPI seed provenance audit\n\n'
        'No analysis was rerun or modified. Original L3 → L2 → L1 links, saved physiological custom EV paths, full saved FSFs/model matrices and small timecourse files are indexed here. '
        'The scan covers original project N132 group outputs and original insula-like PPI directories for manuscript participants; accepted/excluded runs remain explicitly marked. '
        'Current and Git-historical scripts plus source log excerpts are evidence, not proof that a filename was used in a particular run.\n\n'
        'physiological_EVs.json checks both FEAT custom_timing_files/evN.txt snapshots and current external series against the actual saved design column when direct unfiltered physiological-EV processing supports an affine comparison. '
        'Cached snapshots take priority for numerical reconstruction, while external-versus-cache differences remain recorded. '
        'Mismatch can indicate an overwritten source; unsupported processing is left unresolved. '
        'A numerical candidate match is strongest when both current original-BOLD mean extraction and saved-design membership agree. '
        'Only exact FSF-referenced BOLD paths are used, never regenerated substitutes. Candidate masks must share the exact voxel lattice; field-of-view embedding preserves voxel centers without interpolation. '
        'Binary/integer mask means are computed with float64 accumulation in a single serial 4D read per selected original image. '
        'The sample is selected by task/seed family before numerical matching, with at most '+str(args.max_bolds)+' original BOLD images. '
        'Missing BOLD, equivalent mask supports, stale external series, or lack of contemporaneous extraction logs may prevent an exact historical seed-file attribution.\n\n'
        'source_manifest.json maps original paths/Git commits to exact content-addressed members of split source-files-*.tar.gz archives. '
        'Archives contain small text/model sources and 3D mask candidates, never participant 4D images. '
        'candidate_timecourse_comparisons.json records raw and affine agreements separately; raw numerical agreement alone does not prove historical filename identity. '
        'Use analysis_chain.json to identify which L1 seed variants actually feed saved group analyses, rather than assuming every directory was used.\n')
    print(json.dumps(state,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--publish',action='store_true');p.add_argument('--max-bolds',type=int,default=4);p.add_argument('--brief',action='store_true');args=p.parse_args()
    RECORDS.mkdir(exist_ok=True)
    if args.brief:
        choices=sorted(RECORDS.glob('linux2-PPI-seed-audit-*'))
        if not choices:print('No PPI seed audit record.');return
        record=choices[-1];state=json.loads((record/'state.json').read_text());print('Record:',record.name);print('State:',state['state']);print('Saved designs:',state.get('saved_designs_read','see index'));print('Physiological EVs:',state.get('physiological_EVs',0));print('Model checks:',state.get('model_column_status',{}));print('Original BOLD checks:',state.get('numerical_reproduction','pending'));print('Error:',state.get('error','none'))
        if (record/'publication-error.txt').exists():print('Publication:',(record/'publication-error.txt').read_text())
        return
    core.require(0<=args.max_bolds<=8,'Use 0-8 serial original-BOLD checks')
    (ROOT/'code/ppi_seed_audit/output').mkdir(exist_ok=True)
    with (ROOT/'code/ppi_seed_audit/output/audit.lock').open('a') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('PPI seed audit already running.');return
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6];record=RECORDS/('linux2-PPI-seed-audit-'+stamp);record.mkdir()
        state=dict(state='collecting_scripts_and_seed_candidates',max_original_BOLD_reads=args.max_bolds);save(record/'state.json',state)
        try:run(record,args,state)
        except Exception as exc:
            state.update(state='PPI_seed_audit_failed',error=str(exc));save(record/'state.json',state);(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        if args.publish:
            try:publish(record,'Audit original PPI seed provenance, saved designs and physiological timecourses')
            except Exception as exc:(record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Record:',record,flush=True)
        raise SystemExit(0 if state['state']=='PPI_seed_provenance_audit_complete_review_required' else 1)


if __name__=='__main__':main()

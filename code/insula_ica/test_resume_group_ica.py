from contextlib import ExitStack, nullcontext
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import run_group_ica as g
import resume_group_ica as r
import group_ica_vn_pilot as vn


class ResumeTests(unittest.TestCase):
    def test_original_commands_remain_identical_to_the_real_linux2_plan(self):
        record=Path(__file__).parent/'records/linux2-paired-group-ica-20261028T215022Z-679a3e'
        root=Path(json.loads((record/'state.json').read_text())['output_root'])
        for row in json.loads((record/'ica-plan.json').read_text()):
            self.assertEqual(g.melodic_command(row['mode'],row['roi'],row['dimensionality'],root),row['argv'])

    def test_retry_changes_iteration_budget_and_start_without_changing_concat_or_tolerance(self):
        args=g.melodic_command('joint_bold','right',5,Path('/analysis'),output_folder=Path('/analysis/retry'),seed=r.SEEDS[1],maxit=r.MAXIT)
        self.assertEqual(args[args.index('-a')+1],'concat')
        self.assertIn('--eps=0.00005',args);self.assertIn('--nl=pow3',args)
        self.assertIn('--maxit=10000',args);self.assertNotIn('--maxrestart=10',args)
        self.assertIn('--vn',args)
        native=g.melodic_command('joint_bold','right',5,Path('/analysis'),variance_normalization=True)
        self.assertNotIn('--vn',native);self.assertIn('--update_mask',native)

    def test_explicit_restart_recovers_nonconvergence_and_retains_failed_files(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);record=root/'record';(record/'attempts').mkdir(parents=True)
            seeds=[];failed=[]
            def fit(mode,roi,dim,root,masks,slices,reference,**kw):
                folder=kw['output_folder'];folder.mkdir(parents=True);seeds.append(kw['seed'])
                if len(seeds)==1:
                    log=folder.parent/(folder.name+'.log');log.write_text('No convergence after 10000 steps\n')
                    failed.append((log,log.read_bytes()))
                    raise ValueError('MELODIC produced no maps')
                return dict(mode=mode,roi=roi,dimensionality=dim,oriented_maps=str(folder/'oriented_IC.nii.gz'))
            with patch.object(g,'fit_one',side_effect=fit):
                result=r.retry_fit(('joint_bold','right',5),root,None,None,None,record)
            self.assertEqual(result['state'],'complete');self.assertEqual(seeds,list(r.SEEDS[:2]))
            self.assertEqual(failed[0][0].read_bytes(),failed[0][1])
            self.assertEqual(len(list((record/'attempts').glob('*.json'))),2)

    def test_validation_failure_is_not_retried_as_an_optimizer_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);record=root/'record';(record/'attempts').mkdir(parents=True)
            with patch.object(g,'fit_one',side_effect=ValueError('MELODIC changed fixed mask')) as fit:
                with self.assertRaisesRegex(ValueError,'changed fixed mask'):
                    r.retry_fit(('joint_bold','right',5),root,None,None,None,record)
                self.assertEqual(fit.call_count,1)

    def test_unknown_partial_retry_is_preserved_without_restarting_it(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder=root/'ica_retries/joint_bold/right/dim-5'/f'seed-{r.SEEDS[0]}_maxit-{r.MAXIT}'
            folder.mkdir(parents=True);(folder/'unknown.txt').write_text('inspect me')
            with patch.object(g,'fit_one') as fit:
                with self.assertRaisesRegex(ValueError,'Unverified partial retry'):
                    r.retry_fit(('joint_bold','right',5),root,None,None,None,root)
                fit.assert_not_called()
            self.assertEqual((folder/'unknown.txt').read_text(),'inspect me')

    def test_resume_runs_only_eleven_unfinished_fits_without_preprocessing_or_input_generation(self):
        with tempfile.TemporaryDirectory() as d,ExitStack() as p:
            root=Path(d).resolve();source=root/'source';source.mkdir();record=root/'record';record.mkdir()
            for name in ('attempts','figures','errors'):(record/name).mkdir()
            (root/'masks').mkdir();reference=nib.Nifti1Image(np.ones((5,5,5),dtype='float32'),np.eye(4))
            masks={roi:np.ones(reference.shape,bool) for roi in g.core.ROIS}
            for roi in g.core.ROIS:g.core.save_image(root/'masks'/(roi+'.nii.gz'),masks[roi].astype('uint8'),reference)
            verified={(mode,roi):'verified' for mode in g.inputs.MODES for roi in g.core.ROIS}
            failed={('joint_bold','left',7),('joint_bold','right',5),('joint_bold','right',7),('native_film','right',5),('native_film','right',7)}
            originals=[]
            def result(mode,roi,dim,folder,seed=g.SEED,maxit=2000,custom=False,variance_normalization=False):
                folder.mkdir(parents=True);actual=dim or 4
                data=np.random.default_rng(actual).normal(size=(*reference.shape,actual)).astype('float32')
                maps=folder/'oriented_IC.nii.gz';g.core.save_image(maps,data,reference)
                kwargs=dict(output_folder=folder,seed=seed,maxit=maxit,variance_normalization=variance_normalization) if custom else {}
                value=dict(mode=mode,roi=roi,dimensionality=dim,actual_components=actual,oriented_maps=str(maps),
                    seed=seed,max_iterations=maxit,argv=g.melodic_command(mode,roi,dim,root,**kwargs),
                    mask_sha256=g.core.sha(root/'masks'/(roi+'.nii.gz')),verified_inputs_sha256='verified',map_identity=g.core.identity(maps))
                g.save(folder/'complete.json',value);return value
            for mode,roi,dim in g.plans():
                if dim and (mode,roi,dim) not in failed:
                    value=result(mode,roi,dim,root/'ica'/mode/roi/f'dim-{dim}')
                    originals.append((Path(value['oriented_maps']),value['map_identity']))
            prepared=dict(runs={f'run{k}':dict(state='complete') for k in range(490)})
            p.enter_context(patch.object(r,'verify_prepared',return_value=(prepared,source,verified,masks,None,reference)))
            p.enter_context(patch.object(g.core,'software',return_value=dict(fsl='fixture')))
            memory=p.enter_context(patch.object(g.inputs.MemoryGate,'reserve',side_effect=lambda amount:nullcontext()))
            def fit(mode,roi,dim,root,masks,slices,ref,**kwargs):
                return result(mode,roi,dim,kwargs['output_folder'],kwargs['seed'],kwargs['maxit'],True,kwargs['variance_normalization'])
            mock=p.enter_context(patch.object(g,'fit_one',side_effect=fit))
            p.enter_context(patch.object(g,'sessions',side_effect=AssertionError('4D files must not be rebuilt')))
            p.enter_context(patch.object(g.inputs,'ensure_native',side_effect=AssertionError('FILM must not run')))
            def review(results,record,masks,ref,state):
                self.assertEqual(len(results),24);state['state']='paired_ICA_complete_review_pending';(record/'README.md').write_text('Review')
            p.enter_context(patch.object(g,'review_outputs',side_effect=review))
            def pilot(root,record,state,*args):state['native_variance_normalization_pilot']={m:dict(state='complete') for m in g.inputs.MODES}
            p.enter_context(patch.object(vn,'run_pilot',side_effect=pilot))
            state=dict(runs={},fits={});r.run(source,root,record,state,11)
            self.assertEqual(mock.call_count,11)
            self.assertEqual(sum(c.args[2]==0 for c in mock.call_args_list),6)
            self.assertEqual(len(state['fits']),24);self.assertEqual(state['state'],'paired_ICA_complete_review_pending')
            self.assertEqual({c.args[0] for c in memory.call_args_list},{2,16})
            for path,identity in originals:self.assertEqual(g.core.identity(path),identity)
            chosen=json.loads((root/'selected_ica_fits.json').read_text())['fits']
            self.assertEqual(len(chosen),24)
            for plan in failed:self.assertIn('/ica_retries/',chosen[f'{plan[0]}/{plan[1]}/dim-{plan[2]}']['oriented_maps'])


if __name__=='__main__':unittest.main()

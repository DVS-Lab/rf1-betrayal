"""Numerical safeguards and an optional real-FSL synthetic integration pilot.

python -m unittest discover -s code/insula_ica -p 'test_*.py' -v
python code/insula_ica/test_insula_ica.py --integration /absolute/scratch/path
"""
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
import numpy as np
import nibabel as nib
import insula_ica as p


class NumericalTests(unittest.TestCase):
    def test_joint_projection_preserves_orthogonality(self):
        rng = np.random.default_rng(11)
        n = rng.normal(size=(130, 4))
        n[:, 1] += n[:, 0]
        q = p.basis(np.column_stack([np.ones(130), n, n[:, 0]]))
        task = n[:, :2] + rng.normal(size=(130, 2))
        truth = np.array([[2.1], [-.8]])
        y = task @ truth + n @ rng.normal(size=(4, 1))
        yr, xr = p.residualize(q, y), p.residualize(q, task)
        fit = np.linalg.lstsq(xr, yr, rcond=None)[0]
        np.testing.assert_allclose(fit, truth, atol=1e-12)
        np.testing.assert_allclose(q.T @ (yr-xr@fit), 0, atol=1e-12)
        # Regressing an unprojected task design gives a different answer.
        self.assertGreater(np.linalg.norm(np.linalg.lstsq(task, yr, rcond=None)[0]-truth), .1)

    def test_partial_connectivity_matches_simultaneous_model(self):
        rng = np.random.default_rng(12)
        n = rng.normal(size=(150, 5))
        t = rng.normal(size=(150, 5)) + n
        e = rng.normal(size=(150, 3))
        y = t @ np.eye(5) + n @ np.ones((5, 5)) + e @ np.ones((3, 5))
        q = p.basis(np.column_stack([np.ones(150), n, e]))
        actual = np.linalg.pinv(p.residualize(q, t)) @ p.residualize(q, y)
        np.testing.assert_allclose(actual, np.eye(5), atol=1e-12)

    def test_contrast_estimability(self):
        x = np.column_stack([np.arange(120), np.arange(120)])
        with self.assertRaisesRegex(ValueError, 'not estimable'):
            p.contrast_fit(x, np.ones((120, 5)), np.array([1, -1]), 1, 30, 1e6)

    def test_spatial_sign_does_not_use_task_effect(self):
        data = np.zeros((2, 2, 2, 2))
        data[0, 0, 0] = [-5, 4]
        data[1, 1, 1] = [2, -3]
        oriented, signs = p.orient_maps(data.copy(), np.ones((2, 2, 2), bool))
        self.assertEqual(signs, [-1, 1])
        self.assertEqual(oriented[0, 0, 0, 0], 5)

    def test_thirty_test_family(self):
        np.testing.assert_allclose(p.bonferroni([.001, .1]), [.03, 1])

    def test_hemispheres_follow_mni_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            # Reversed voxel x-axis: left must come from world coordinates.
            a = np.diag([-2., 2., 2., 1.]); a[0, 3] = 8
            atlas = np.full((10, 4, 4), 2, dtype='uint8')
            path = folder/'atlas.nii.gz'; nib.save(nib.Nifti1Image(atlas, a), path)
            xml = folder/'atlas.xml'; xml.write_text('<atlas><data><label index="1">Insular Cortex</label></data></atlas>')
            masks, _, _ = p.anatomical_masks(path, xml, nib.load(path))
            self.assertTrue(masks['left'][5:].all())
            self.assertFalse(masks['left'][:5].any())
            self.assertFalse((masks['left'] & masks['right']).any())


def integration(folder):
    """Exercise real feat_model + MELODIC + both regressions on fake NIfTIs."""
    require = p.require
    require(not folder.exists(), 'Integration destination must be new')
    original_root = p.ROOT
    folder.mkdir(parents=True)
    p.ROOT = folder
    (folder/'code').mkdir()
    (folder/'code/sublist_n132.txt').write_text('10001\n10002\n10003\n')
    rng = np.random.default_rng(2026)
    affine = np.diag([4., 4., 4., 1.]); affine[:3, 3] = [-48, -40, -24]
    shape, nvol, tr = (25, 24, 19), 180, 1.615
    upstream = folder/'upstream'
    # A lightweight local Git checkout is sufficient for source-hash provenance.
    import subprocess
    subprocess.run(['git', 'init', '-q', str(folder)], check=True)
    subprocess.run(['git', '-C', str(folder), '-c', 'user.name=Synthetic', '-c', 'user.email=test@example.invalid', 'commit', '--allow-empty', '-qm', 'Synthetic fixture'], check=True)
    subprocess.run(['git', 'init', '-q', str(upstream)], check=True)
    subprocess.run(['git', '-C', str(upstream), '-c', 'user.name=Synthetic', '-c', 'user.email=test@example.invalid', 'commit', '--allow-empty', '-qm', 'Synthetic fixture'], check=True)
    records = []
    shared_spatial = np.zeros((*shape, 9), dtype='float32')
    labels = rng.choice(9, size=shape, p=np.arange(9, 0, -1)/45)
    for k in range(9):
        shared_spatial[..., k] = (labels == k) * rng.choice([-1., 1.], size=shape)
    for sub in ('10001', '10002', '10003'):
        qr = {'sub_id': sub, 'include_final': 'TRUE'}
        for task, (model, _, _) in p.TASKS.items():
            qr[f'{task}_run1_full_ok'] = 'TRUE'; qr[f'{task}_run2_full_ok'] = 'FALSE'
            key = f'sub-{sub}_ses-01_task-{task}_run-1'
            func = upstream/f'derivatives/fmriprep-24/sub-{sub}/ses-01/func'; func.mkdir(parents=True, exist_ok=True)
            base = func/f'{key}_part-mag'
            bold = Path(str(base)+'_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz')
            mask = Path(str(base)+'_space-MNI152NLin6Asym_desc-brain_mask.nii.gz')
            nuisance = rng.normal(size=(nvol, 13))
            cols = p.BASE_CONFOUNDS + ['non_steady_state_outlier00']
            nss = np.zeros(nvol); nss[0] = 1
            confvalues = np.column_stack([nuisance, nss])
            fconf = Path(str(base)+'_desc-confounds_timeseries.tsv')
            p.write_table(fconf, [{c:float(v) for c,v in zip(cols, r)} for r in confvalues])
            conf = upstream/f'derivatives/fsl/confounds_tedana-24/sub-{sub}/{key}_desc-TedanaPlusConfounds.tsv'
            conf.parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(conf, confvalues, delimiter='\t')
            noise = rng.normal(size=(*shape, nvol)).astype('float32')
            spatial = shared_spatial
            latent = rng.laplace(size=(nvol, 9)).astype('float32')
            data = 1000 + 2 * (spatial @ latent.T) + noise + nuisance[:, 0].astype('float32')
            image = nib.Nifti1Image(data.astype('float32'), affine)
            image.header.set_zooms((4.,4.,4.,tr)); image.header.set_xyzt_units('mm','sec')
            nib.save(image, bold); nib.save(nib.Nifti1Image(np.ones(shape,'uint8'), affine),mask)
            p.write_json(Path(str(bold).removesuffix('.nii.gz')+'.json'), {'RepetitionTime':tr})
            fsf = folder/f'derivatives/fsl/sub-{sub}/ses-01/L1_task-{task}_model-{model}_type-act_run-1_sm-5.feat/design.fsf'
            fsf.parent.mkdir(parents=True, exist_ok=True)
            template = (original_root/f'templates/L1_task-{task}_model-{model}_type-act.fsf').read_text()
            evdir = folder/f'events/{key}';evdir.parent.mkdir(parents=True, exist_ok=True)
            nev = 10 if task=='trust' else 11
            # Distinct random onsets, task events remain in the generated design.
            for i in range(1, nev):
                ev = folder/f'ev_{key}_{i}.txt'
                onsets = np.sort(rng.uniform(2, nvol*tr-20, 12))
                np.savetxt(ev, np.column_stack([onsets, np.full(12,2), rng.uniform(.5,2,12)]))
                template = __import__('re').sub(r'^set fmri\(custom'+str(i)+r'\).*$', f'set fmri(custom{i}) "{ev}"', template, flags=__import__('re').M)
            replacements = {'OUTPUT':str(fsf.parent), 'NVOLUMES':str(nvol),'DATA':str(bold),'CONFOUNDEVS':str(conf),
                            'SMOOTH':'5','SHAPE_EV':'10','MISSED_TRIAL':str(folder/'absent_miss.txt')}
            for k,v in replacements.items():template=template.replace(k,v)
            fsf.write_text(template)
        records.append(qr)
    qc = folder/'code/master_subject_table_full.csv'
    import csv
    with qc.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    c=json.loads((original_root/'code/insula_ica/config.example.json').read_text())
    c.update(upstream_root=str(upstream), master_qc=str(qc), output_root=str(folder/'output'),expected_cohort_size=3)
    c['provenance_review']={'signal_description':'SYNTHETIC BOLD ONLY','evidence':'synthetic integration fixture',
        'tedana_noise_not_already_removed':True,'ugr_pmods_increase_with_fairness':True,'ugr_coding_evidence':'synthetic increasing amplitudes'}
    template_path=folder/'template.nii.gz'
    nib.save(nib.Nifti1Image(rng.normal(size=shape),affine),template_path)
    c['network_templates']=[{'name':'synthetic reference','path':str(template_path),'citation':'SYNTHETIC TEST ONLY'}]
    require(p.audit(c)==0, 'Synthetic audit failed')
    p.prepare(c);p.ica(c,'pilot');p.regress(c,'pilot');p.summarize(c,'pilot')
    key = 'sub-10001_ses-01_task-trust_run-1'
    p.run(['fsl_glm', '-i', folder/'output/prepared'/f'{key}_clean.nii.gz',
           '-d', folder/'output/pilot/left/dim-5/oriented_IC.nii.gz',
           '-m', folder/'output/audit/left_mask.nii.gz', '--demean', '-o', folder/'fsl_stage1.txt'], folder/'crosscheck')
    fsl_tc = np.loadtxt(folder/'fsl_stage1.txt')
    actual_tc = np.loadtxt(folder/'output/pilot/left/runs'/key/'timecourses.tsv')
    fsl_tc /= fsl_tc.std(axis=0, ddof=1)
    np.testing.assert_allclose(fsl_tc, actual_tc, atol=1e-4)
    table=p.read_table(folder/'output/pilot/group_task_effects.tsv')
    require(len(table)==30, 'Wrong multiplicity family size')
    require(all(int(r['n'])==3 for r in table), 'Runs incorrectly treated as participants')
    # Repeat against the complete cohort with a separate ICA basis/output area.
    p.ica(c,'full');p.regress(c,'full');p.summarize(c,'full')
    print('SYNTHETIC INTEGRATION PASSED:',folder/'output/full/results_summary.md')


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--integration':
        integration(Path(sys.argv[2]).resolve())
    else:
        unittest.main()

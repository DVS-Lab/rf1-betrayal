"""Generalized spatial regression, independent ROI means and safe resume."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

import nibabel as nib
import numpy as np
from scipy import stats

import task_dimension_sensitivity as d


class DimensionTests(unittest.TestCase):
    def test_spatial_models_recover_three_and_seven_sources_with_shared_offset(self):
        rng = np.random.default_rng(31)
        for dim in (3, 7):
            maps = rng.normal(size=(70, dim)); maps[:, 1] += .7*maps[:, 0]
            tc = rng.normal(size=(dim, 90)); offset = rng.normal(size=90)+1000
            operator = d.spatial_weights(maps)
            np.testing.assert_allclose(operator@(maps@tc+offset), tc, atol=1e-11)
            np.testing.assert_allclose(operator@np.ones(70), 0, atol=1e-12)
        with self.assertRaisesRegex(ValueError, 'Nonestimable'): d.spatial_weights(np.ones((20, 7)))

    def fixture(self, rng):
        shape = (6, 6, 4); masks = {}
        for roi in d.n.g.core.ROIS:
            mask = np.ones(shape, bool)
            if roi == 'left': mask[3:] = False
            if roi == 'right': mask[:3] = False
            masks[roi] = mask
        bases = {}
        for mode in d.n.g.inputs.MODES:
            for roi in masks:
                for dim in (3, 7):
                    maps = rng.normal(size=(int(masks[roi].sum()), dim))
                    bases[mode, roi, dim] = dict(weights=d.spatial_weights(maps))
        return shape, masks, bases

    def test_mean_is_actual_masked_BOLD_average_once_per_ROI_not_mean_of_ICs(self):
        rng = np.random.default_rng(5); shape, masks, bases = self.fixture(rng)
        data = 1000+rng.normal(size=(*shape, 80)); data[..., 0] += 30
        actual = d.extract_responses(data, bases, masks); layout = d.column_layout(bases, masks)
        self.assertEqual(actual.shape, (80, 63)); self.assertEqual(len(layout), 63)
        mean_columns = [i for i, v in enumerate(layout) if v['series'] == 'insula_mean']
        self.assertEqual(len(mean_columns), 3)
        for i in mean_columns:
            np.testing.assert_allclose(actual[:, i], data[masks[layout[i]['roi']]].mean(0))
        # A spatially uniform response goes to the ROI means and not any IC.
        uniform = np.broadcast_to(np.arange(80), (*shape, 80))
        extracted = d.extract_responses(uniform, bases, masks)
        np.testing.assert_allclose(extracted[:, :60], 0, atol=1e-12)
        np.testing.assert_allclose(extracted[:, 60:], np.repeat(np.arange(80)[:, None], 3, axis=1))

    def test_real_NIfTI_run_fits_all_responses_and_resume_avoids_4D_read(self):
        rng = np.random.default_rng(6); shape, masks, bases = self.fixture(rng); t = 100
        data = (1000+rng.normal(size=(*shape, t))).astype('float32')
        x = rng.normal(size=(t, 13)); cs = np.eye(13)[:3]; labels = list(d.n.BAR_LABELS[2:])
        row = dict(subject='test', task='trust', run=1)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); path = root/'input.nii.gz'; nib.save(nib.Nifti1Image(data, np.eye(4)), path)
            loader = Mock(side_effect=lambda:np.asarray(nib.load(path).dataobj, dtype='float32'))
            signature = dict(policy=d.POLICY, input=d.n.g.core.identity(path))
            folder = root/'run'; layout = d.column_layout(bases, masks)
            rows, cache = d.run_cached_model(row, folder, signature, x, cs, labels, layout, bases, masks, loader)
            self.assertEqual(len(rows), 63*3); self.assertEqual(loader.call_count, 1)
            self.assertFalse(cache['fit_reused'])
            actual = np.load(folder/'fit/standardized_betas.npy')
            zx = (x-x.mean(0))/x.std(0, ddof=1); y = np.load(folder/'stage1_responses.npy')
            zy = (y-y.mean(0))/y.std(0, ddof=1)
            np.testing.assert_allclose(actual, np.linalg.lstsq(zx, zy, rcond=1e-10)[0], atol=1e-12)
            again, cache = d.run_cached_model(row, folder, signature, x, cs, labels, layout, bases, masks, loader)
            self.assertEqual(loader.call_count, 1); self.assertEqual(rows, again); self.assertTrue(cache['fit_reused'])
            (folder/'stage1_responses.npy').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed'):
                d.run_cached_model(row, folder, signature, x, cs, labels, layout, bases, masks, loader)

    def test_interrupted_temporal_fit_reuses_complete_extraction_and_preserves_partial(self):
        rng = np.random.default_rng(7); shape, masks, bases = self.fixture(rng); t = 90
        data = 1000+rng.normal(size=(*shape, t)); x = rng.normal(size=(t, 12)); cs = np.eye(12)[:2]
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp); signature = dict(fixture='unchanged')
            first = Mock(return_value=data)
            d.cached_extract(folder, signature, first, bases, masks)
            (folder/'fit').mkdir(); (folder/'fit/partial.txt').write_text('retain')
            forbidden = Mock(side_effect=AssertionError('4D must not be reread'))
            _, cache = d.run_cached_model(dict(subject='test', task='ugr', run=1), folder, signature,
                x, cs, list(d.n.BAR_LABELS[:2]), d.column_layout(bases, masks), bases, masks, forbidden)
            self.assertTrue(cache['extraction_reused']); self.assertFalse(cache['fit_reused'])
            self.assertEqual(forbidden.call_count, 0)
            self.assertEqual(len(list(folder.glob('fit.interrupted-*/partial.txt'))), 1)

    def test_aggregate_all_dimensions_and_unique_means_reuses_5D_statistics(self):
        rng = np.random.default_rng(8); subjects = [str(k) for k in range(132)]; effects = []; original = []
        selections = [(mode, roi, dim, 'component', k) for mode in d.n.g.inputs.MODES
                      for roi in d.n.g.core.ROIS for dim in (3, 5, 7) for k in range(1, dim+1)]
        selections.extend(('insula_mean', roi, 0, 'insula_mean', 0) for roi in d.n.g.core.ROIS)
        for mode, roi, dim, series, k in selections:
            for label in d.n.BAR_LABELS:
                task = 'ugr' if label.startswith('ugr') else 'trust'; values = rng.normal(size=132)
                for subject, value in zip(subjects, values):
                    for run, offset in ((1, -.2), (2, .2)):
                        effects.append(dict(mode=mode, roi=roi, dimensionality=dim, series=series, component=k,
                            subject=subject, task=task, contrast=label, run=run, effect=float(value+offset)))
                if dim == 5:
                    mean = float(values.mean()); se = float(stats.sem(values)); delta = float(stats.t.ppf(.975, 131)*se)
                    original.append(dict(mode=mode, roi=roi, component=k, task=task, contrast=label,
                        effect=mean, se=se, ci_low=mean-delta, ci_high=mean+delta))
        participants, groups = d.aggregate(effects, subjects)
        self.assertEqual(len(groups), 465); self.assertEqual(len(participants), 465*132)
        self.assertEqual(sum(r['series'] == 'insula_mean' for r in groups), 15)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp); d.n.g.core.write_table(source/'group_network_task_effects.tsv', original)
            verified = d.verify_five_reuse(groups, source)
            self.assertEqual(verified['group_estimates_reused'], 150); self.assertLess(verified['max_absolute_summary_difference'], 1e-12)
        with self.assertRaisesRegex(ValueError, 'Duplicate'): d.aggregate(effects+[effects[0]], subjects)
        with self.assertRaisesRegex(ValueError, 'Incomplete'): d.aggregate(effects[:-2], subjects)

    def test_frozen_sources_cannot_be_silently_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); path = root/'frozen.json'; d.freeze_json(path, dict(x=1)); before = path.read_bytes()
            with self.assertRaisesRegex(ValueError, 'changed'): d.freeze_json(path, dict(x=2))
            self.assertEqual(path.read_bytes(), before)
            source = root/'source'; source.write_bytes(b'first'); target = root/'published'
            d.copy_frozen(source, target); source.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed'): d.copy_frozen(source, target)
            self.assertEqual(target.read_bytes(), b'first')

    def test_new_task_inputs_must_match_completed_five_component_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'complete.json'
            inputs = [dict(path='/model', sha256='original')]
            old = dict(fsf='fsf', FEAT_input='bold', tr=1.615, inputs=inputs)
            path.write_text(json.dumps(dict(signature=old)))
            new = dict(fsf='fsf', FEAT_input='bold', tr=1.615, model_inputs=inputs)
            d.verify_task_source(new, path)
            with self.assertRaisesRegex(ValueError, 'BOLD/FSF/TR'):
                d.verify_task_source(dict(new, FEAT_input='different'), path)
            with self.assertRaisesRegex(ValueError, 'design/nuisance/mask'):
                d.verify_task_source(dict(new, model_inputs=[dict(path='/model', sha256='different')]), path)


if __name__ == '__main__': unittest.main()

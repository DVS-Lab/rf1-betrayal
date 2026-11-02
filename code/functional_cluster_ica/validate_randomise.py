#!/usr/bin/env python3
"""Small native-FSL interface check; synthetic data, not cohort inference."""
import argparse
import os
from pathlib import Path
import tempfile
import nibabel as nib
import numpy as np
import wholebrain as w


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    rng=np.random.default_rng(19)
    with tempfile.TemporaryDirectory(prefix='rf1-dr-randomise-') as tmp:
        root=Path(tmp);shape=(7,7,7);n=16;data=rng.normal(size=(*shape,n)).astype('float32');data[2:5,2:5,2:5]+=1.5
        merged=root/'data.nii.gz';mask=root/'mask.nii.gz';design=root/'ones.mat';contrasts=root/'tails.con'
        nib.save(nib.Nifti1Image(data,np.eye(4)),merged);nib.save(nib.Nifti1Image(np.ones(shape,dtype='uint8'),np.eye(4)),mask)
        w.vest(design,np.ones((n,1)));w.vest(contrasts,[[1],[-1]],True)
        result=w.randomise_job(merged,root/'randomise',mask,design,contrasts,100,3)
        positive=nib.load(root/'randomise/dr_tstat1.nii.gz').get_fdata();negative=nib.load(root/'randomise/dr_tstat2.nii.gz').get_fdata()
        np.testing.assert_allclose(positive,-negative,rtol=1e-5,atol=1e-6)
        reused=w.randomise_job(merged,root/'randomise',mask,design,contrasts,100,3);assert reused['reused']
        args.output.parent.mkdir(parents=True,exist_ok=True)
        w.f.save(args.output,dict(test='Native randomise explicit one-sample design, positive/negative contrasts and 3D TFCE',data='Synthetic 7x7x7 volume, 16 observations, seeded noise plus positive patch; not participant data',permutations=100,production_permutations=5000,fsl_version=(Path(os.environ['FSLDIR'])/'etc/fslversion').read_text().strip(),TFCE_both_tails_finite_probability_maps=True,tstat2_equals_negative_tstat1=True,completed_checkpoint_reused=reused['reused'],command_template=w.randomise_command('132_participant_maps.nii.gz','dr','inference_mask.nii.gz','ones.mat','tails.con',5000)))
        print('Native randomise passed; completed checkpoint reused:',args.output)


if __name__=='__main__':main()

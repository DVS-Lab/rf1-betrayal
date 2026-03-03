# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

ug=$maindir/derivatives/fsl/L3-act/L3_model-3_task-ugr_type-act-n132-cov-ones-flame1/L3_model-3_task-ugr_n132-cov-ones-cope-11_cname-offer-unfairness_pmod-flame1.gfeat/cope1.feat/cluster_mask_zstat2.nii.gz
trust=$maindir/derivatives/fsl/L3-act/L3_model-01_task-trust_type-act-n132-cov-ones-flame1/L3_model-01_task-trust_n132-cov-ones-cope-10_cname-rec-def-flame1.gfeat/cope1.feat/cluster_mask_zstat2.nii.gz

fslmaths $ug -mas $trust "$maindir/masks/ugr-trust_intersection"
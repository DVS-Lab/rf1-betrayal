# Completed unsmoothed ICA exploration

Fixed-fit spatial figures and primary component maps are available now. Automatic ICA continues independently. Whole-brain connectivity and task projections are pending unless state.json reports their completion. No network/task findings are inferred from convergence alone. See policy.json and component_anatomy.tsv.

- [native_film, left, 3 components](figures/native_film_left_dim-3.png)
- [joint_bold, left, 3 components](figures/joint_bold_left_dim-3.png)
- [joint_bold, right, 3 components](figures/joint_bold_right_dim-3.png)
- [native_film, right, 3 components](figures/native_film_right_dim-3.png)
- [native_film, bilateral, 3 components](figures/native_film_bilateral_dim-3.png)
- [joint_bold, bilateral, 3 components](figures/joint_bold_bilateral_dim-3.png)
- [joint_bold, right, 5 components](figures/joint_bold_right_dim-5.png)
- [native_film, left, 5 components](figures/native_film_left_dim-5.png)
- [joint_bold, left, 5 components](figures/joint_bold_left_dim-5.png)
- [native_film, right, 5 components](figures/native_film_right_dim-5.png)
- [native_film, bilateral, 5 components](figures/native_film_bilateral_dim-5.png)
- [joint_bold, bilateral, 5 components](figures/joint_bold_bilateral_dim-5.png)
- [native_film, left, 7 components](figures/native_film_left_dim-7.png)
- [joint_bold, left, 7 components](figures/joint_bold_left_dim-7.png)
- [native_film, right, 7 components](figures/native_film_right_dim-7.png)
- [joint_bold, right, 7 components](figures/joint_bold_right_dim-7.png)
- [native_film, bilateral, 7 components](figures/native_film_bilateral_dim-7.png)
- [joint_bold, bilateral, 7 components](figures/joint_bold_bilateral_dim-7.png)


## Task-projection review ready

All 60 planned projections are available; 8 have Bonferroni-60 p<0.05. Same-cohort exploratory estimates, not independent confirmation. Direction is as stored in each task contrast. Connectivity remains pending at this checkpoint.

- [Joint-BOLD task projections](figures/joint_bold_task_projections.png)
- [Native-FILM task projections](figures/native_film_task_projections.png)
- [All task estimates and CIs](group_task_projections.tsv)


Connectivity and all 60 task projections completed for N=132 / 490 runs. 8/60 task tests have p Bonferroni-60 <0.05. All estimates, CIs and null results are in group_task_projections.tsv. Directions use the stored contrasts; no unfairness-direction sign was inferred. Connectivity images are descriptive means, not thresholded group inference. Atlas regional targets are anatomical summaries, not independent reference-network identification. Native-FILM maps describe innovation-scale coupling; joint-BOLD maps describe original-scale residual coupling. These are same-cohort exploratory associations; no unique betrayal encoding, causal role, independent validation, or formal between-task/hemisphere effect is established.

- [joint_bold, left: anatomy and outside-ROI connectivity](figures/joint_bold_left_anatomy_connectivity.png)
- [joint_bold, right: anatomy and outside-ROI connectivity](figures/joint_bold_right_anatomy_connectivity.png)
- [joint_bold, bilateral: anatomy and outside-ROI connectivity](figures/joint_bold_bilateral_anatomy_connectivity.png)
- [native_film, left: anatomy and outside-ROI connectivity](figures/native_film_left_anatomy_connectivity.png)
- [native_film, right: anatomy and outside-ROI connectivity](figures/native_film_right_anatomy_connectivity.png)
- [native_film, bilateral: anatomy and outside-ROI connectivity](figures/native_film_bilateral_anatomy_connectivity.png)

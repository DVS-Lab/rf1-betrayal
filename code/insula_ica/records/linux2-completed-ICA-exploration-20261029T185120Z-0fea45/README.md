# Completed unsmoothed ICA exploration

Fixed-fit spatial figures and primary component maps are available now. Automatic ICA continues independently. Background connectivity is complete; task effects require the replacement timecourse analysis. No network/task findings are inferred from convergence alone. See policy.json and component_anatomy.tsv.

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


## Task analysis withdrawn (2026-10-29)

The previous task-analysis outputs and conclusions have been removed. Use
`run_network_timecourse_analysis.py` to extract network timecourses from
task-preserved BOLD and fit the task and nuisance model jointly.

All 490 background-connectivity runs are complete. Retained connectivity maps
are descriptive background coefficients, not task-dependent connectivity.

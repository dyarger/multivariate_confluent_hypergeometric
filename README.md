# Multivariate Confluent Hypergeometric covariances

This repository houses code accompanying the paper "On Valid Multivariate Generalizations of the Confluent Hypergeometric Covariance Function."
We break the files into sections below.

Source files: mch_source.R contains files that create matrices of evaluated covariances, as well as log-likelihood functions for different multivariate covariances used in the simulation and data analysis. 
mch.cpp provides relevant C++ functions, including evaluations of confluent hypergeometric functions and computation of the cross-covariances using the discrete Hilbert transform. 
fft_source.R implements computation of cross-covariances through 1d and 2d fast Fourier transforms. 

Plotting/exploration: mch_check_condition.R and the shiny application check_condition summarize the results in the paper when applied to a p=2 example. This includes visualizing how sharp are bounds in the paper, and what the normalized spectral density looks like for the multivariate CH model. 
The file sim_2d.R simulates CH and Matern processes plotted in the paper. 
The files plot_imaginary.R, plot_general_form.R, and plot_2d.R plot spectral cross-covariances, including making plots presented in Section 3.4 of the paper. 

Simulation: mch_simu_all.R runs each of the simulation studies in succession and provides results in tables and figures.

Data analysis: get_soccom_data.R and get_core_data.R extract Argo data for BGC floats and Core Argo floats respectively. 
These two files require downloading versions of the SOCCOM and Argo database. 
The data analysis is then implemented in argo_analysis.R and argo_analysis_cv.R. The first file contains the optimization of parameters and prediction, while the second implements cross-validation. 
The second file should be run twice to reproduce results in the paper, with leave_out_type <- 2 as well as leave_out_type <- 'one'

Data/results: contains selected .RData files from the simulation and data analysis of moderate to low storage size. In principle, results from the paper may be reproduced based on these files, without requiring downloading the Argo and SOCCOM data and running get_soccom_data.R and get_core_data.R. 
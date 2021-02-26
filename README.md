# Simulation study of Local Projections, VARs, and related estimators

Matlab code for large-scale simulation studies of impulse response estimators, including Local Projections (LPs), Vector Autoregressions (VARs), and several variants of these

**Reference:**
[Li, Dake](https://github.com/dake-li), [Mikkel Plagborg-Møller](https://scholar.princeton.edu/mikkelpm), and [Christian K. Wolf](https://www.christiankwolf.com/) (2021), "Local Projections vs. VARs: Lessons From Thousands of DGPs" (:page_facing_up:[paper](Documents/lp_var_simul.pdf), :bar_chart:[supplement](Documents/lp_var_simul_supplement.pdf))

Tested in: Matlab R2020a on Windows 10 PC (64-bit)

**Acknowledgements:**
We rely on penalized LP code by [Regis Barnichon & Christian Brownless](https://github.com/ctbrownlees/MATLAB-package-lproj), as well as VAR model averaging code by [Bruce Hansen](https://www.ssc.wisc.edu/~bhansen/progs/var.html). We have slightly modified both of these sets of codes to improve their run-time without affecting their numerical output.

## Contents

**[Documents](Documents):** Paper, supplement, and documentation
- [lp_var_simul.pdf](Documents/lp_var_simul.pdf): Main paper
- [lp_var_simul_supplement.pdf](Documents/lp_var_simul_supplement.pdf): Online supplement

**[Analytical_Illustration](Analytical_Illustration):** Plots for simple analytical illustration
- [plot_tradeoff.m](Analytical_Illustration/plot_tradeoff.m): Figure 1 in the paper
- [plot_indiff.m](Analytical_Illustration/plot_indiff.m): Figure 2 in the paper

**[Estimation_Routines](Estimation_Routines):** General-purpose impulse response estimation functions
- [BVAR_est.m](Estimation_Routines/BVAR_est.m): Bayesian VAR
- [LP_est.m](Estimation_Routines/LP_est.m): Least-squares LP
- [LP_shrink_est.m](Estimation_Routines/LP_shrink_est.m): Penalized LP ([Barnichon & Brownlees, 2019](https://www.mitpressjournals.org/doi/full/10.1162/rest_a_00778))
- [SVAR_est.m](Estimation_Routines/SVAR_est.m): Least-squares VAR (can also do bias correction)
- [SVAR_IV_est.m](Estimation_Routines/SVAR_IV_est.m): Least-squares SVAR-IV
- [VAR_avg_est.m](Estimation_Routines/VAR_avg_est.m): VAR model averaging ([Hansen, 2016](https://www.ssc.wisc.edu/~bhansen/papers/var.html))

**[DFM](DFM):** Simulation study based on encompassing Dynamic Factor Model (DFM)
- [run_dfm.m](DFM/run_dfm.m): Main file for executing simulations
- [Reporting](DFM/Reporting): Folder with files that produces results figures and tables
  - [run_plot_dgp.m](DFM/Reporting/run_plot_dgp.m): Plots and tables of DGP summary statistics (Table 1 and Figure 3 in the paper)
  - [run_plot_loss.m](DFM/Reporting/run_plot_loss.m): Plots of bias, standard deviation, median bias, and interquartile range (Figures 4-5 and 10-11 in the paper)
  - [run_plot_tradeoff.m](DFM/Reporting/run_plot_tradeoff.m): Plots of head-to-head loss function comparisons and best method (Figures 6-9 in the paper)
  - [settings_shared.m](DFM/Reporting/settings_shared.m): Shared settings for plotting functions
- [Settings](DFM/Settings): Folder with simulation settings
- [Subroutines](DFM/Subroutines): Folder with data for calibrating the DFM and various functions for simulation and computing summary statistics

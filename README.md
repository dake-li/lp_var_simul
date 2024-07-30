# Simulation study of Local Projections, VARs, and related estimators

Matlab code for large-scale simulation studies of impulse response estimators, including Local Projections (LPs), Vector Autoregressions (VARs), and several variants of these

**Reference:**
[Li, Dake](https://github.com/dake-li), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), and [Christian K. Wolf](https://www.christiankwolf.com/) (2024), "Local Projections vs. VARs: Lessons From Thousands of DGPs" (:page_facing_up:[paper](Documents/lp_var_simul.pdf), :bar_chart:[supplement](Documents/lp_var_simul_supplement.pdf))

Tested in: Matlab R2023a on Windows 10 PC (64-bit)

## Contents

**[Documents](Documents):** Paper, supplement, and documentation
- [lp_var_simul.pdf](Documents/lp_var_simul.pdf): Main paper
- [lp_var_simul_supplement.pdf](Documents/lp_var_simul_supplement.pdf): Online supplement
- [lp_var_simul_companion.pdf](Documents/lp_var_simul_companion.pdf): Technical documentation

**[Estimation_Routines](Estimation_Routines):** General-purpose impulse response estimation functions
- [BVAR_est.m](Estimation_Routines/BVAR_est.m): Bayesian VAR ([Giannone, Lenza & Primiceri, 2015](https://doi.org/10.1162/REST_a_00483))
- [LP_est.m](Estimation_Routines/LP_est.m): Least-squares LP (with optional bias correction as in [Herbst & Johanssen, 2023](https://edherbst.net/research/pdfs/hj_lp-revised.pdf))
- [LP_shrink_est.m](Estimation_Routines/LP_shrink_est.m): Penalized LP ([Barnichon & Brownlees, 2019](https://www.mitpressjournals.org/doi/full/10.1162/rest_a_00778))
- [SVAR_est.m](Estimation_Routines/SVAR_est.m): Least-squares VAR (with optional bias correction as in [Pope, 1990](https://doi.org/10.1111/j.1467-9892.1990.tb00056.x))
- [SVAR_IV_est.m](Estimation_Routines/SVAR_IV_est.m): Least-squares SVAR-IV
- [VAR_avg_est.m](Estimation_Routines/VAR_avg_est.m): VAR model averaging ([Hansen, 2016](https://www.ssc.wisc.edu/~bhansen/papers/var.html))

**[DFM](DFM):** Simulation study based on encompassing Dynamic Factor Model (DFM)
- [run_dfm.m](DFM/run_dfm.m): Main file for executing simulations
- [run_combine.m](DFM/run_combine.m): Combines simulation data files from multiple DGPs
- [Reporting](DFM/Reporting): Folder with files that produce results figures and tables
  - [run_plot_dgp.m](DFM/Reporting/run_plot_dgp.m): Plots and tables of DGP summary statistics (Table 1 and Figure 1 in the paper)
  - [run_plot_loss.m](DFM/Reporting/run_plot_loss.m): Plots of bias, standard deviation, median bias, and interquartile range (Figures 2-3 and 10-11 in the paper)
  - [run_plot_tradeoff.m](DFM/Reporting/run_plot_tradeoff.m): Plots of head-to-head loss function comparisons and best method (Figures 4-9 in the paper)
  - [settings_shared.m](DFM/Reporting/settings_shared.m): Shared settings for plotting functions
- [Settings](DFM/Settings): Folder with simulation settings
- [Subroutines](DFM/Subroutines): Folder with data for calibrating the DFM and various functions for simulation and computing summary statistics
  - [SW_DFM_Estimation](DFM/Subroutines/SW_DFM_Estimation): DFM code and data adapted from [Lazarus, Lewis, Stock & Watson (2018)](https://doi.org/10.1080/07350015.2018.1506926)

## Replication

### Main results in this paper

1. **Estimate IRFs from simulated data**: Run the following scripts to pick 6000 DGPs (under observed-shock identification), repeat 5000 Monte Carlo simulations for each DGP, and test multiple estimators for each simulation. This step provides raw IRF estimates of this simulation study.

    - In ``Settings/shared.m``, set ``settings.specifications.random_n_spec = 100``, and ``settings.simul.n_MC = 5000``.
    - In ``run_dfm.m``, first set ``estimand_type = 'ObsShock'``, ``lag_type = 4``, and ``mode_type = 1``.
    - After the setup above, run ``run_dfm.m`` 60 times, by varying the following, to iterate through 6000 DGPs:
      - ``dgp_type`` from ``'G'`` (for fiscal policy type) to ``'MP'`` (for monetary policy type);
      - ``spec_id`` from 1 to 30 (for 30 distinct seeds, where each seed draws 100 random DGPs).
    - Finally, raw IRF estimates will be saved in ``Results/`` directory. (**Warning**: File sizes might be super large.)

2. **Summarize key statistics**: Run the following scripts to obtain summary statistics of raw IRF estimates across 5000 simulations. This step largely reduces dimensionality of results.

    - In ``run_combine.m``, first set ``spec_id_array = [1:30]``, ``dgp_type = 'G'`` (or ``'MP'``). Besides, specify ``estimand_type``, ``lag_type`` and ``mode_type`` to be consistent with Step 1.
    - Finally, run ``run_combine.m`` once to summarize 5000 simulations into 9 summary statistics for each DGP, to simplify results on IRF estimates.
    - These summary statistics will also be saved in ``Results/`` directory.

3. **Examine heterogeneity of selected DGPs and true IRFs**: Run the following scripts to check heterogeneous properties of selected DGPs, and their true IRFs which are derived analytically (See ``Section 3.4``).

    - First in ``Reporting/run_plot_dgp.m``, set ``mode_select = 1``, ``lags_select = 2``, and ``exper_select_group = {[2,5]}``.
    - Then run ``Reporting/run_plot_dgp.m`` to summarize selected DGPs (i.e. ``Table 1``) and plot examples of true IRFs (i.e. ``Figure 1``). (**Warning**: In ``Table 1``, the row of "IV first stage F-statistic" will be computed in Step 5.)
    - Outputs are all saved in ``Reporting/fig/`` directory.

4. **Evaluate loss for observed-shock estimators**: Run the following scripts to show bias & variance profiles of each estimator under observed shock identification, and compare their loss at different bias weights and target horizons (See ``Section 5.1 ~ 5.3``).

    - Always set ``mode_select = 1``, ``lags_select = 2``, and ``exper_select_group = {[2,5]}`` below.
    - Run ``Reporting/run_plot_loss.m`` to get bias and variance profiles separately for each estimator (i.e. ``Figure 2 & 3``).
    - Run ``Reporting/run_plot_tradeoff.m`` to get head-to-head loss comparison btw. two estimators (i.e. ``Figure 4 & 5``, ``Figure 7 ~ 9``).
    - Run python code ``Reporting/plot_best_method.ipynb`` (with ``folder`` set to the output directory in Step 3), to depict the optimal estimator given different bias weights and target horizons (i.e. ``Figure 6``).
    - Outputs are all saved in ``Reporting/fig/`` directory.

5. **Evaluate loss for IV estimators**: Run the following scripts to get bias & variance profiles for estimators under IV identification (i.e. ``Figure 10 & 11``).

    - Redo Step 1 ~ 4, but change ``estimand_type = 'IV'`` in Step 1/2, and ``exper_select_group = {[1,4]}``  in Step 3/4.

### Additional results in the online appendix

6. **Examples of estimated IRFs**: To plot examples of IRF estimates in ``Appendix E``:

    - First finish Step 1.
    - Then in ``Subroutines/run_plot_irf_estimate.m``, set ``spec_id = 1``, ``dgp_type = 'G'``. Besides, specify ``estimand_type``, ``lag_type`` and ``mode_type`` to be consistent with Step 1.
    - Run ``Subroutines/run_plot_irf_estimate.m`` to plot examples of IRF estimates (i.e. ``Figure E.1 ~ E.7``)
    - Outputs are saved in ``Results/`` directory.

7. **Further evaluations for IV estimators**: To further check bias & variance profiles of IV estimators, as in ``Appendix F.1``:

    - ``Figure F.1 & F.2`` should have been generated in the outputs of Step 5.
    - For ``Figure F.3 & F.4``, repeat Step 5, but specify ``DGP_select`` to 2 (low degree of invertibility) or 3 (high degree of invertibility) in ``Reporting/run_plot_loss.m``.

8. **Robustness checks for observed-shock estimators**: Run the following scripts to revisit bias & variance trade-off in some extended cases for observed-shock estimators (See ``Appendix F.2 ~ F.10``).

    - Each of the sub-steps below requires repeating Step 1 ~ 4, but with slight adjustments.
    - For stationary DGPs (``Table F.1``, ``Figure F.5 ~ F.7``), change ``mode_type = 6`` in Step 1/2, and ``mode_select = 6`` in Step 3/4.
    - For recursive identification (``Table F.2``, ``Figure F.8 ~ F.10``), change ``estimand_type = 'Recursive'`` in Step 1/2, and ``exper_select_group = {[3,6]}`` in Step 3/4.
    - For salient observables (``Figure F.11 ~ F.13``), change ``mode_type = 4`` in Step 1/2, and ``mode_select = 4`` in Step 3/4. (**Warning**: In Step 1, code might not run when ``spec_id`` gets large, because the number of available DGPs could have been exhausted.)
    - For 90th percentile loss (``Figure F.14 ~ F.16``), change ``loss_quant = 0.9`` in Step 4.
    - For fiscal and monetary shocks (``Figure F.17 ~ F.20``), change ``exper_select_group = {[2]}`` (to display fiscal shocks only) or ``exper_select_group = {[5]}`` (to display monetary shocks only) in Step 3/4.
    - For longer estimation lag length (``Figure F.21 ~ F.23``), change ``lag_type = 8`` in Step 1/2, and ``lags_select = 3`` in Step 3/4.
    - For smaller sample size (``Figure F.24 ~ F.26``), change ``mode_type = 2`` in Step 1/2, and ``mode_select = 2`` in Step 3/4.
    - For larger sample size and estimation lag length ( ``Figure F.27 ~ F.29``), change ``mode_type = 3`` & ``lag_type = 12`` in Step 1/2, coupled with ``mode_select = 3`` & ``lags_select = 4`` in Step 3/4.
    - For more observables ( ``Figure F.30 ~ F.32``), change ``mode_type = 5`` in Step 1/2, and ``mode_select = 5`` in Step 3/4. Also redo Step 5 for IV estimators.

9. **Splitting by variable categories**: To replicate regression analysis in ``Appendix F.11``:

    - ``Table F.3`` should have been generated in the outputs after finishing Step 1 ~ 4.

## Acknowledgements

We rely on BVAR code by [Domenico Giannone, Michele Lenza & Giorgio Primiceri](http://faculty.wcas.northwestern.edu/gep575/GLPreplicationWeb.zip), penalized LP code by [Regis Barnichon & Christian Brownlees](https://drive.google.com/drive/folders/1Fjzw-U3hjIl467KXywRqeQod2jdHOmDo?usp=sharing), as well as VAR model averaging code by [Bruce Hansen](https://www.ssc.wisc.edu/~bhansen/progs/var.html). We have slightly modified these sets of code to improve their run-time without affecting their numerical output. We also use Dynamic Factor Model code and data by [Eben Lazarus, Daniel Lewis, Jim Stock & Mark Watson](http://www.princeton.edu/~mwatson/ddisk/LLSW_ReplicationFiles_071418.zip).

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049), and Wolf does the same for [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).


%
% OMEN-SED 1.1 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al. (2018) GMD paper + Pika et al. (2020) GMD paper with RCM approximation
% Email: dominik.huelse@ucr.edu or dominik.huelse@gmx.de

% 1.0 updates from 0.9: included pyrite precipitation (implicitly) & gammaH2S dependent on BW oxygenation in Fortran code
% 1.1 updates from 1.0: included the nG approximation of a Reactive Continuum Model (RCM) as described in Pika et al. (2020)
%			- for now we kept the original representation of the 2G-model, but this could be replaced if the nG gives the same results

% 1) Example for executing the model:
% How to calculate a single sediment column with default parameters and plotting results:
% Execute the model from the directory including the MATLAB source code files by typing 
% >> benthic_test.run_OMEN to run the default 2G-model
% >> benthic_test.run_OMEN_RCM to run the RCM approximation
% 1) default sediment-water interface boundary conditions and the parameters for the RCM approximation are specified 
%    in benthic_test.default_swi() and in the properties of benthic_main.m
%    if using the 2G model OM degradation constants need to be changed in benthic_zTOC.m
% 2) the subroutines for the different species are called in benthic_test.test_benthic(1,swi):
% 3) results saved in res and are plotted with benthic_test.plot_column()

% 2) OUTPUT:
% 2.1) Sediment profiles: will be automatically calculated and plotted by
% benthic_test.plot_column(). Two EPS files are saved in the working directory.
% 
% 2.2) Other results like O2, NO3, SO4 penetration depths (i.e. zox, zNO3, zSO4) and simulated
% sediment-water interface fluxes (i.e. flxswiO2, flxswiNO3, flxswiNH4, flxswiSO4,flxswiH2S, flxswiP, flxswiALK, flxswiDIC) 
% are printed in the MATLAB terminal.


% Files of this directory:
%
%   benthic_main    	- Global properties for benthic model (e.g. water-depth, sediment characteristics, stoichiometric factors, functions to calculate sedimentation & bioturbation rate)
%                     	Adjust these parameters here if needed
%   benthic_test    	- Fcts. to execute OMEN-SED, plot output figures and other results and functions for test cases for OMEN-SED & for the RCM approximation
%                     	  benthic_test.default_swi(): define sediment-water interface boundary consitions (e.g. bottom water concentrations)
%                    	  benthic_test.run_OMEN(): execute OMEN-SED
%                     	  benthic_test.plot_column(res, debug, swi, str_date): plot results (res) & save .EPS files ending in 'str_date'  
%			  benthic_test.RCM(): approximation of the RCM with parameters as specified in benthic_test.default_swi()
%   benthic_zTOC    	- Organic matter - 2 Fractions: change the degradation rate constants in this file!
%                     	  also includes the functions dealing with matching across bioturbated boundary (i.e. prepfg_l12 & calcfg_l*)
%   benthic_utils   	- Other utility functions for Generic Boundary Condition Matching 
%   benthic_zALK    	- Solve Alkalinity (including definition of diffusion coefficients)
%   benthic_zDIC    	- Solve DIC (including definition of diffusion coefficients)
%   benthic_zH2S    	- Solve H2S (including definition of diffusion coefficients)
%   benthic_zNH4    	- Solve NH4 (including definition of diffusion coefficients)
%   benthic_zNO3    	- Solve NO3 (including definition of diffusion coefficients)
%   benthic_zO2     	- Solve O2 (including definition of diffusion coefficients)
%   benthic_zPO4_M  	- Solve PO4 and FeP (or M) (including definition of diffusion coefficients)
%   benthic_zSO4    	- Solve SO4 (including definition of diffusion coefficients)
%   fzero_vec       	- Vectorized single-variable nonlinear zero finding based on fzero
%   benthic_zTOC_RCM	- Organic matter - nG fractions as an approximation of an RCM: the type of approximation (i.e. #fractions, values for the parameters a & nu are set in benthic_test.default_swi()) 
%			  and the RCM approximation is computed in benthic_test.RCM()
%                     	  also includes the functions dealing with matching across bioturbated boundary (i.e. prepfg_l12 & calcfg_l*) for the RCM approximation


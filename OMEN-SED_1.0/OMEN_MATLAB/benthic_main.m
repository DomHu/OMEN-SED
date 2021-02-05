% OMEN-SED 1.0 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al (2017) GMD paper

% benthic_main.m
% Global properties for benthic model

classdef benthic_main < handle
    % Global properties for benthic model
    
    properties
        
        ncl;                                    % number of sediment columns
        usescalarcode = true;                   % use scalar code
        
        tol_const = 1e-18;                      % non-zero constant to avoid numerical errors (e.g. division by zero)
        
        %sediment characteristics
        rho_sed=2.5;                            % sediment density (g/cm3)
        wdepth=50.0;                           % water depth (m)
        w;                                      % burial velocity  (cm/yr) - calculated by internal fct. sedrate()
        z0  = 0;                                % surface
        zbio=10.0;                              % bioturbation depth (cm)
        zinf=100;                               %Inifinity (cm)
        loc_BW_O2_anoxia=20E-009;
        Dbio;                                   % bioturbation coefficient (cm2/yr) - calculated by internal fct. biorate()
        por=0.85;                               % porosity (-)
        tort=3.0;                               %tortuosity (-)
        irrigationFactor=1.0;                   %irrigation factor (-)
        dispFactor;                             %dispersion factor (-)
        
        %stoichiometric factors
        X_C;                                    % Carbon Redfield stoichiometry
        Y_N;                                    % Nitrogen Redfield stoichiometry
        Z_P;                                    % Phosphorous Redfield stoichiometry
        SD;                                     % volume factor solid->dissolved phase
        OC;                                     % O2/C (mol/mol)
        NC1;                                    % N/C first TOC fraction (mol/mol)
        NC2;                                    % N/C second TOC fraction (mol/mol)
        PC1;                                    % P/C first TOC fraction (mol/mol)
        PC2;                                    % P/C second TOC fraction (mol/mol)
        FeIIIC;                                 % FeIII/C (mol/mol)
        SO4C;                                   % SO4/C (mol/mol)
        O2Fe2;                                  % O2/Fe2 ratio for oxidation of Fe2 (mol/mol)
        O2H2S;                                  % O2/H2S ratio for oxidation of H2S (mol/mol)
        FeIIIH2S;                            	% Mol of FeIII to oxidize 1 mol H2S
        DICC1;                                  % DIC/C until zSO4 (mol/mol)
        DICC2;                                  % DIC/C below zSO4 (mol/mol)
        MC;                                     % CH4/C (mol/mol)
        gamma=0.05;                           	% fraction of NH4 that is oxidised in oxic layer
        gammaH2S=0.05;                         	% fraction of H2S that is oxidised in oxic layer
        gammaH2SFe=0.0;                       	% fraction of H2S that is oxidised with FeIII  (assume after zFeIII only sulfate reduction, as we don't know zSO4 yet)
        gammaFe2;                               % fraction of Fe2 that is oxidised in oxic layer, moved to Fe2-routine as calculated with S. vd Velde's fit to Cox and BW [O2])
        gammaFe_pp;                         	% fraction of Fe2 that is precipitated as pyrite, moved to Fe2-routine as calculated with S. vd Velde's fit to Cox and BW [O2])
        gammaFeS=0.0;                         	% fraction of H2S that is precipitated as pyrite
        gammaCH4=0.99;                         	% fraction of CH4 that is oxidised at SO4
        satSO4=0.0;                           	% SO4 saturation
        NO3CR;                                  % NO3 consumed by Denitrification
        % Alkalinity production/consumption from:
        ALKROX;                                	% Aerobic degradation
        ALKRNIT;                               	% Nitrification
        ALKRDEN;                             	% Denitrification
        ALKRFeIII;                              % FeIII reduction
        ALKRSUL;                              	% Sulfate reduction
        ALKRH2S;                               	% H2S oxydation
        ALKRFe2;                               	% Fe2 oxydation
        ALKRFeS;                               	% pyrite precipitation
        ALKRMET;                               	% Methanogenesis
        ALKRAOM;                               	% AOM
        
        
        zoxgf = 0.0;                            % cm, rolloff NH4, H2S oxidation for small zox depth (was 0.1)
        
        % Diagnostic output from root finder
        %fzerooptions;
        %fzerooptions = optimset('Display','iter');
        %fzerooptions = optimset('Display','final');
        %fzerooptions = optimset('TolX',0.001);
        fzerooptions = optimset('TolX',100*eps);        % 2.2e-14
        tol_Fe3 = 10^4*eps;
        fzerooptions_Fe3 = optimset('TolX',10^4*eps);   % 2.2e-12
%        fzerooptions_Fe3 = optimset('Display','iter');
    end
    
    methods
        function obj = benthic_main(ncl, wdepth)
            % set default values for the sediment columns
            if nargin > 0
                obj.ncl = ncl;
                obj.wdepth = obj.wdepth*ones(1,obj.ncl);
                obj.zbio = obj.zbio*ones(1,obj.ncl);
                obj.Dbio = obj.Dbio*ones(1,obj.ncl);
                obj.zinf = obj.zinf*ones(1,obj.ncl);
                obj.z0 = obj.z0*ones(1,obj.ncl);
            else
                obj.ncl = 1;
            end
            
            if nargin > 1
                obj.wdepth = wdepth;
            end
            
            obj.usescalarcode = (obj.ncl == 1);
            
            obj.w=benthic_main.sedrate(obj.wdepth);
            obj.Dbio=benthic_main.biorate(obj.wdepth);
            obj.dispFactor=obj.por.^(obj.tort-1.0).*obj.irrigationFactor;	%dispersion factor (-)
            obj.SD=(1-obj.por)./obj.por;
            
            obj.X_C=106;                        % Carbon Redfield stoichiometry
            obj.Y_N=16;                        	% Nitrogen Redfield stoichiometry
            obj.Z_P=1;                        	% Phosphorous Redfield stoichiometry
            obj.OC=(138/106)*obj.SD;          	%O2/C (mol/mol) (was 1.0*obj.SD;)
            obj.NC1= 16.0/106.0*obj.SD;         % N/C first TOC fraction: 16/106 (mol/mol)  % Bohlen 0.067*obj.SD;
            obj.NC2= 16.0/106.0*obj.SD;        	% N/C second TOC fraction (mol/mol)
            obj.PC1=1/106*obj.SD;               % P/C first TOC fraction  1/106 (mol/mol)
            obj.PC2=1/106*obj.SD;               % P/C second TOC fraction 1/106 (mol/mol)
            obj.FeIIIC=4.0;                     % FeIII/C (mol/mol)
            obj.SO4C=(138.0/212.0)*obj.SD;      % SO4/C (mol/mol) (was 0.5*obj.SD;)
            obj.O2Fe2=1.0/4.0;                  % Mol of O2 to oxidize 1 mol Fe2
            obj.O2H2S=2.0;                      % Mol of O2 to oxidize 1 mol H2S
            obj.FeIIIH2S=8.0/obj.SD;            % Mol of FeIII to oxidize 1 mol H2S
            obj.DICC1=1.0*obj.SD;             	% DIC/C until zSO4 (mol/mol)
            obj.DICC2=0.5*obj.SD;             	% DIC/C below zSO4 (mol/mol)
            obj.MC=0.5*obj.SD;                	% CH4/C (mol/mol)
            obj.NO3CR=(94.4/106)*obj.SD;      	% NO3 consumed by Denitrification
            
            %       For explicit Nitrogen (i.e. stand-alone OMEN-SED):
            obj.ALKROX=(obj.Y_N-2*obj.Z_P)/obj.X_C*obj.SD;  	% Aerobic degradation
            obj.ALKRNIT=-2.0;                                   % Nitrification
            obj.ALKRDEN=(4*obj.X_C+3*obj.Y_N-10*obj.Z_P)/(5*obj.X_C)*obj.SD;       	% Denitrification
            obj.ALKRSUL=(obj.X_C+obj.Y_N-2*obj.Z_P)/obj.X_C*obj.SD;                 % Sulfate reduction (was (obj.X_C+obj.Y_N)/obj.X_C*obj.SD;)
            
            obj.ALKRH2S=-2.0;                               % H2S oxydation
            obj.ALKRFeS=-2.0;                               % pyrite precipitation
            obj.ALKRMET= (obj.Y_N-2*obj.Z_P)/obj.X_C;       % Methanogenesis
            obj.ALKRAOM= 2.0;                               % AOM
            
            % % % % 	For implicit Nitrogen (e.g. coupled to GENIE) - change these values:
            % % % %             obj.ALKROX=-(obj.Y_N)/obj.X_C*obj.SD;                        % -(obj.Y_N+2*obj.Z_P)/obj.X_C*obj.SD was +15      % Aerobic degradation
            % % % %             obj.ALKRNIT=0.0;  % no ALK                                             % Nitrification
            % % % %             obj.ALKRDEN=0.0;       %  462/530 was 93.4;    % Denitrification
            % % % %             obj.ALKRSUL=(obj.X_C+obj.Y_N)/obj.X_C*obj.SD;                % 122/106  was 15      	% Sulfato reduction
            
        end
    end
    
    methods(Static)
        
        function w = sedrate(wdepth)
            % sedimentation rate, cm/yr (after Middelburg et al. (1997)) 
            w_Middelburg = 10.0.^(-0.87478367-0.00043512*wdepth)*3.3;
            % sedimentation rate, cm/yr (after Burwicz et al. (2011)) 
            w1 = 0.117;
            w2 = 0.006;
            z1 = 200;
            z2 = 4000;
            c1 = 3;
            c2 = 10;
            w = w1/(1+(wdepth/z1)^c1) + w2/(1+(wdepth/z2)^c2);
%            w = 60/1000;    % To compare with Dale et al. (2015)
        end
        
        function Dbio = biorate(wdepth)
            % bioturbation coeff, cm^2/yr (after Middelburg et al. (1997)) 
            Dbio= 5.2*(10.0^(0.7624-0.0003972*wdepth));	
        end
        
        
        
    end
    
end


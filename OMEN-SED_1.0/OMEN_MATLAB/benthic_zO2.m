classdef benthic_zO2
    % Solve O2
    
    properties
        qdispO2=348.62172;          % O2 diffusion coefficient in water (cm2/yr)
        adispO2=14.08608;         	% O2 linear coefficient for temperature dependence (cm2/yr/oC)
        DO21;                       % O2 diffusion coefficient in bioturbated layer (cm2/yr)
        DO22;                    	% O2 diffusion coefficient in non-bioturbated layer (cm2/yr)
        reac1;                      % reactive terms for O2 consumption - labile frac
        reac2;                      % reactive terms for O2 consumption - refractory frac
    end
    
    methods
        function obj = benthic_zO2(bsd, swi)
            obj.DO21=(obj.qdispO2+obj.adispO2*swi.T)*bsd.dispFactor+bsd.Dbio;
            obj.DO22=(obj.qdispO2+obj.adispO2*swi.T)*bsd.dispFactor;
            
            
            %reactive terms: OM degradation (-) and nitrification (-)
            obj.reac1=-bsd.OC-2*bsd.gamma*bsd.NC1;
            obj.reac2=-bsd.OC-2*bsd.gamma*bsd.NC2;
        end
        
        function r = calc(obj, bsd, swi, r)
            % iteratively solve for zox
            
            %          -ve                        +ve
            fun=@(zox) obj.calcbc(zox,bsd,swi,r,1)  + obj.calcFO2(zox,bsd, swi, r);
            
            % Test for eg zero oxygen at swi
            fun0 = fun(1e-9); % >=0 for eg zero oxygen at swi
            fun1 = fun(1e-12); % >=0 for eg zero oxygen at swi
            funzinf = fun(bsd.zinf); 
            
            % Try zero flux at zinf and see if we have any O2 left
            [flxzox, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r, 2);
            
            if bsd.usescalarcode
                if fun0 >= 0   % i.e. zero oxygen at swi bc O2 flows into sediments (so -)
                    r.zox = 0;
                    bctype = 1;       	% BC: zero concentration
                    conczinf = 0.0;
                elseif conczinf >=0    	% still O2 at zinf -> zox = zinf
                    r.zox = bsd.zinf;
                    bctype = 2;      	% BC: zero flux
                else                    % search zox in the interval
                    bctype = 1;         % BC: zero concentration
                    r.zox=fzero(fun,[1e-9 bsd.zinf],bsd.fzerooptions);
                    conczinf = 0.0;
                end
            else % same logic, in vector form
                lzinf = concinf >=0;  % true for each x(i) that still has O2 at zinf
                
                zox=fzero_vec(fun,1e-10, bsd.zinf,bsd.fzerooptions);
                lz0 = fun0 >= 0;      % true for each x(i) with ~zero oxygen at swi
                
                bctype = 1*~lzinf + 2*lzinf;
                r.zox  = 0*lz0 + bsd.zinf.*lzinf + zox.*(~lz0 & ~lzinf); % project out appropriate solution
            end
            
            [flxzox, conczox, flxswiO2, r] = obj.calcbc(r.zox, bsd, swi, r, bctype);
            %format longEng
            r.flxswiO2_DOU = flxswiO2;
            %advflux0 = bsd.por.*bsd.w.*(swi.O20)
            %advfluxinf = bsd.por.*bsd.w.*(conczinf)
            flxswiO2 = flxswiO2 - bsd.por.*bsd.w.*(swi.O20-conczinf);
            r.flxzox = flxzox;
            r.conczox = conczox;
            r.flxswiO2 = flxswiO2;
            
            % OUTPUT FOR FLUX of reduxed substances at z_ox
            Fredzox = obj.calcFO2(r.zox,bsd, swi, r);
%            fprintf('zox = %g Approx F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  r.zox, Fredzox);
            
        end
        
        function [flxzox, conczox, flxswi,r] = calcbc(obj, zox, bsd, swi, r, bctype)
            % calculate solution for given zox
            if(swi.TwoG_OM_model)
                % Preparation: sort out solution-matching across bioturbation boundary (if necessary)
                rO2.ls = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, obj.reac2, 0, bsd.z0, zox, obj.DO21, obj.DO22);
                
                % basis functions at upper boundary
                [ e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0] ...
                    = r.zTOC.calcfg_l12(bsd.z0, bsd, swi, r, obj.reac1, obj.reac2, 0, rO2.ls);
                % ... and lower boundary
                [ e_zox, dedz_zox, f_zox, dfdz_zox, g_zox, dgdz_zox] ...
                    = r.zTOC.calcfg_l12(zox, bsd, swi, r, obj.reac1, obj.reac2, 0, rO2.ls);
            else
                % Preparation: sort out solution-matching across bioturbation boundary (if necessary)
                rO2.ls = r.zTOC_RCM.prepfg_l12(bsd, swi, r, obj.reac1, 0, bsd.z0, zox, obj.DO21, obj.DO22);
                
                % basis functions at upper boundary
                [ e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0] ...
                    = r.zTOC_RCM.calcfg_l12(bsd.z0, bsd, swi, r, obj.reac1, 0, rO2.ls);
                % ... and lower boundary
                [ e_zox, dedz_zox, f_zox, dfdz_zox, g_zox, dgdz_zox] ...
                    = r.zTOC_RCM.calcfg_l12(zox, bsd, swi, r, obj.reac1, 0, rO2.ls);
                
                
            end
            
            
            % Solve for AO2, BO2 given boundary conditions (expressed in terms of transformed soln)
            
            % Case 1 zero concentration at zox
            % AO2*e_zox   +   BO2*f_zox  + g_zox = 0;
            % AO2*e_0     +   BO2*f_0     + g_0  = swi.O20;
            
            % | e_zox f_zox |  |AO2|   = | -g_zox         |
            % | e_0   f_0   |  |BO2|     | swi.O20 - gt_0 |
            
            [ bctype1_AO2, bctype1_BO2]      = benthic_utils.solve2eqn(e_zox, f_zox, e_0, f_0, -g_zox, swi.O20 - g_0);
            
            % Case  2 zero flux at zox
            % AO2*dedz_zox +  BO2*dfz_zox + dgz_zox = 0;
            % AO2*e_0     +   BO2*f_0     + g_0     = swi.O20;
            [ bctype2_AO2, bctype2_BO2]      = benthic_utils.solve2eqn(dedz_zox, dfdz_zox, e_0, f_0, -dgdz_zox, swi.O20 - g_0);
            
            % Choose type of solution requested (vectorized form)
            rO2.AO2 = (bctype==1).*bctype1_AO2 + (bctype==2).*bctype2_AO2;
            rO2.BO2 = (bctype==1).*bctype1_BO2 + (bctype==2).*bctype2_BO2;
            
            Dzox = (zox < bsd.zbio).*obj.DO21 + (zox >= bsd.zbio).*obj.DO22;
            flxzox =  Dzox.*(rO2.AO2.*dedz_zox+rO2.BO2.*dfdz_zox + dgdz_zox);                    % no por factor as this is per cm^2 pore area
            conczox = rO2.AO2.*e_zox+rO2.BO2.*f_zox + g_zox;
            
            Dswi = (0 < bsd.zbio).*obj.DO21 + (0 >= bsd.zbio).*obj.DO22;
            flxswi = bsd.por.*(Dswi.*(rO2.AO2.*dedz_0+rO2.BO2.*dfdz_0 + dgdz_0)); % - bsd.w.*swi.O20);   % por fac so this is per cm^2 water column
            
            r.zxf = zox./(bsd.zoxgf + zox); % roll off oxidation at low zox
            
            r.rO2 = rO2;
            
        end
        
        
        
        function FO2 = calcFO2(obj, zox, bsd, swi, r)
            % Oxydation of reduced species at zox (NEED A RATIO for ODU! and add NH4 adsporption!
            % O2H2S = 2.0 mole of O2 to oxidize 1 mol H2S
            tmpreac1=bsd.gammaH2S*bsd.O2H2S*bsd.SO4C+2*bsd.gamma*bsd.NC1;
            tmpreac2=bsd.gammaH2S*bsd.O2H2S*bsd.SO4C+2*bsd.gamma*bsd.NC2;
            % tmpreac1=0.2;
            % tmpreac2=0.2;
            %tmpreac1=bsd.OC+2*bsd.gamma*bsd.NC1;
            %tmpreac2=bsd.OC+2*bsd.gamma*bsd.NC2;
            %FLUX of NH4 and Reduced species from ZOX to ZINF
            
            %            FO2 = 0.0; % no secondary redox!
            if(swi.TwoG_OM_model)
                FO2 = zox./(bsd.zoxgf + zox).*r.zTOC.calcReac(zox, bsd.zinf, tmpreac1, tmpreac2, bsd, swi, r) + (1-bsd.gammaFe_pp)*bsd.gammaFe2*bsd.O2Fe2.*swi.Flux_FeIII0;  % taking into account approximated Fe2-flux to zox, just as 'gammaFe2 * O2Fe2* Influx_Fe3+'
                % FO2 = zox./(bsd.zoxgf + zox).*r.zTOC.calcReac(zox, bsd.zinf, tmpreac1, tmpreac2, bsd, swi, r);                                  % not taking into account iron re-oxidation
            else
                if(swi.Iron)
                FO2 = zox./(bsd.zoxgf + zox).*r.zTOC_RCM.calcReac(zox, bsd.zinf, tmpreac1, bsd, swi, r) + (1-bsd.gammaFe_pp)*bsd.gammaFe2*bsd.O2Fe2.*swi.Flux_FeIII0;        
                % taking into account approximated Fe2-flux to zox, just as 'gammaFe2 * O2Fe2 * Influx_Fe3+'
                else
                    FO2 = zox./(bsd.zoxgf + zox).*r.zTOC_RCM.calcReac(zox, bsd.zinf, tmpreac1, bsd, swi, r);        % not taking into account iron re-oxidation
                end
            end
            % NB (1-bsd.por)/bsd.por  has been included in OC etc stoich factors, so this is flux / cm^2 pore area
            
        end
    
    function FO2 = calcFO2_exact(obj, zox, bsd, swi, r)
    % Oxydation of reduced species at zox (CALCULATED WITH ACTUAL PENETRATION DEPTHS)
    
    tmpreac1_N=2*bsd.gamma*bsd.NC1;
    tmpreac2_N=2*bsd.gamma*bsd.NC2;
    tmpreac1_S=bsd.gammaH2S*(1-bsd.gammaFeS)*bsd.O2H2S*bsd.SO4C;
    tmpreac2_S=bsd.gammaH2S*(1-bsd.gammaFeS)*bsd.O2H2S*bsd.SO4C;
    
    %tmpreac1=bsd.OC+2*bsd.gamma*bsd.NC1;
    %tmpreac2=bsd.OC+2*bsd.gamma*bsd.NC2;
    %FLUX of NH4 and Reduced species from ZOX to ZINF
    
    if(swi.TwoG_OM_model)
        FO2 = zox./(bsd.zoxgf + zox).*(r.zTOC.calcReac(r.zno3, bsd.zinf, tmpreac1_N, tmpreac2_N, bsd, swi, r) ...
            + r.zTOC.calcReac(r.zno3, r.zso4, tmpreac1_S, tmpreac2_S, bsd, swi, r));
        % NB (1-bsd.por)/bsd.por  has been included in OC etc stoich factors, so this is flux / cm^2 pore area
    else
        FO2 = zox./(bsd.zoxgf + zox).*(r.zTOC_RCM.calcReac(r.zno3, bsd.zinf, tmpreac1_N, bsd, swi, r) ...
            + r.zTOC_RCM.calcReac(r.zno3, r.zso4, tmpreac1_S, bsd, swi, r));
        % NB (1-bsd.por)/bsd.por  has been included in OC etc stoich factors, so this is flux / cm^2 pore area
        
    end
    end
    
    function [O2, flxO2, flxO2D, flxO2adv] = calcO2(obj, z, bsd, swi, r)
    % Calculate O2 conc and flux at depth z from solution
    if z <= r.zox    % <= so handle case of fully oxic with zox = zinf
        % basis functions at z
        if(swi.TwoG_OM_model)
            [ e, dedz, f, dfdz, g, dgdz] ...
                = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1, obj.reac2, 0, r.rO2.ls);
        else
            [ e, dedz, f, dfdz, g, dgdz] ...
                = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, obj.reac1, 0, r.rO2.ls);
        end
        if z < bsd.zbio  % < so handle zbio = 0
            D = obj.DO21;
        else
            D = obj.DO22;
        end
        
        O2 = r.rO2.AO2.*e + r.rO2.BO2.*f + g;
        flxO2D = bsd.por.*D.*(r.rO2.AO2.*dedz+r.rO2.BO2.*dfdz + dgdz); % diffusive component
        flxO2adv = - bsd.por.*bsd.w*O2;                                % advective component
        flxO2 = flxO2D + flxO2adv;                                     % total
        
    else
        O2 = 0;
        flxO2 = 0;
        flxO2D = 0;
        flxO2adv = 0;
    end
    end
    
    
end

end


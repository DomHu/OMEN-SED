classdef benthic_zSO4
    % Solve SO4
    
    properties
        qdispSO4=157.68;            % SO4 diffusion coefficient in water (cm2/yr)
        adispSO4=7.884;             % SO4 linear coefficient for temperature dependence (cm2/yr/oC)
        DSO41;                      % SO4 diffusion coefficient in bioturbated layer (cm2/yr)
        DSO42;                      % SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        reac1;
% %         reac2;
    end
    
    methods
        function obj = benthic_zSO4(bsd, swi)
            obj.DSO41=(obj.qdispSO4+obj.adispSO4*swi.T).*bsd.dispFactor+bsd.Dbio;  	% SO4 diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DSO42=(obj.qdispSO4+obj.adispSO4*swi.T).*bsd.dispFactor;          	% SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
            
            %reactive terms: OM degradation
            obj.reac1=-bsd.SO4C;
%             obj.reac2=-bsd.SO4C;
            
        end
        
        function r = calc(obj, bsd, swi, r)
            % Iteratively solve for zso4
            
            % try zero flux at zinf and see if we have any SO4 left, also
            % calculate [SO4] at zinf for advective loss
            [flxzso4, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r, 2);
            
            
            if r.zno3 == bsd.zinf
                r.zso4 = bsd.zinf;
                bctype = 2;
            else
                
                fun=@(zso4)-obj.calcbc(zso4,bsd,swi,r,1) - obj.calcFSO4(zso4,bsd, swi, r);
                
                %             % try zero flux at zinf and see if we have any SO4 left
                %             [flxzso4, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r, 2);
                if bsd.usescalarcode
                    if conczinf >=0
                        r.zso4 = bsd.zinf;
                        bctype = 2;
                    else
                        bctype = 1;
                        conczinf = 0.0;
                        funzno3=fun(r.zno3);
                        funzinf=fun(bsd.zinf);
                        r.zso4=fzero(fun,[max(r.zno3, 1e-10), bsd.zinf],bsd.fzerooptions);
                    end
                else  % vectorized version
                    bctype = (conczinf < 0)*1 + (conczso4>=0)*2;
                    zso4=fzero_vec(fun,max(r.zno3, 1e-10), bsd.zinf,bsd.fzerooptions);
                    r.zso4 = (bctype==1).*zso4 + (bctype==2).*bsd.zinf;
                end
                
            end
            [flxzso4, conczso4, flxswiSO4, r] = obj.calcbc(r.zso4, bsd, swi, r, bctype);    % Dom18.05.2016: not necessary for bctype 2 (done in line 32 already)
            
            flxswiSO4 = flxswiSO4 - bsd.por.*bsd.w.*(swi.SO40-conczinf);
            if(abs(flxswiSO4) <= bsd.tol_const)
                flxswiSO4 = 0.0;
            end
            
            r.flxzso4 = flxzso4;
            r.conczso4 = conczso4;
            r.flxswiSO4 = flxswiSO4;
        end
        
        function [flxzso4, conczso4, flxswi,r] = calcbc(obj, zso4, bsd, swi, r, bctype)
            % Calculate trial solution for given zso4, matching boundary conditions from layer-by-layer solutions
            
            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zox, passive diffn
            %      ls =      prepfg_l12( bsd, swi, r, reac1,   ktemp, zU, zL, D1,        D2)
            rSO4.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, 0,       0,     0, r.zox, obj.DSO41, obj.DSO42);
            % layer 2: zox < z < zno3, passive diffn
            rSO4.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, 0,       0,  r.zox, r.zno3, obj.DSO41, obj.DSO42);
            % layer 3: zno3 < z < zso4, SO4 consumption by OM oxidation
            rSO4.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, 0, r.zno3, zso4, obj.DSO41, obj.DSO42);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 3 zso4
            [ e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4] ...
                = r.zTOC.calcfg_l12(zso4, bsd, swi, r, obj.reac1,  0, rSO4.ls3);
            
            % Match at zno3, layer 2 - layer 3 (continuity and flux)
            % basis functions at bottom of layer 2
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r,     0,             0, rSO4.ls2);
            % ... and top of layer 3
            [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r, obj.reac1,  0, rSO4.ls3);
            % match solutions at zno3 - continuous concentration and flux
            [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                0, 0);
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
            % flux of H2S to oxic interface (Source of SO4)
            % NB: include methane region as AOM will produce sulphide as well..
            %            FH2S = 0.0; %r.zTOC.calcReac(r.zno3, zso4, bsd.SO4C, bsd.SO4C, bsd, swi, r) + 0.0; % no secondary redox!
            FH2S = r.zTOC.calcReac(r.zno3, zso4, bsd.SO4C, bsd, swi, r) ... % MULTIPLY BY 1/POR ????
                + bsd.gammaCH4.*r.zTOC.calcReac(zso4, bsd.zinf, bsd.MC, bsd, swi, r); % Dominik 25.02.2016
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0, 0, rSO4.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0, 0, rSO4.ls2);
            % transform to use coeffs from l3
            [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
            
            % match solutions at zox - continuous concentration, flux discontinuity from H2S ox
            D = (r.zox <= bsd.zbio).*obj.DSO41 + (r.zox > bsd.zbio).*obj.DSO42;
            
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, -r.zxf.*bsd.gammaH2S*FH2S./D);
            %Dom 09.02.2016: is there a ...*gammaH2S*FH2S... missing?
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, 0 , 0, rSO4.ls1);
            % transform to use coeffs from l3
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            
            
            % Find solutions for two possible types of lower bc
            %  case 1  zero concentration at zso4
            % Solve for ASO4, BSO4 given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
            % ASO4*e3_zso4   +  BSO4*f3_zso4  + g3_zso4 = 0;
            % ASO4*e1_0     +   BSO4*f1_0     + g1_0  = swi.SO40;
            
            % | e3_zso4 f3_zso4 |  |ASO4|   = | -g3_zso4       |
            % | e1_0     f1_0   |  |BSO4|     | swi.SO40 - g1_0 |
            [ bctype1_A3, bctype1_B3]      = benthic_utils.solve2eqn(e3_zso4, f3_zso4, e1_0, f1_0, -g3_zso4, swi.SO40 - g1_0);
            
            %  case  2 zero flux at zso4
            % ASO4*de3dz_zso4   +  BSO4*dfdz3_zso4  + dgdz3_zso4 = 0;
            % ASO4*e1_0         +   BSO4*f1_0       + g1_0       = swi.SO40;
            [ bctype2_A3, bctype2_B3]      = benthic_utils.solve2eqn(dedz3_zso4, dfdz3_zso4, e1_0, f1_0, -dgdz3_zso4, swi.SO40 - g1_0);
            
            % Choose type of solution requested (vectorized form)
            rSO4.A3 = (bctype==1).*bctype1_A3 + (bctype==2).*bctype2_A3;
            rSO4.B3 = (bctype==1).*bctype1_B3 + (bctype==2).*bctype2_B3;
            
            % calculate conc and flux at zso4
            conczso4 = rSO4.A3.*e3_zso4+rSO4.B3.*f3_zso4 + g3_zso4;
            D = (zso4 <= bsd.zbio).*obj.DSO41 + (zso4 > bsd.zbio).*obj.DSO42;
            flxzso4 = D.*(rSO4.A3.*dedz3_zso4+rSO4.B3.*dfdz3_zso4 + dgdz3_zso4);        % includes 1/por ie flux per (cm^2 pore area)
            
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            flxswi = bsd.por.*(obj.DSO41.*(rSO4.A3.*dedz1_0+rSO4.B3.*dfdz1_0 + dgdz1_0)); % - bsd.w.*swi.SO40);   % NB: use A3, B3 as these are _xformed_ layer 1 basis functions
            
            % save coeffs for layers 2 and 1
            rSO4.A2 = zno3.a.*rSO4.A3 + zno3.b.*rSO4.B3 + zno3.e;
            rSO4.B2 = zno3.c.*rSO4.A3 + zno3.d.*rSO4.B3 + zno3.f;
            
            rSO4.A1 = zox.a.*rSO4.A3 + zox.b.*rSO4.B3 + zox.e;
            rSO4.B1 = zox.c.*rSO4.A3 + zox.d.*rSO4.B3 + zox.f;
            
            
            r.rSO4 = rSO4;
            
        end
        
        
        function FSO4 = calcFSO4(obj, zso4, bsd, swi, r)
            % Calculate SO4 consumption below zso4, by organic matter and indirectly via methane oxidation
            
            tmpreac1    = bsd.MC.*bsd.gammaCH4;
%            tmpreac2    = bsd.MC.*bsd.gammaCH4;
            %            FSO4 = 0.0;    % no secondary redox!
            FSO4 = r.zTOC.calcReac(zso4, bsd.zinf, tmpreac1, bsd, swi, r);
            % TODO confirm (1-bsd.por)*  has been added (to k1 & k2 ?)
        end
        
        
        function [SO4, flxSO4] = calcSO4(obj, z, bsd, swi, r)
            % Calculate SO4 concentration and flux at depth z from solution
            
            rSO4 = r.rSO4;
            if z <= r.zso4
                if z <= bsd.zbio
                    D = obj.DSO41;
                else
                    D = obj.DSO42;
                end
                
                if z <= r.zox   % layer 1
                    [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 ,  0, rSO4.ls1);
                    SO4     = r.rSO4.A1.*e + r.rSO4.B1.*f + g;
                    flxSO4  = D.*(r.rSO4.A1.*dedz+r.rSO4.B1.*dfdz + dgdz);
                elseif z <= r.zno3 % layer 2
                    [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 ,  0, rSO4.ls2);
                    SO4     = r.rSO4.A2.*e + r.rSO4.B2.*f + g;
                    flxSO4  = D.*(r.rSO4.A2.*dedz+r.rSO4.B2.*dfdz + dgdz);
                else
                    [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1, 0, rSO4.ls3);
                    SO4     = r.rSO4.A3.*e + r.rSO4.B3.*f + g;
                    flxSO4  = D.*(r.rSO4.A3.*dedz+r.rSO4.B3.*dfdz + dgdz);
                end
                
            else
                
                SO4 = 0;
                flxSO4 = 0;
            end
        end
        
    end
    
end


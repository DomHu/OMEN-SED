classdef benthic_zFe2
    % Solve Fe2
    
    properties
        qdispFe2=107.538;  % 307.476; %        % Fe2 diffusion coefficient in water (cm2/yr)
        adispFe2=4.768; % 9.636;  %           % Fe2 linear coefficient for temperature dependence (cm2/yr/oC)
        DFe21;                      % Fe2 diffusion coefficient in bioturbated layer (cm2/yr)
        DFe22;                      % Fe2 diffusion coefficient in non-bioturbated layer (cm2/yr)
        FeIrr_scale=0.2;            % Irrigation of Fe2+ is scaled to 20% compared to other solutes due to its high affinity for oxidation on burrow walls

        
        reac1;
        reac2;
        
    end
    
    methods
        function obj = benthic_zFe2(bsd, swi, Cox)
            obj.DFe21=(obj.qdispFe2+obj.adispFe2*swi.T).*bsd.dispFactor.*obj.FeIrr_scale+bsd.Dbio;  	% Fe2 diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DFe22=(obj.qdispFe2+obj.adispFe2*swi.T).*bsd.dispFactor.*obj.FeIrr_scale;          	% Fe2 diffusion coefficient in non-bioturbated layer (cm2/yr)
            
            %reactive terms: OM degradation
            obj.reac1=bsd.FeIIIC*bsd.SD;    % the *SD is needed here as Fe2 is dissolved
            obj.reac2=bsd.FeIIIC*bsd.SD;

            
        end
        
        function r = calc(obj, bsd, swi, r)
            % Calculate Fe2
            
            if(swi.TwoG_OM_model)
                % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zox, passive diffn
            %      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)
            rFe2.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,         0,     0, r.zox, obj.DFe21, obj.DFe22);
            % layer 2: zox < z < zno3, passive diffn
            rFe2.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,         0,  r.zox, r.zno3, obj.DFe21, obj.DFe22);
            % layer 3: zno3 < z < zfeIII, Fe2 prod. by OM oxidation
            rFe2.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, obj.reac2, 0, r.zno3, r.zfeIII, obj.DFe21, obj.DFe22);
            % layer 4: zfeIII < z < zinf, passive diffn
            rFe2.ls4 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,         0, r.zfeIII, bsd.zinf, obj.DFe21, obj.DFe22);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 4 zinf
            [ e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf] ...
                = r.zTOC.calcfg_l12(bsd.zinf, bsd, swi, r, 0, 0, 0, rFe2.ls4);
            
            % Match at zfeIII, layer 3 - layer 4 (continuity and flux with AOM production)
            % basis functions at bottom of layer 3
            [ e3_zfeIII, dedz3_zfeIII, f3_zfeIII, dfdz3_zfeIII, g3_zfeIII, dgdz3_zfeIII] ...
                = r.zTOC.calcfg_l12(r.zfeIII, bsd, swi, r,  obj.reac1, obj.reac2, 0, rFe2.ls3);
            % ... and top of layer 4
            [ e4_zfeIII, dedz4_zfeIII, f4_zfeIII, dfdz4_zfeIII, g4_zfeIII, dgdz4_zfeIII] ...
                = r.zTOC.calcfg_l12(r.zfeIII, bsd, swi, r,  0,  0, 0, rFe2.ls4);
            %flux of Fe2 produced below zFeIII - oxidation of H2S by FeIII (Source of Fe2)
            zfeIIIFFe2 = 0.0;  % no Fe2 production below zFeIII -- DH - TODO:  Check 24.07.20
%            zfeIIIFFe2 = r.zTOC.calcReac(r.zfeIII, bsd.zinf, bsd.MC, bsd.MC, bsd, swi, r); % MULTIPLY BY 1/POR ????
            % match solutions at zfeIII - continuous concentration and flux
           	D  = (r.zfeIII <= bsd.zbio).*obj.DFe21 + (r.zfeIII > bsd.zbio).*obj.DFe22;
            [zfeIII.a, zfeIII.b, zfeIII.c, zfeIII.d, zfeIII.e, zfeIII.f] = benthic_utils.matchsoln(e3_zfeIII, f3_zfeIII, g3_zfeIII, dedz3_zfeIII, dfdz3_zfeIII, dgdz3_zfeIII, ...
                e4_zfeIII, f4_zfeIII, g4_zfeIII, dedz4_zfeIII, dfdz4_zfeIII, dgdz4_zfeIII, ...
                0, zfeIIIFFe2./D);
%                0, -bsd.gammaCH4.*zfeIIIFFe2./D);  % use - for source of Fe2

            % Match at zno3, layer 2 - layer 3 (continuity and flux)
            % basis functions at bottom of layer 2
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r,     0,            0, 0, rFe2.ls2);
            % ... and top of layer 3
            [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r, obj.reac1, obj.reac2, 0, rFe2.ls3);
            % ... transformed to use coeffs from l4
            [e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3] = benthic_utils.xformsoln(e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                zfeIII.a , zfeIII.b , zfeIII.c , zfeIII.d , zfeIII.e ,zfeIII.f);
            % match solutions at zno3 - continuous concentration and flux
            [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                0, 0);
            
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from Fe2 source)
            %flux of Fe2 to oxic interface (from all sources of Fe2 below)
%            zoxFFe2 = 0.0;  % no secondary redox! -- DH - TODO: Check 24.07.20
            zoxFFe2 = r.zTOC.calcReac(r.zno3, r.zfeIII, obj.reac1, obj.reac2, bsd, swi, r);   
%                + r.zTOC.calcReac(r.zfeIII, bsd.zinf, bsd.MC, bsd.MC, bsd, swi, r); % Dominik 25.02.2016
            
            % Dom 24.02.2016: actually should be 2 integrals for Fe2 produced: SO4-reduction + AOM (see documentation, but has the same reac const = 0.5) :
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0 , 0 , 0, rFe2.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0, 0, 0, rFe2.ls2);
            % transform to use coeffs from l4
            [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
            
            % match solutions at zox - continuous concentration, flux discontinuity from Fe2 ox
            
            D = (r.zox <= bsd.zbio).*obj.DFe21 + (r.zox > bsd.zbio).*obj.DFe22;
            
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, r.zxf.*bsd.gammaFe2.*zoxFFe2./D);
                % with pyrite:      0, r.zxf.*bsd.gammaFe2.*zoxFFe2./D);
               	% without pyrite:   0, r.zxf.*bsd.gammaFe2.*zoxFFe2./D);
                
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, 0 , 0 , 0, rFe2.ls1);
            % transform to use coeffs from l4
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            else
                % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zox, passive diffn
            rFe2.ls1 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, 0,         0,         0, r.zox, obj.DFe21, obj.DFe22);
            % layer 2: zox < z < zno3, passive diffn
            rFe2.ls2 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, 0,         0,  r.zox, r.zno3, obj.DFe21, obj.DFe22);
            % layer 3: zno3 < z < zfeIII, Fe2 prod. by OM oxidation
            rFe2.ls3 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, obj.reac1, 0, r.zno3, r.zfeIII, obj.DFe21, obj.DFe22);
            % layer 4: zfeIII < z < zinf, passive diffn
            rFe2.ls4 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, 0,         0, r.zfeIII, bsd.zinf, obj.DFe21, obj.DFe22);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 4 zinf
            [ e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf] ...
                = r.zTOC_RCM.calcfg_l12(bsd.zinf, bsd, swi, r, 0, 0, rFe2.ls4);
            
            % Match at zfeIII, layer 3 - layer 4 (continuity and flux with AOM production)
            % basis functions at bottom of layer 3
            [ e3_zfeIII, dedz3_zfeIII, f3_zfeIII, dfdz3_zfeIII, g3_zfeIII, dgdz3_zfeIII] ...
                = r.zTOC_RCM.calcfg_l12(r.zfeIII, bsd, swi, r,  obj.reac1, 0, rFe2.ls3);
            % ... and top of layer 4
            [ e4_zfeIII, dedz4_zfeIII, f4_zfeIII, dfdz4_zfeIII, g4_zfeIII, dgdz4_zfeIII] ...
                = r.zTOC_RCM.calcfg_l12(r.zfeIII, bsd, swi, r,  0,  0, rFe2.ls4);
             %flux of Fe2 consumed at zFeIII, reacts with flux of H2S from below is precipitated as pyrite (Sink of Fe2)
            zfeIIIFFe2 = r.zTOC_RCM.calcReac(r.zno3, r.zfeIII, obj.reac1, bsd, swi, r)*bsd.gammaFe_pp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NEW IDEA 06.05.2022
            % now also add Fe2+ from SMI!   
            % But this would need to be accounted for in H2S and SO4 as
            % well... get's too messy, just use lower Fe3+ input which
            % results in lower Fe2+ return and higher H2S return as it is
            % not oxidized
%             zfeIIIFFe2_reac = r.zTOC_RCM.calcReac(r.zno3, r.zfeIII, obj.reac1, bsd, swi, r);
%             zfeIIIFFe2_SMI = swi.Flux_Fe2_from_SMI;
%             zfeIIIFFe2 = (zfeIIIFFe2_reac + zfeIIIFFe2_SMI)*bsd.gammaFe_pp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%            zfeIIIFFe2 = r.swi.Flux_FeIII0*1.0;  % 
%            zfeIIIFH2S = r.zTOC_RCM.calcReac(r.zfeIII, bsd.zinf, bsd.SO4C, bsd, swi, r); % assume zso4 = zinf
            %flux of Fe2 produced below zFeIII - oxidation of H2S by FeIII (Source of Fe2)
%            zfeIIIFFe2 = r.zTOC_RCM.calcReac(r.zfeIII, bsd.zinf, bsd.MC, bsd.MC, bsd, swi, r); % MULTIPLY BY 1/POR ????
            % match solutions at zfeIII - continuous concentration and flux
            D  = (r.zfeIII <= bsd.zbio).*obj.DFe21 + (r.zfeIII > bsd.zbio).*obj.DFe22;
            [zfeIII.a, zfeIII.b, zfeIII.c, zfeIII.d, zfeIII.e, zfeIII.f] = benthic_utils.matchsoln(e3_zfeIII, f3_zfeIII, g3_zfeIII, dedz3_zfeIII, dfdz3_zfeIII, dgdz3_zfeIII, ...
                e4_zfeIII, f4_zfeIII, g4_zfeIII, dedz4_zfeIII, dfdz4_zfeIII, dgdz4_zfeIII, ...
                0, zfeIIIFFe2./D);
%                0, -bsd.gammaCH4.*zfeIIIFFe2./D);  % use - for source of Fe2

            % Match at zno3, layer 2 - layer 3 (continuity and flux)
            % basis functions at bottom of layer 2
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC_RCM.calcfg_l12(r.zno3, bsd, swi, r,     0,            0, rFe2.ls2);
            % ... and top of layer 3
            [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                = r.zTOC_RCM.calcfg_l12(r.zno3, bsd, swi, r, obj.reac1, 0, rFe2.ls3);
            % ... transformed to use coeffs from l4
            [e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3] = benthic_utils.xformsoln(e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                zfeIII.a , zfeIII.b , zfeIII.c , zfeIII.d , zfeIII.e ,zfeIII.f);
            % match solutions at zno3 - continuous concentration and flux
            [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                0, 0);
            
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from Fe2 source)
            % flux of Fe2 to oxic interface (from all sources of Fe2 below)
            % NB: include methane region as AOM will produce sulphide as well..
%            zoxFFe2 = 0.0;  % no secondary redox! -- DH - TODO: Check 24.07.20
            zoxFFe2 = r.zTOC_RCM.calcReac(r.zno3, r.zfeIII, obj.reac1, bsd, swi, r)*(1-bsd.gammaFe_pp);   
%                + r.zTOC_RCM.calcReac(r.zfeIII, bsd.zinf, bsd.MC, bsd.MC, bsd, swi, r); % Dominik 25.02.2016
            
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC_RCM.calcfg_l12(r.zox, bsd, swi, r, 0 , 0, rFe2.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC_RCM.calcfg_l12(r.zox, bsd, swi, r, 0, 0, rFe2.ls2);
            % transform to use coeffs from l4
            [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
            
            % match solutions at zox - continuous concentration, flux discontinuity from Fe2 ox
            
            D = (r.zox <= bsd.zbio).*obj.DFe21 + (r.zox > bsd.zbio).*obj.DFe22;
            
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, r.zxf.*bsd.gammaFe2*zoxFFe2./D);
                % with pyrite:      0, r.zxf.*bsd.gammaFe2*zoxFFe2./D);
               	% without pyrite:   0, r.zxf.*bsd.gammaFe2.*zoxFFe2./D);
                
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC_RCM.calcfg_l12(0, bsd, swi, r, 0 , 0, rFe2.ls1);
            % transform to use coeffs from l4
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);                
            end
            
            % Solve for AFe2, BFe2 given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
            % AFe2*dedz4_zinf   +  BFe2*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
            % AFe2*e1_0     +   BFe2*f1_0     + g1_0  = swi.Fe20;
            
            % | dedz4_zinf dfdz4_zinf |  |AFe2|   = | -dgz4_zinf       |
            % | e1_0     f1_0         |  |BFe2|     | swi.Fe20 - g1_0 |
            
            [ rFe2.A4, rFe2.B4]      = benthic_utils.solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, swi.Fe20 - g1_0);
            
            %             % calculate conc and flux at zfeIII
            %             r.conczfeIIIh2s = rFe2.A4.*e4_zfeIII+rFe2.B4.*f4_zfeIII + g4_zfeIII;
            % calculate concentration at zinf
            r.conczinfFe2 = rFe2.A4.*e4_zinf+rFe2.B4.*f4_zinf + g4_zinf;
         	
%             if(r.conczinfFe2 < 0)
%                 warning('NEGATIVE conczinfFe2 -- set to zero!');
%                 r.conczinfFe2 = 0.0;
%            end
            
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            r.flxswiFe2 = bsd.por.*(obj.DFe21.*(rFe2.A4.*dedz1_0+rFe2.B4.*dfdz1_0 + dgdz1_0) - bsd.w.*(swi.Fe20 - r.conczinfFe2));   % NB: use A4, B4 as these are _xformed_ layer 1 basis functions

            if(r.flxswiFe2 < 0)
                warning('NEGATIVE flxswiFe2 -- set to zero!');
                r.flxswiFe2 = 0.0;
           end                
            % save coeffs for layers 3, 2 and 1
            rFe2.A3 = zfeIII.a.*rFe2.A4 + zfeIII.b.*rFe2.B4 + zfeIII.e;
            rFe2.B3 = zfeIII.c.*rFe2.A4 + zfeIII.d.*rFe2.B4 + zfeIII.f;
            
            rFe2.A2 = zno3.a.*rFe2.A4 + zno3.b.*rFe2.B4 + zno3.e;
            rFe2.B2 = zno3.c.*rFe2.A4 + zno3.d.*rFe2.B4 + zno3.f;
            
            rFe2.A1 = zox.a.*rFe2.A4 + zox.b.*rFe2.B4 + zox.e;
            rFe2.B1 = zox.c.*rFe2.A4 + zox.d.*rFe2.B4 + zox.f;
            
            
            r.rFe2 = rFe2;
            
            
        end
        
        
        
        function [Fe2, flxFe2, e, dedz, f, dfdz, g, dgdz] = calcFe2(obj, z, bsd, swi, r)
            % Calculate Fe2 concentration and flux at depth z from solution
            
            rFe2 = r.rFe2;
            
            if z <= bsd.zbio
                D = obj.DFe21;
            else
                D = obj.DFe22;
            end
            
            if(swi.TwoG_OM_model)
            if z <= r.zox   % layer 1
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0 , 0, rFe2.ls1);
                Fe2     = r.rFe2.A1.*e + r.rFe2.B1.*f + g;
                flxFe2  = D.*(r.rFe2.A1.*dedz+r.rFe2.B1.*dfdz + dgdz);
            elseif z <= r.zno3 % layer 2
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0 , 0, rFe2.ls2);
                Fe2     = r.rFe2.A2.*e + r.rFe2.B2.*f + g;
                flxFe2  = D.*(r.rFe2.A2.*dedz+r.rFe2.B2.*dfdz + dgdz);
            elseif z <= r.zfeIII
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1, obj.reac2 , 0, rFe2.ls3);
                Fe2     = r.rFe2.A3.*e + r.rFe2.B3.*f + g;
                flxFe2  = D.*(r.rFe2.A3.*dedz+r.rFe2.B3.*dfdz + dgdz);
            else
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0,       0 , 0, rFe2.ls4);
                Fe2     = r.rFe2.A4.*e + r.rFe2.B4.*f + g;
                flxFe2  = D.*(r.rFe2.A4.*dedz+r.rFe2.B4.*dfdz + dgdz);
            end
            else
              	if z <= r.zox   % layer 1
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0 , 0, rFe2.ls1);
                Fe2     = r.rFe2.A1.*e + r.rFe2.B1.*f + g;
                flxFe2  = D.*(r.rFe2.A1.*dedz+r.rFe2.B1.*dfdz + dgdz);
            elseif z <= r.zno3 % layer 2
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0 , 0, rFe2.ls2);
                Fe2     = r.rFe2.A2.*e + r.rFe2.B2.*f + g;
                flxFe2  = D.*(r.rFe2.A2.*dedz+r.rFe2.B2.*dfdz + dgdz);
            elseif z <= r.zfeIII
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, obj.reac1, 0, rFe2.ls3);
                Fe2     = r.rFe2.A3.*e + r.rFe2.B3.*f + g;
                flxFe2  = D.*(r.rFe2.A3.*dedz+r.rFe2.B3.*dfdz + dgdz);
            else
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0, 0, rFe2.ls4);
                Fe2     = r.rFe2.A4.*e + r.rFe2.B4.*f + g;
                flxFe2  = D.*(r.rFe2.A4.*dedz+r.rFe2.B4.*dfdz + dgdz);
            end
        end
            
        end
        
    end
    
end


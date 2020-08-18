classdef benthic_zALK
    % Solve ALK
    
    properties
        qdispALK=151.69;        	% ALK diffusion coefficient in water (cm2/yr)
        adispALK=7.93;          	% ALK linear coefficient for temperature dependence (cm2/yr/oC)
        DALK1;                      % ALK diffusion coefficient in bioturbated layer (cm2/yr)
        DALK2;                      % ALK diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        KNH4=1.3;                   % Adsorption coefficient (same in ocix and anoxic layer) (-)
        
        reac11;                     % z < zox:  Nitrification (-2) Aerobic degradation (+15/106)
        reac12;                     % z < zox:  Nitrification (-2) Aerobic degradation (+15/106)
        reac21;                     % zox < z < zno3: Denitrification (+93.4/106)
        reac22;                     % zox < z < zno3: Denitrification (+93.4/106)
        reac3;                      % zno3 < z < zso4: Sulfate reduction (+15/106)
        reac4;                      % zso4 < z < zinf: Methanogenesis (+14/106)
    end
    
    methods
        function obj = benthic_zALK(bsd, swi)
            obj.DALK1=(obj.qdispALK+obj.adispALK*swi.T).*bsd.dispFactor+bsd.Dbio;   	% ALK diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DALK2=(obj.qdispALK+obj.adispALK*swi.T).*bsd.dispFactor;              	% ALK diffusion coefficient in non-bioturbated layer (cm2/yr)
            
            
            %reactive terms: OM degradation
            % Following is too high!
            obj.reac11=bsd.gamma*bsd.NC1/(1+obj.KNH4)*bsd.ALKRNIT+bsd.ALKROX;       % z < zox:  Nitrification (-2) Aerobic degradation (+15/106)
            obj.reac12=bsd.gamma*bsd.NC2/(1+obj.KNH4)*bsd.ALKRNIT+bsd.ALKROX;       % z < zox:  Nitrification (-2) Aerobic degradation (+15/106)
            obj.reac21=bsd.ALKRDEN;                                               	% zox < z < zno3: Denitrification (+93.4/106)
            obj.reac22=bsd.ALKRDEN;                                                	% zox < z < zno3: Denitrification (+93.4/106)
            obj.reac3=bsd.ALKRSUL;                                                  % zno3 < z < zso4: Sulfate reduction (+15/106)
            obj.reac4=bsd.ALKRMET;                                                  % zso4 < z < zinf: Methanogenesis (+14/106)
            
        end
        
        function r = calc(obj, bsd, swi, r)
            % Calculate ALK
            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zox, Nitrification (-) Aerobic degradation (+)
            %      ls =      prepfg_l12( bsd, swi, r, reac1,           ktemp,     zU,   zL,       D1,      D2)
            rALK.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac11,     0,     0, r.zox, obj.DALK1, obj.DALK2);
            % layer 2: zox < z < zno3, Denitrification (+)
            rALK.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac21,     0,  r.zox, r.zno3, obj.DALK1, obj.DALK2);
            % layer 3: zno3 < z < zso4, Sulfate reduction (+)
            rALK.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac3,      0, r.zno3, r.zso4, obj.DALK1, obj.DALK2);
            % layer 4: zso4 < z < zinf, Methanogenesis (+)
            rALK.ls4 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac4,      0, r.zso4, bsd.zinf, obj.DALK1, obj.DALK2);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 4 zinf
            [ e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf] ...
                = r.zTOC.calcfg_l12(bsd.zinf, bsd, swi, r, obj.reac4, 0, rALK.ls4);
            
            % Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
            % basis functions at bottom of layer 3
            [ e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4] ...
                = r.zTOC.calcfg_l12(r.zso4,  bsd, swi, r,  obj.reac3, 0, rALK.ls3);
            % ... and top of layer 4
            [ e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4] ...1
                = r.zTOC.calcfg_l12(r.zso4, bsd, swi, r,   obj.reac4, 0, rALK.ls4);
            %flux of ALK produced by AOM interface (Source of ALK)
            %            zso4FALK = 0.0;     % no secondary redox!
            zso4FALK = bsd.ALKRAOM.*bsd.gammaCH4.*r.zTOC.calcReac(r.zso4, bsd.zinf, bsd.SD,  bsd, swi, r); % MULTIPLY BY 1/POR ????
            % match solutions at zso4 - continuous concentration and flux
            [zso4.a, zso4.b, zso4.c, zso4.d, zso4.e, zso4.f] = benthic_utils.matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, ...
                e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, ...
                0, -zso4FALK./obj.DALK2);   % DOM TODO
            %Dom 24.02.2016: No *(1-gammaCH4)*zso4FALK... missing here!
            % Match at zno3, layer 2 - layer 3 (continuity and flux)
            % basis functions at bottom of layer 2
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r,     obj.reac21, 0, rALK.ls2);
            % ... and top of layer 3
            [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r, obj.reac3, 0, rALK.ls3);
            % ... transformed to use coeffs from l4
            [e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3] = benthic_utils.xformsoln(e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                zso4.a , zso4.b , zso4.c , zso4.d , zso4.e ,zso4.f);
            % match solutions at zno3 - continuous concentration and flux
            [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                0, 0);
            
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from ALK source)
            %flux of ALK to oxic interface (from all sources of ALK below)
            % from NH4 and H2S
            %            zoxFALK = 0.0;     % no secondary redox!
            zoxFALK = bsd.ALKRNIT.*bsd.gamma.*r.zTOC.calcReac(r.zno3, bsd.zinf, bsd.NC1/(1+obj.KNH4), bsd, swi, r) ... % MULTIPLY BY 1/POR ????
                + bsd.ALKRH2S.*bsd.gammaH2S.*r.zTOC.calcReac(r.zno3, r.zso4, bsd.SO4C, bsd, swi, r); % Dominik 25.02.2016
            %Dom 24.02.2016: actually should be 2 integrals for ALK produced: ALK-reduction + AOM (see documentation, but has the same reac const = 0.5) :
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, obj.reac11,  0, rALK.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, obj.reac21,  0, rALK.ls2);
            % transform to use coeffs from l4
            [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
            
            % match solutions at zox - continuous concentration, flux discontinuity from ALK ox
            
            D = (r.zox <= bsd.zbio).*obj.DALK1 + (r.zox > bsd.zbio).*obj.DALK2;
            
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, -r.zxf.*zoxFALK./D); % DOM TODO gammaH2S???
            % Dominik 24.02.2016 think it should be -r.zxf.*(1-bsd.gammaALK).*zoxFALK./D -> but changes profile significantly!
            % Dominik 24.02.2016 was r.zxf.*zoxFALK./D    need gammaALK here!
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, obj.reac11,      0, rALK.ls1);
            % transform to use coeffs from l4
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            
            
            % Solve for AALK, BALK given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
            % AALK*dedz4_zinf   +  BALK*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
            % AALK*e1_0     +   BALK*f1_0     + g1_0  = swi.ALK0;
            
            % | dedz4_zinf dfdz4_zinf |  |AALK|   = | -dgz4_zinf       |
            % | e1_0     f1_0         |  |BALK|     | swi.ALK0 - g1_0 |
            
            [ rALK.A4, rALK.B4]      = benthic_utils.solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, swi.ALK0 - g1_0);
            
            % calculate concentration at zinf
            r.conczinfALK = rALK.A4.*e4_zinf+rALK.B4.*f4_zinf + g4_zinf;
            
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            r.flxswiALK = bsd.por.*(obj.DALK1.*(rALK.A4.*dedz1_0+rALK.B4.*dfdz1_0 + dgdz1_0) - bsd.w.*(swi.ALK0 - r.conczinfALK));   % NB: use A4, B4 as these are _xformed_ layer 1 basis functions
            
            % save coeffs for layers 3, 2 and 1
            rALK.A3 = zso4.a.*rALK.A4 + zso4.b.*rALK.B4 + zso4.e;
            rALK.B3 = zso4.c.*rALK.A4 + zso4.d.*rALK.B4 + zso4.f;
            
            rALK.A2 = zno3.a.*rALK.A4 + zno3.b.*rALK.B4 + zno3.e;
            rALK.B2 = zno3.c.*rALK.A4 + zno3.d.*rALK.B4 + zno3.f;
            
            rALK.A1 = zox.a.*rALK.A4 + zox.b.*rALK.B4 + zox.e;
            rALK.B1 = zox.c.*rALK.A4 + zox.d.*rALK.B4 + zox.f;
            
            
            r.rALK = rALK;
            
            
        end
        
        
        
        function [ALK, flxALK, e, dedz, f, dfdz, g, dgdz] = calcALK(obj, z, bsd, swi, r)
            % Calculate ALK concentration and flux at depth z from solution
            
            rALK = r.rALK;
            
            if z <= bsd.zbio
                D = obj.DALK1;
            else
                D = obj.DALK2;
            end
            
            if z <= r.zox   % layer 1
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac11, 0, rALK.ls1);
                ALK     = r.rALK.A1.*e + r.rALK.B1.*f + g;
                flxALK  = D.*(r.rALK.A1.*dedz+r.rALK.B1.*dfdz + dgdz);
            elseif z <= r.zno3 % layer 2
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac21, 0, rALK.ls2);
                ALK     = r.rALK.A2.*e + r.rALK.B2.*f + g;
                flxALK  = D.*(r.rALK.A2.*dedz+r.rALK.B2.*dfdz + dgdz);
            elseif z <= r.zso4
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac3,  0, rALK.ls3);
                ALK     = r.rALK.A3.*e + r.rALK.B3.*f + g;
                flxALK  = D.*(r.rALK.A3.*dedz+r.rALK.B3.*dfdz + dgdz);
            else
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac4,  0, rALK.ls4);
                ALK     = r.rALK.A4.*e + r.rALK.B4.*f + g;
                flxALK  = D.*(r.rALK.A4.*dedz+r.rALK.B4.*dfdz + dgdz);
            end
        end
        
    end
    
end


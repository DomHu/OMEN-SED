classdef benthic_zH2S
    % Solve H2S
    
    properties
        qdispH2S=307.476;           % H2S diffusion coefficient in water (cm2/yr)
        adispH2S=9.636;             % H2S linear coefficient for temperature dependence (cm2/yr/oC)
        DH2S1;                      % H2S diffusion coefficient in bioturbated layer (cm2/yr)
        DH2S2;                      % H2S diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        reac1;
%        reac2;
    end
    
    methods
        function obj = benthic_zH2S(bsd, swi)
            obj.DH2S1=(obj.qdispH2S+obj.adispH2S*swi.T).*bsd.dispFactor+bsd.Dbio;  	% H2S diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DH2S2=(obj.qdispH2S+obj.adispH2S*swi.T).*bsd.dispFactor;          	% H2S diffusion coefficient in non-bioturbated layer (cm2/yr)
            
            %reactive terms: OM degradation
            obj.reac1=bsd.SO4C;
 %           obj.reac2=bsd.SO4C;
            
        end
        
        function r = calc(obj, bsd, swi, r)
            % Calculate H2S
            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zox, passive diffn
            %      ls =      prepfg_l12( bsd, swi, r, reac1,     ktemp, zU, zL, D1,        D2)
            rH2S.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,     0, r.zox, obj.DH2S1, obj.DH2S2);
            % layer 2: zox < z < zno3, passive diffn
            rH2S.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,  r.zox, r.zno3, obj.DH2S1, obj.DH2S2);
            % layer 3: zno3 < z < zso4, H2S prod. by OM oxidation
            rH2S.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, 0, r.zno3, r.zso4, obj.DH2S1, obj.DH2S2);
            % layer 4: zso4 < z < zinf, passive diffn
            rH2S.ls4 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0, r.zso4, bsd.zinf, obj.DH2S1, obj.DH2S2);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 4 zinf
            [ e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf] ...
                = r.zTOC.calcfg_l12(bsd.zinf, bsd, swi, r, 0, 0, rH2S.ls4);
            
            % Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
            % basis functions at bottom of layer 3
            [ e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4] ...
                = r.zTOC.calcfg_l12(r.zso4, bsd, swi, r,  obj.reac1, 0, rH2S.ls3);
            % ... and top of layer 4
            [ e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4] ...
                = r.zTOC.calcfg_l12(r.zso4, bsd, swi, r,   0, 0, rH2S.ls4);
            %flux of H2S produced by AOM interface (Source of H2S)
            %            zso4FH2S = 0.0;  % no secondary redox!
            zso4FH2S = r.zTOC.calcReac(r.zso4, bsd.zinf, bsd.MC, bsd, swi, r); % MULTIPLY BY 1/POR ????
            % match solutions at zso4 - continuous concentration and flux
            [zso4.a, zso4.b, zso4.c, zso4.d, zso4.e, zso4.f] = benthic_utils.matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, ...
                e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, ...
                0, -bsd.gammaCH4.*zso4FH2S./obj.DH2S2);
            %Dom 24.02.2016: No *(1-gammaCH4)*zso4FH2S... missing here!
            % Match at zno3, layer 2 - layer 3 (continuity and flux)
            % basis functions at bottom of layer 2
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r,     0,     0, rH2S.ls2);
            % ... and top of layer 3
            [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r, obj.reac1, 0, rH2S.ls3);
            % ... transformed to use coeffs from l4
            [e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3] = benthic_utils.xformsoln(e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                zso4.a , zso4.b , zso4.c , zso4.d , zso4.e ,zso4.f);
            % match solutions at zno3 - continuous concentration and flux
            [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                0, 0);
            
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
            %flux of H2S to oxic interface (from all sources of H2S below)
            % NB: include methane region as AOM will produce sulphide as well..
            % zoxFH2S = 0.0; %r.zTOC.calcReac(r.zno3, r.zso4, bsd.SO4C, bsd.SO4C, bsd, swi, r) + 0.0; % no secondary redox!
            zoxFH2S = r.zTOC.calcReac(r.zno3, r.zso4, bsd.SO4C, bsd, swi, r) ... % MULTIPLY BY 1/POR ????
                + r.zTOC.calcReac(r.zso4, bsd.zinf, bsd.MC, bsd, swi, r); % Dominik 25.02.2016
            
            % Dom 24.02.2016: actually should be 2 integrals for H2S produced: SO4-reduction + AOM (see documentation, but has the same reac const = 0.5) :
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0 , 0, rH2S.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0, 0, rH2S.ls2);
            % transform to use coeffs from l4
            [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
            
            % match solutions at zox - continuous concentration, flux discontinuity from H2S ox
            
            D = (r.zox <= bsd.zbio).*obj.DH2S1 + (r.zox > bsd.zbio).*obj.DH2S2;
            
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, r.zxf.*bsd.gammaH2S.*zoxFH2S./D);
            % Dominik 24.02.2016 think it should be -r.zxf.*(1-bsd.gammaH2S).*zoxFH2S./D -> but changes profile significantly!
            % Dominik 24.02.2016 was r.zxf.*zoxFH2S./D    need gammaH2S here!
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, 0 , 0, rH2S.ls1);
            % transform to use coeffs from l4
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            
            
            % Solve for AH2S, BH2S given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
            % AH2S*dedz4_zinf   +  BH2S*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
            % AH2S*e1_0     +   BH2S*f1_0     + g1_0  = swi.H2S0;
            
            % | dedz4_zinf dfdz4_zinf |  |AH2S|   = | -dgz4_zinf       |
            % | e1_0     f1_0         |  |BH2S|     | swi.H2S0 - g1_0 |
            
            [ rH2S.A4, rH2S.B4]      = benthic_utils.solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, swi.H2S0 - g1_0);
            
            %             % calculate conc and flux at zso4
            %             r.conczso4h2s = rH2S.A4.*e4_zso4+rH2S.B4.*f4_zso4 + g4_zso4;
            % calculate concentration at zinf
            r.conczinfH2S = rH2S.A4.*e4_zinf+rH2S.B4.*f4_zinf + g4_zinf;
            
            
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            r.flxswiH2S = bsd.por.*(obj.DH2S1.*(rH2S.A4.*dedz1_0+rH2S.B4.*dfdz1_0 + dgdz1_0) - bsd.w.*(swi.H2S0 - r.conczinfH2S));   % NB: use A4, B4 as these are _xformed_ layer 1 basis functions
            
            % save coeffs for layers 3, 2 and 1
            rH2S.A3 = zso4.a.*rH2S.A4 + zso4.b.*rH2S.B4 + zso4.e;
            rH2S.B3 = zso4.c.*rH2S.A4 + zso4.d.*rH2S.B4 + zso4.f;
            
            rH2S.A2 = zno3.a.*rH2S.A4 + zno3.b.*rH2S.B4 + zno3.e;
            rH2S.B2 = zno3.c.*rH2S.A4 + zno3.d.*rH2S.B4 + zno3.f;
            
            rH2S.A1 = zox.a.*rH2S.A4 + zox.b.*rH2S.B4 + zox.e;
            rH2S.B1 = zox.c.*rH2S.A4 + zox.d.*rH2S.B4 + zox.f;
            
            
            r.rH2S = rH2S;
            
            
        end
        
        
        
        function [H2S, flxH2S, e, dedz, f, dfdz, g, dgdz] = calcH2S(obj, z, bsd, swi, r)
            % Calculate H2S concentration and flux at depth z from solution
            
            rH2S = r.rH2S;
            
            if z <= bsd.zbio
                D = obj.DH2S1;
            else
                D = obj.DH2S2;
            end
            
            if z <= r.zox   % layer 1
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0, rH2S.ls1);
                H2S     = r.rH2S.A1.*e + r.rH2S.B1.*f + g;
                flxH2S  = D.*(r.rH2S.A1.*dedz+r.rH2S.B1.*dfdz + dgdz);
            elseif z <= r.zno3 % layer 2
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0, rH2S.ls2);
                H2S     = r.rH2S.A2.*e + r.rH2S.B2.*f + g;
                flxH2S  = D.*(r.rH2S.A2.*dedz+r.rH2S.B2.*dfdz + dgdz);
            elseif z <= r.zso4
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1,  0, rH2S.ls3);
                H2S     = r.rH2S.A3.*e + r.rH2S.B3.*f + g;
                flxH2S  = D.*(r.rH2S.A3.*dedz+r.rH2S.B3.*dfdz + dgdz);
            else
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0, rH2S.ls4);
                H2S     = r.rH2S.A4.*e + r.rH2S.B4.*f + g;
                flxH2S  = D.*(r.rH2S.A4.*dedz+r.rH2S.B4.*dfdz + dgdz);
            end
        end
        
    end
    
end


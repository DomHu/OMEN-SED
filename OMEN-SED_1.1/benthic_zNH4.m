classdef benthic_zNH4
    % Solve NH4
    
    properties
        qdispNH4=309.0528;       	% NH4 diffusion coefficient in water (cm2/yr)
        adispNH4=12.2640;         	% NH4 linear coefficient for temperature dependence (cm2/yr/oC)
        DNH41;                      % NH4 diffusion coefficient in bioturbated layer (cm2/yr)
        DNH42;                      % NH4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        KNH4=1.3;                   % Adsorption coefficient (same in ocix and anoxic layer) (-)
        
    end
    
    methods
        function obj = benthic_zNH4(bsd, swi)
            obj.DNH41=((obj.qdispNH4+obj.adispNH4*swi.T).*bsd.dispFactor+bsd.Dbio)/(1+obj.KNH4); 	% NH4 diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DNH42=((obj.qdispNH4+obj.adispNH4*swi.T).*bsd.dispFactor)/(1+obj.KNH4);           	% NH4 diffusion coefficient in non-bioturbated layer (cm2/yr)
            
        end
        
        function r = calc(obj, bsd, swi, r)
            
            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zox, NH4 prod (remaining after oxidation)
            %      ls =      prepfg_l12( bsd, swi, r,        reac1,                     ,ktemp, zU, zL,       D1,        D2)
            rNH4.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, (1-bsd.gamma)*bsd.NC1/(1+obj.KNH4),0, 0, r.zox, obj.DNH41, obj.DNH42);
            % Dominik 23.02.2016: NC1 and NC2 missing in react1/2 ???
            % layer 2: zox < z < zno3, passive diffn TODO NH4 from denitrification?
            rNH4.ls2 = r.zTOC.prepfg_l12(bsd, swi, r,          0                        ,0, r.zox, r.zno3, obj.DNH41, obj.DNH42);
            % layer 3: zno3 < z < zinf, NH4 production
            rNH4.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, bsd.NC1/(1+obj.KNH4)              , 0, r.zno3, bsd.zinf, obj.DNH41, obj.DNH42);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 3 zinf
            [ e3_zinf, dedz3_zinf, f3_zinf, dfdz3_zinf, g3_zinf, dgdz3_zinf] ...
                = r.zTOC.calcfg_l12(bsd.zinf, bsd, swi, r, bsd.NC1/(1+obj.KNH4)         , 0, rNH4.ls3);
            
            % Match at zno3, layer 2 - layer 3 (continuity and flux)
            % basis functions at bottom of layer 2
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r,     0                          ,0, rNH4.ls2);
            % ... and top of layer 3
            [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r, bsd.NC1/(1+obj.KNH4)           , 0, rNH4.ls3);
            % match solutions at zno3 - continuous concentration and flux
            [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                0, 0);
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 sink)
            %flux of NH4 to oxic interface  TODO NH4 prod by denitrification?
            FNH4 = r.zTOC.calcReac(r.zno3, bsd.zinf, bsd.NC1/(1+obj.KNH4), bsd, swi, r); % MULTIPLY BY 1/POR ????
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r,(1-bsd.gamma)*bsd.NC1/(1+obj.KNH4), 0, rNH4.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0                                , 0, rNH4.ls2);
            % transform to use coeffs from l3
            [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
            % match solutions at zox - continuous concentration, flux discontinuity from NH4 ox
            
            D  = (r.zox <= bsd.zbio).*obj.DNH41 + (r.zox > bsd.zbio).*obj.DNH42;
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, r.zxf.*bsd.gamma.*FNH4./D);
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, (1-bsd.gamma)*bsd.NC1/(1+obj.KNH4)   , 0, rNH4.ls1);
            % transform to use coeffs from l3
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            
            
            % Solve for ANH4, BNH4 given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
            % ANH4*dedz3_zinf   +  BNH4*dfdz3_zinf  + dgdz3_zinf = 0;
            % ANH4*e1_0     +   BNH4*f1_0     + g1_0  = swi.NH40;
            
            % | dedz3_zinf dfdz3_zinf |  |ANH4|   = | -dgdz3_zinf       |
            % | e1_0     f1_0         |  |BNH4|     | swi.NH40 - g1_0 |
            
            [ rNH4.A3, rNH4.B3]      = benthic_utils.solve2eqn(dedz3_zinf, dfdz3_zinf, e1_0, f1_0, -dgdz3_zinf, swi.NH40 - g1_0);
            
            
            % calculate concentration at zinf
            r.conczinfNH4 = rNH4.A3.*e3_zinf+rNH4.B3.*f3_zinf + g3_zinf;
            
            % flux at swi - DO include por so this is per cm^2 water column area
            r.flxswiNH4 = bsd.por.*(obj.DNH41.*(rNH4.A3.*dedz1_0+rNH4.B3.*dfdz1_0 + dgdz1_0) - bsd.w.*(swi.NH40 - r.conczinfNH4));
            
            % save coeffs for layers 2 and 1
            rNH4.A2 = zno3.a.*rNH4.A3 + zno3.b.*rNH4.B3 + zno3.e;
            rNH4.B2 = zno3.c.*rNH4.A3 + zno3.d.*rNH4.B3 + zno3.f;
            
            rNH4.A1 = zox.a.*rNH4.A3 + zox.b.*rNH4.B3 + zox.e;
            rNH4.B1 = zox.c.*rNH4.A3 + zox.d.*rNH4.B3 + zox.f;
            
            
            r.rNH4 = rNH4;
            
        end
        
        
        function [NH4, flxNH4] = calcNH4(obj, z, bsd, swi, r)
            % Calculate NH4 concentration and flux at depth z from solution
            
            rNH4 = r.rNH4;
            
            if z <= bsd.zbio
                D = obj.DNH41;
            else
                D = obj.DNH42;
            end
            
            if z <= r.zox   % layer 1
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, (1-bsd.gamma)*bsd.NC1/(1+obj.KNH4) , 0, rNH4.ls1);
                NH4     = r.rNH4.A1.*e + r.rNH4.B1.*f + g;
                flxNH4  = D.*(r.rNH4.A1.*dedz+r.rNH4.B1.*dfdz + dgdz);
            elseif z <= r.zno3 % layer 2
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0                                  , 0, rNH4.ls2);
                NH4     = r.rNH4.A2.*e + r.rNH4.B2.*f + g;
                flxNH4  = D.*(r.rNH4.A2.*dedz+r.rNH4.B2.*dfdz + dgdz);
            else
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, bsd.NC1/(1+obj.KNH4)               , 0, rNH4.ls3);
                NH4     = r.rNH4.A3.*e + r.rNH4.B3.*f + g;
                flxNH4  = D.*(r.rNH4.A3.*dedz+r.rNH4.B3.*dfdz + dgdz);
            end
            
            
            
        end
        
    end
    
end


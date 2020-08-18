classdef benthic_zDIC
    % Solve DIC
    
    properties
        qdispDIC=151.69;        % DIC diffusion coefficient in water (cm2/yr)
        adispDIC=7.93;      	% DIC linear coefficient for temperature dependence (cm2/yr/oC)
        DDIC1;              	% DIC diffusion coefficient in bioturbated layer (cm2/yr)
        DDIC2;              	% DIC diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        reac1;
        reac2;
    end
    
    methods
        function obj = benthic_zDIC(bsd, swi)
            obj.DDIC1=(obj.qdispDIC+obj.adispDIC*swi.T).*bsd.dispFactor+bsd.Dbio;   	% DIC diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DDIC2=(obj.qdispDIC+obj.adispDIC*swi.T).*bsd.dispFactor;            	% DIC diffusion coefficient in non-bioturbated layer (cm2/yr)
            
            
            %reactive terms: OM degradation
            obj.reac1=bsd.DICC1;                       %DIC/C until zSO4 (mol/mol)
            obj.reac2=bsd.DICC2;                       %DIC/C below zSO4 (mol/mol)
            
        end
        
        function r = calc(obj, bsd, swi, r)
            % Calculate DIC
            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            % layer 1: 0 < z < zso4, DIC production by OM degradation
            %      ls =      prepfg_l12( bsd, swi, r, reac1    , ktemp, zU, zL,    D1,        D2)
            rDIC.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, 0,     0, r.zso4, obj.DDIC1, obj.DDIC2);
            % layer 2: zso4 < z < zinf, DIC production by OM degradation (Methanogenesis) -> different production rate
            rDIC.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac2, 0,   r.zso4, bsd.zinf, obj.DDIC1, obj.DDIC2);
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 2 zinf
            [ e2_zinf, dedz2_zinf, f2_zinf, dfdz2_zinf, g2_zinf, dgdz2_zinf] ...
                = r.zTOC.calcfg_l12(bsd.zinf, bsd, swi, r, obj.reac2, 0, rDIC.ls2);
            
            % Match at zso4, layer 1 - layer 2 (continuity and flux with AOM production)
            % basis functions at bottom of layer 1
            [ e1_zso4, dedz1_zso4, f1_zso4, dfdz1_zso4, g1_zso4, dgdz1_zso4] ...
                = r.zTOC.calcfg_l12(r.zso4, bsd, swi, r,  obj.reac1, 0, rDIC.ls1);
            % ... and top of layer 2
            [ e2_zso4, dedz2_zso4, f2_zso4, dfdz2_zso4, g2_zso4, dgdz2_zso4] ...
                = r.zTOC.calcfg_l12(r.zso4, bsd, swi, r,  obj.reac2, 0, rDIC.ls2);
            %flux of DIC produced by AOM interface (Source of DIC)
            zso4FDIC = r.zTOC.calcReac(r.zso4, bsd.zinf, bsd.MC, bsd, swi, r); % MULTIPLY BY 1/POR ????
            % match solutions at zso4 - continuous concentration and flux
            [zso4.a, zso4.b, zso4.c, zso4.d, zso4.e, zso4.f] = benthic_utils.matchsoln(e1_zso4, f1_zso4, g1_zso4, dedz1_zso4, dfdz1_zso4, dgdz1_zso4, ...
                e2_zso4, f2_zso4, g2_zso4, dedz2_zso4, dfdz2_zso4, dgdz2_zso4, ...
                0, -bsd.gammaCH4.*zso4FDIC./obj.DDIC2);
            %Dom 24.02.2016: No *(1-gammaCH4)*zso4FDIC... missing here!
            
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, obj.reac1 , 0, rDIC.ls1);
            % transform to use coeffs from l2
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zso4.a , zso4.b , zso4.c , zso4.d , zso4.e ,zso4.f);
            
            
            % Solve for ADIC, BDIC given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
            % ADIC*dedz2_zinf   +  BDIC*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
            % ADIC*e1_0     +   BDIC*f1_0     + g1_0  = swi.DIC0;
            
            % | dedz2_zinf dfdz2_zinf |  |ADIC|   = | -dgz4_zinf       |
            % | e1_0     f1_0         |  |BDIC|     | swi.DIC0 - g1_0 |
            
            [ rDIC.A2, rDIC.B2]      = benthic_utils.solve2eqn(dedz2_zinf, dfdz2_zinf, e1_0, f1_0, -dgdz2_zinf, swi.DIC0 - g1_0);
            
            
            % calculate concentration at zinf
            r.conczinfDIC = rDIC.A2.*e2_zinf+rDIC.B2.*f2_zinf + g2_zinf;
            
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            r.flxswiDIC = bsd.por.*(obj.DDIC1.*(rDIC.A2.*dedz1_0+rDIC.B2.*dfdz1_0 + dgdz1_0) - bsd.w.*(swi.DIC0 - r.conczinfDIC));   % NB: use A2, B2 as these are _xformed_ layer 1 basis functions
            
            % save coeffs for layers 1
            rDIC.A1 = zso4.a.*rDIC.A2 + zso4.b.*rDIC.B2 + zso4.e;
            rDIC.B1 = zso4.c.*rDIC.A2 + zso4.d.*rDIC.B2 + zso4.f;
            
            %             rDIC.A2 = zno3.a.*rDIC.A4 + zno3.b.*rDIC.B4 + zno3.e;
            %             rDIC.B2 = zno3.c.*rDIC.A4 + zno3.d.*rDIC.B4 + zno3.f;
            %
            %             rDIC.A1 = zox.a.*rDIC.A4 + zox.b.*rDIC.B4 + zox.e;
            %             rDIC.B1 = zox.c.*rDIC.A4 + zox.d.*rDIC.B4 + zox.f;
            
            
            r.rDIC = rDIC;
            
            
        end
        
        
        
        function [DIC, flxDIC, e, dedz, f, dfdz, g, dgdz] = calcDIC(obj, z, bsd, swi, r)
            % Calculate DIC concentration and flux at depth z from solution
            
            rDIC = r.rDIC;
            
            if z <= bsd.zbio
                D = obj.DDIC1;
            else
                D = obj.DDIC2;
            end
            
            if z <= r.zso4   % layer 1
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1 , 0, rDIC.ls1);
                DIC     = r.rDIC.A1.*e + r.rDIC.B1.*f + g;
                flxDIC  = D.*(r.rDIC.A1.*dedz+r.rDIC.B1.*dfdz + dgdz);
            else
                [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac2 ,  0, rDIC.ls2);
                DIC     = r.rDIC.A2.*e + r.rDIC.B2.*f + g;
                flxDIC  = D.*(r.rDIC.A2.*dedz+r.rDIC.B2.*dfdz + dgdz);
            end
        end
        
    end
    
end


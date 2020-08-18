classdef benthic_zNO3
    % Solve NO3
    
    properties
        qdispNO3=308.42208;  	% NO3 diffusion coefficient in water (cm2/yr)
        adispNO3=12.2640;     	% NO3 linear coefficient for temperature dependence (cm2/yr/oC)
        DN1;                    % NO3 diffusion coefficient in bioturbated layer (cm2/yr)
        DN2;                  	% NO3 diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        KNH4=1.3;               % Adsorption coefficient (same in ocix and anoxic layer) (-)
    end
    
    methods
        function obj = benthic_zNO3(bsd, swi)
            obj.DN1=(obj.qdispNO3+obj.adispNO3*swi.T).*bsd.dispFactor+bsd.Dbio;
            obj.DN2=(obj.qdispNO3+obj.adispNO3*swi.T).*bsd.dispFactor;
            
            
        end
        
        function r = calc(obj, bsd, swi, r)
            % Iteratively solve for zno3
            
            % Try zero flux at zinf and see if we have any NO3 left - in
            % the rare case of zox < zinf but zNO3 = zinf, also
            % calculate [NO3] at zinf for advective loss
            [flxzinf, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r,2);
            
            
            if r.zox == bsd.zinf
                r.zno3 = bsd.zinf;
                bctype = 2;
            else
                
                fun=@(zno3)-obj.calcbc(zno3,bsd,swi,r,1);
                
                %                 % Try zero flux at zinf and see if we have any NO3 left - in
                %                 % the rare case of zox < zinf but zNO3 = zinf
                %                 [flxzinf, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r,2);
                
                if bsd.usescalarcode
                    if conczinf > 0     % Dom 30062016: change this as can be neg as well (NO3 different) -> check for zox == zinf
                        r.zno3 = bsd.zinf;
                        bctype = 2; % BC: zero flux
                    else
                        r.zno3=fzero(fun, [max(r.zox, 1e-10), bsd.zinf] ,bsd.fzerooptions);
                        bctype = 1; % BC: zero concentration
                        conczinf = 0.0;
                    end
                else  % same logic, in vector form
                    zno3=fzero_vec(fun,max(r.zox, 1e-10), bsd.zinf,bsd.fzerooptions);
                    bctype =         2*(conczinf>0) +     1*(conczinf<=0);
                    r.zno3 = bsd.zinf.*(conczinf>0) + zno3.*(conczinf<=0);
                end
                
            end
            [flxzno3, conczno3, flxswiNO3, r] = obj.calcbc(r.zno3, bsd, swi, r, bctype);
            
            flxswiNO3 = flxswiNO3 - bsd.por.*bsd.w.*(swi.NO30-conczinf);
            %             if(abs(flxswiNO3) <= bsd.tol_const)
            %                 flxswiNO3 = 0.0
            %             end
            
            r.flxzno3 = flxzno3;  % should be zero - not if zno3 = zinf
            r.conczno3 = conczno3; % may be non-zero if eg fully oxic sediment
            r.flxswiNO3 = flxswiNO3;
        end
        
        function [flxzno3, conczno3, flxswi,r] = calcbc(obj, zno3, bsd, swi, r, bctype)
            % Calculate trial solution for given zno3, matching boundary conditions from layer-by-layer solutions
            
            
            % Preparation: for each layer, sort out solution - matching across bioturbation boundary (if necessary)
            % layer 1: 0 < z < zox, nitrification
            %      ls =      prepfg_l12( bsd, swi, r, reac1,            ,ktemp, zU, zL, D1,        D2)
            rNO3.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, bsd.gamma.*bsd.NC1,  0,  0, r.zox, obj.DN1, obj.DN2);
            % layer 2: zox < z < zno3, denitrification
            rNO3.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, -bsd.NO3CR        ,  0,  r.zox, zno3, obj.DN1, obj.DN2);
            
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 2 zno3
            
            [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                = r.zTOC.calcfg_l12(zno3, bsd, swi, r,     -bsd.NO3CR   , 0, rNO3.ls2);
            
            
            % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 -> NO3 source)
            
            % basis functions at bottom of layer 1
            [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, bsd.gamma.*bsd.NC1, 0, rNO3.ls1);
            % basis functions at top of layer 2
            [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, -bsd.NO3CR        , 0, rNO3.ls2);
            
            %flux of NH4 to zox  TODO NH4 production by denitrification?
            FNH4 = r.zTOC.calcReac(zno3, bsd.zinf, bsd.NC1/(1+obj.KNH4) ,  bsd, swi, r); % MULTIPLY BY 1/POR ????
            % match solutions at zox - continuous concentration, flux discontinuity from H2S ox
            D = (r.zox <= bsd.zbio).*obj.DN1 + (r.zox > bsd.zbio).*obj.DN2;
            
            [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                0, -r.zxf.*bsd.gamma.*FNH4./D);
            
            % Solution at swi, top of layer 1
            [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                = r.zTOC.calcfg_l12(0, bsd, swi, r, bsd.gamma.*bsd.NC1, 0, rNO3.ls1);
            % transform to use coeffs from l2
            [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            
            
            % Solve for ANO3, BNO3 given boundary conditions (expressed in terms of transformed basis fns, layer 2 A, B)
            
            % Case 1 zero concentration at zno3
            % ANO3*e2_zno3   +  BNO3*f2_zno3  + g2_zno3 = 0;
            % ANO3*e1_0     +   BNO3*f1_0     + g1_0  = swi.NO30;
            
            % | e2_zno3 f2_zno3 |  |ANO3|   = | -g2_zno3       |
            % | e1_0     f1_0   |  |BNO3|     | swi.NO30 - g1_0 |
            
            [ bctype1_A2, bctype1_B2]      = benthic_utils.solve2eqn(e2_zno3, f2_zno3, e1_0, f1_0, -g2_zno3, swi.NO30 - g1_0);
            
            % Case  2 zero flux at zno3
            % ANO3*de2dz_zno3   +  BNO3*dfdz2_zno3  + dgdz2_zno3 = 0;
            % ANO3*e1_0         +   BNO3*f1_0       + g1_0       = swi.NO30;
            [ bctype2_A2, bctype2_B2]      = benthic_utils.solve2eqn(dedz2_zno3, dfdz2_zno3, e1_0, f1_0, -dgdz2_zno3, swi.NO30 - g1_0);
            
            % Choose type of solution requested (vectorized form)
            rNO3.A2 = (bctype==1).*bctype1_A2 + (bctype==2).*bctype2_A2;
            rNO3.B2 = (bctype==1).*bctype1_B2 + (bctype==2).*bctype2_B2;
            
            % calculate flux at zno3
            D = (zno3 <= bsd.zbio).*obj.DN1 + (zno3 > bsd.zbio).*obj.DN2;
            
            flxzno3 = D.*(rNO3.A2.*dedz2_zno3+rNO3.B2.*dfdz2_zno3 + dgdz2_zno3);        % includes 1/por ie flux per (cm^2 pore area)
            conczno3 = rNO3.A2.*e2_zno3+rNO3.B2.*f2_zno3 + g2_zno3;
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            flxswi = bsd.por.*(obj.DN1.*(rNO3.A2.*dedz1_0+rNO3.B2.*dfdz1_0 + dgdz1_0)); % - bsd.w.*swi.NO30);   % NB: use A2, B2 as these are _xformed_ layer 1 basis functions
            
            % save coeffs for layer 1
            rNO3.A1 = zox.a.*rNO3.A2 + zox.b.*rNO3.B2 + zox.e;
            rNO3.B1 = zox.c.*rNO3.A2 + zox.d.*rNO3.B2 + zox.f;
            
            
            r.rNO3 = rNO3;
            
        end
        
        
        function [NO3, flxNO3] = calcNO3(obj, z, bsd, swi, r)
            % Calculate NO3 concentration and flux at depth z from solution
            
            rNO3 = r.rNO3;
            if z <= r.zno3
                if z <= bsd.zbio
                    D = obj.DN1;
                else
                    D = obj.DN2;
                end
                
                if z <= r.zox   % layer 1
                    [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, bsd.gamma.*bsd.NC1 , 0, rNO3.ls1);
                    NO3     = r.rNO3.A1.*e + r.rNO3.B1.*f + g;
                    flxNO3  = D.*(r.rNO3.A1.*dedz+r.rNO3.B1.*dfdz + dgdz);
                elseif z <= r.zno3 % layer 2
                    [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, -bsd.NO3CR         , 0, rNO3.ls2);
                    NO3     = r.rNO3.A2.*e + r.rNO3.B2.*f + g;
                    flxNO3  = D.*(r.rNO3.A2.*dedz+r.rNO3.B2.*dfdz + dgdz);
                end
                
            else
                
                NO3 = 0;
                flxNO3 = 0;
            end
        end
        
    end
    
end


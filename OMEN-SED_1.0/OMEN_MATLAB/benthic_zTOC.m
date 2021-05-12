classdef benthic_zTOC < handle
    %Organic matter- 2 Fractions
    %
    % In bioturbated layer 1, each phase i (i=1,2) evolves according to:
    %  dCi/dt = DC * d^2 Ci/ dt^2 - w * dCi/dt - ki
    % in non-bioturbated layer 2,
    %  dCi/di = -w * dCi/dt - ki
    
    properties
        DC1;                    % TOC diffusion coefficient (cm2/yr)
        
        k1= 0.01;             	% TOC degradation rate constnat (1/yr)
        k2= 0.0001;             % TOC degradation rate constant (1/yr)
        
    end
    
    methods
        function obj = benthic_zTOC(bsd)
            % Constructior - set up constants etc
            
            obj.DC1 = bsd.Dbio;
            
        end
        
        function res = calc(obj, bsd, swi, res)
            % Calculate solution for TOC vs z
            
            rTOC.a11=(bsd.w-sqrt(bsd.w.^2+4.*obj.DC1.*obj.k1))./(2.*obj.DC1);
            rTOC.b11=(bsd.w+sqrt(bsd.w.^2+4.*obj.DC1.*obj.k1))./(2.*obj.DC1);
            rTOC.a21=(-obj.k1./bsd.w);
            % calculate bioturbated SWI
            % comment C01 calculation for calculating observed profiles &
            % sensitivity analysis (as we have concentrations C0i for this)
%             swi.C01 = (swi.Fnonbio1*(-rTOC.a11*exp(rTOC.a11*bsd.zbio)+rTOC.b11*exp(rTOC.b11*bsd.zbio)))/(-obj.DC1*rTOC.b11*rTOC.a11*exp(rTOC.b11*bsd.zbio) + obj.DC1*rTOC.b11*rTOC.a11*exp(rTOC.a11*bsd.zbio) + ...
%                 obj.DC1*rTOC.b11*rTOC.a11*bsd.por*exp(rTOC.b11*bsd.zbio) - obj.DC1*rTOC.b11*rTOC.a11*bsd.por*exp(rTOC.a11*bsd.zbio) - bsd.w*rTOC.a11*exp(rTOC.a11*bsd.zbio) + bsd.w*rTOC.b11*exp(rTOC.b11*bsd.zbio) + ...
%                 bsd.w*bsd.por*rTOC.a11*exp(rTOC.a11*bsd.zbio) - bsd.w*bsd.por*rTOC.b11*exp(rTOC.b11*bsd.zbio));
            res.swi.C01 = swi.C01;
            
            rTOC.A11=-(swi.C01.*rTOC.b11.*exp(rTOC.b11.*bsd.zbio))./(rTOC.a11.*exp(rTOC.a11.*bsd.zbio)-rTOC.b11.*exp(rTOC.b11.*bsd.zbio)+bsd.tol_const);
            rTOC.A21=(rTOC.A11.*(exp(rTOC.a11.*bsd.zbio)-exp(rTOC.b11.*bsd.zbio))+swi.C01.*exp(rTOC.b11.*bsd.zbio))./(exp(rTOC.a21.*bsd.zbio)+bsd.tol_const);
            
            rTOC.a12=(bsd.w-sqrt(bsd.w.^2+4.*obj.DC1.*obj.k2))./(2.*obj.DC1);
            rTOC.b12=(bsd.w+sqrt(bsd.w.^2+4.*obj.DC1.*obj.k2))./(2.*obj.DC1);
            rTOC.a22=(-obj.k2./bsd.w);
            % calculate bioturbated SWI
            % comment C02 calculation for sensitivity analysis (as we have concentrations C0i for this)
%             swi.C02 = (swi.Fnonbio2*(-rTOC.a12*exp(rTOC.a12*bsd.zbio)+rTOC.b12*exp(rTOC.b12*bsd.zbio)))/(-obj.DC1*rTOC.b12*rTOC.a12*exp(rTOC.b12*bsd.zbio) + obj.DC1*rTOC.b12*rTOC.a12*exp(rTOC.a12*bsd.zbio) + ...
%                 obj.DC1*rTOC.b12*rTOC.a12*bsd.por*exp(rTOC.b12*bsd.zbio) - obj.DC1*rTOC.b12*rTOC.a12*bsd.por*exp(rTOC.a12*bsd.zbio) - bsd.w*rTOC.a12*exp(rTOC.a12*bsd.zbio) + bsd.w*rTOC.b12*exp(rTOC.b12*bsd.zbio) + ...
%                 bsd.w*bsd.por*rTOC.a12*exp(rTOC.a12*bsd.zbio) - bsd.w*bsd.por*rTOC.b12*exp(rTOC.b12*bsd.zbio));
            res.swi.C02 = swi.C02;
            
            rTOC.A12=-(swi.C02.*rTOC.b12.*exp(rTOC.b12.*bsd.zbio))./(rTOC.a12.*exp(rTOC.a12.*bsd.zbio)-rTOC.b12.*exp(rTOC.b12.*bsd.zbio)+bsd.tol_const);
            rTOC.A22=(rTOC.A12.*(exp(rTOC.a12.*bsd.zbio)-exp(rTOC.b12.*bsd.zbio))+swi.C02.*exp(rTOC.b12.*bsd.zbio))./(exp(rTOC.a22.*bsd.zbio)+bsd.tol_const);
            
            res.rTOC = rTOC;
            
            % Sed input flux to upper boundary, per cm^2 water column
            [res.Fswi_TOC, res.Fswi_TOC1, res.Fswi_TOC2] = obj.calcCflx(bsd.z0, bsd, swi, res);
            
            %Flux through lower boundary zinf, per cm^2 water-column
            res.F_TOC1=-(1-bsd.por).*bsd.w.*rTOC.A21.*exp(rTOC.a21.*bsd.zinf);
            res.F_TOC2=-(1-bsd.por).*bsd.w.*rTOC.A22.*exp(rTOC.a22.*bsd.zinf);
            res.F_TOC=res.F_TOC1+res.F_TOC2;
            
            if(false)   % prints for checking TOC integration for SWI-flux test
                Int_nonbio = obj.k1*((rTOC.A21/rTOC.a21*exp(rTOC.a21*bsd.zinf)) - rTOC.A21/rTOC.a21*exp(rTOC.a21*bsd.zbio))
                Int_bio_Maple = obj.k1*( (rTOC.A11*(exp(rTOC.a11*bsd.zbio)/rTOC.a11 - exp(rTOC.b11*bsd.zbio)/rTOC.b11)) + swi.C01/rTOC.b11 * exp(rTOC.b11*bsd.zbio) - (rTOC.A11*(1/rTOC.a11 - 1/rTOC.b11)+ swi.C01/rTOC.b11)) % for C(inf) = 0:   k* (A1/a1*exp(a1*zinf)-A1/a1)
                
                Int_bio_z0 = obj.k1*((rTOC.A11*(1/rTOC.a11 - 1/rTOC.b11)+ swi.C01/rTOC.b11))
                Int_bio_zbio= obj.k1*((rTOC.A11*(exp(rTOC.a11*bsd.zbio)/rTOC.a11 - exp(rTOC.b11*bsd.zbio)/rTOC.b11)) + swi.C01/rTOC.b11 * exp(rTOC.b11*bsd.zbio))
                Int_nonbio_zbio=obj.k1*(rTOC.A21/rTOC.a21*exp(rTOC.a21*bsd.zbio))
                Int_nonbio_zinf= obj.k1*((rTOC.A21/rTOC.a21*exp(rTOC.a21*bsd.zinf)))
                
                Sum_bio_nonbio_ALL = Int_nonbio + Int_bio_Maple
                
                for i=1:1:1001
                    z(i) = (i-1)/10;
                    if(z<=bsd.zbio)
                        C1(i)=rTOC.A11*(exp(rTOC.a11*z(i))-exp(rTOC.b11*z(i)))+swi.C01*exp(rTOC.b11*z(i));
                    else
                        C1(i)=rTOC.A21*exp(rTOC.a21*z(i));
                    end
                end
                %                 figure
                %                 plot(C1,-z)
                %                  hold on
                %                 t=xlim;
                %                 plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                %                 xlabel ('[TOC] in zTOC')
                %                 ylabel('Depth (cm)')
                %                 hold off
            end
            
        end
        
        
        function [C, C1, C2] = calcC(obj, z, bsd, swi, r)
            % Calculate concentration at depth z from solution
            rTOC = r.rTOC;
            if(z<=bsd.zbio)
                C1=rTOC.A11.*(exp(rTOC.a11.*z)-exp(rTOC.b11.*z))+swi.C01.*exp(rTOC.b11.*z);
                C2=rTOC.A12.*(exp(rTOC.a12.*z)-exp(rTOC.b12.*z))+swi.C02.*exp(rTOC.b12.*z);
            else
                C1=rTOC.A21.*exp(rTOC.a21.*z);
                C2=rTOC.A22.*exp(rTOC.a22.*z);
            end
            
            C = C1 + C2;
        end
        
        function [Cflx, C1flx, C2flx] = calcCflx(obj, z, bsd, swi, r)
            % Calculate flux at depth z, from solution
            rTOC = r.rTOC;
            [C, C1, C2] = obj.calcC(z,bsd,swi,r);
            if(z<bsd.zbio)
                dC1dz =  rTOC.A11.*(rTOC.a11.*exp(rTOC.a11.*z)-rTOC.b11.*exp(rTOC.b11.*z))+swi.C01.*rTOC.b11.*exp(rTOC.b11.*z);
                C1flx = - (1-bsd.por).*(-obj.DC1.*dC1dz + bsd.w.*C1);
                dC2dz =  rTOC.A12.*(rTOC.a12.*exp(rTOC.a12.*z)-rTOC.b12.*exp(rTOC.b12.*z))+swi.C02.*rTOC.b12.*exp(rTOC.b12.*z);
                C2flx = - (1-bsd.por).*(-obj.DC1.*dC2dz + bsd.w.*C2);
            else
                C1flx = - (1-bsd.por).*bsd.w.*C1;
                C2flx = - (1-bsd.por).*bsd.w.*C2;
            end
            
            Cflx = C1flx + C2flx;
        end
        
        function [e, dedz, f, dfdz, g, dgdz] ...
                = calcfg_l1(obj, z, bsd, swi, res, reac1, reac2, Dtemp, ktemp)
            % Basis functions for solutes, case z <= zbio
            %
            % reac1, reac2        - mol./mol S released per organic carbon C
            %
            % General solution for solute S is given by
            %  S(z) = A .* e(z) + B .* f(z) + g(z)
            
            
            r = res.rTOC;
            
            e = ones(1,bsd.ncl);
            dedz = zeros(1,bsd.ncl);
            
            b1=bsd.w./Dtemp;
            f=exp(z.*b1);
            dfdz = b1.*exp(z.*b1);
            
            %pfac=1./bsd.por;   % assume org matter already .*(1-bsd.por)
            pfac = 1;          % in fact, already has (1-por)/por
            
            PhiI1   =-pfac.*obj.k1.*(reac1).*r.A11./(Dtemp.*r.a11.^2-bsd.w.*r.a11-ktemp);
            PhiII1  = pfac.*obj.k1.*(reac1).*r.A11./(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp);
            PhiIII1 =-pfac.*obj.k1.*(reac1).*swi.C01./(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp);
            PhiI2   =-pfac.*obj.k2.*(reac2).*r.A12./(Dtemp.*r.a12.^2-bsd.w.*r.a12-ktemp);
            PhiII2  = pfac.*obj.k2.*(reac2).*r.A12./(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp);
            PhiIII2 =-pfac.*obj.k2.*(reac2).*swi.C02./(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp);
            
            ea11z = exp(r.a11.*z);
            eb11z = exp(r.b11.*z);
            ea12z = exp(r.a12.*z);
            eb12z = exp(r.b12.*z);
            
            g =  PhiI1.*ea11z + PhiII1.*eb11z + PhiIII1.*eb11z + ...
                PhiI2.*ea12z + PhiII2.*eb12z + PhiIII2.*eb12z;
            
            
            dgdz = PhiI1.*r.a11.*ea11z + PhiII1.*r.b11.*eb11z + PhiIII1.*r.b11.*eb11z + ...
                PhiI2.*r.a12.*ea12z + PhiII2.*r.b12.*eb12z + PhiIII2.*r.b12.*eb12z;
            
            
        end
        function [ e, dedz, f, dfdz, g, dgdz] ...
                = calcfg_l2(obj, z, bsd, swi, res, reac1, reac2, Dtemp, ktemp)
            % Basis functions for solutes, case z > zbio
            % reac1, reac2        - mol/mol S released per organic carbon C
            %
            % General solution for solute S is given by
            %  S(z) = A .* e(z) + B .* f(z) + g(z)
            
            r = res.rTOC;
            
            e = ones(1,bsd.ncl);
            dedz = zeros(1,bsd.ncl);
            
            b2 = bsd.w./Dtemp;
            f=exp(z.*b2);
            dfdz = b2.*exp(z.*b2);
            
            %pfac=1./bsd.por;   % assume org matter already .*(1-bsd.por)
            pfac = 1;          % in fact, already has (1-por)/por
            
            PhiI1 = -pfac.*obj.k1.*(reac1).*r.A21./(Dtemp.*r.a21.^2-bsd.w.*r.a21-ktemp);
            PhiI2 = -pfac.*obj.k2.*(reac2).*r.A22./(Dtemp.*r.a22.^2-bsd.w.*r.a22-ktemp);
            
            g = PhiI1.*exp(r.a21.*z) + ...
                PhiI2.*exp(r.a22.*z);
            dgdz = PhiI1.*r.a21.*exp(r.a21.*z) + ...
                PhiI2.*r.a22.*exp(r.a22.*z);
            
            
        end
        
        
        
        function [ e, dedz, f, dfdz, g, dgdz] ...
                = calcfg_l12(obj, z, bsd, swi, res, reac1, reac2, ktemp, ls)
            % calculate solution basis functions, for layer which may cross bioturbation boundary
            %
            % reac1, reac2        - mol/mol S released per organic carbon C
            %
            % General solution for solute S is given by
            %  S(z) = A .* e(z) + B .* f(z) + g(z)
            %
            % Where e,f,g are generated by matching solutions across bioturbation boundary (if necessary)
            % Solution properties (matching etc) are input in ls
            % On input, ls should contain fields generated by prepfg_l12
            
            if bsd.usescalarcode
                switch ls.ltype
                    case 1  % bioturbated
                        %                        msg='CASE 1'
                        [e, dedz, f, dfdz, g, dgdz] ...
                            = obj.calcfg_l1( z, bsd, swi, res, reac1, reac2, ls.D1, ktemp);
                    case 2 % not bioturbated
                        %                        msg='CASE 2'
                        [ e, dedz, f, dfdz, g, dgdz] ...
                            = obj.calcfg_l2( z, bsd, swi, res, reac1, reac2, ls.D2, ktemp);
                    case 3 % crossing boundary
                        if z >= bsd.zbio  % below bioturbated regio
                            %                            msg = 'crossing boundary - IFFFFFF'
                            [e, dedz, f, dfdz, g, dgdz] ...
                                = obj.calcfg_l2( z, bsd, swi, res, reac1, reac2, ls.D2, ktemp);
                        else   % above bioturbated region
                            %                            msg = 'crossing boundary - ELSEEE'
                            [e_1, dedz_1, f_1, dfdz_1, g_1, dgdz_1] ...
                                = obj.calcfg_l1( z, bsd, swi, res, reac1, reac2, ls.D1, ktemp);
                            
                            % Now find 'transformed' basis functions such that in layer 1, O2 = A_2*et + B_2*ft + gt  (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
                            [e, f, g, dedz, dfdz, dgdz] = benthic_utils.xformsoln(e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1, ...
                                ls.a , ls.b , ls.c , ls.d , ls.e , ls.f);
                            
                        end
                    otherwise
                        error('unrecognized ltype  %g\n',ls.ltype);
                end
            else % same logic as scalar code, written in vectorised form
                [e_1, dedz_1, f_1, dfdz_1, g_1, dgdz_1] ...
                    = obj.calcfg_l1( z, bsd, swi, res, reac1, reac2, ls.D1, ktemp);
                % Find 'transformed' basis functions (ie basis at z for A,B at zL: transformation = identity if layer fully bioturbated)
                [e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1] = benthic_utils.xformsoln(e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1, ...
                    ls.a , ls.b , ls.c , ls.d , ls.e , ls.f);
                tl = z<=bsd.zbio;  % bioturbated
                e(tl)=e_1(tl); dedz(tl)=dedz_1(tl); f(tl) = f_1(tl); dfdz(tl) = dfdz_1(tl); g(tl)=g_1(tl); dgdz(tl) = dgdz_1(tl);
                
                [e_2, dedz_2, f_2, dfdz_2, g_2, dgdz_2] ...
                    = obj.calcfg_l2( z, bsd, swi, res, reac1, reac2, ls.D2, ktemp);
                tl = z>bsd.zbio; % not bioturbated
                e(tl)=e_2(tl); dedz(tl)=dedz_2(tl); f(tl) = f_2(tl); dfdz(tl) = dfdz_2(tl); g(tl)=g_2(tl); dgdz(tl) = dgdz_2(tl);
            end
        end
        
        function ls = prepfg_l12(obj, bsd, swi, res, reac1, reac2, ktemp, zU, zL, D1, D2)
            % calculate solution matching (if necessary) for layer that may cross a bioturbation boundary
            ls.zU = zU;
            ls.zL = zL;
            ls.D1 = D1;
            ls.D2 = D2;
            
            if bsd.usescalarcode
                if zL <= bsd.zbio  % wholly within bioturbated layer
                    ls.ltype = 1;
                elseif zU >= bsd.zbio % wholly within non-bioturbated layer
                    ls.ltype = 2;
                else             % crossing boundary - sort out solution matching at zbio
                    ls.ltype = 3;
                    
                    [e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1] ...
                        = obj.calcfg_l1(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D1, ktemp);
                    [e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2] ...
                        = obj.calcfg_l2(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D2, ktemp);
                    
                    % match solutions at zbio - continuous concentration and flux
                    [ls.a, ls.b, ls.c, ls.d, ls.e, ls.f] = benthic_utils.matchsoln( ...
                        e_zbio_l1, f_zbio_l1, g_zbio_l1, ls.D1*dedz_zbio_l1, ls.D1*dfdz_zbio_l1, ls.D1*dgdz_zbio_l1, ...
                        e_zbio_l2, f_zbio_l2, g_zbio_l2, ls.D2*dedz_zbio_l2, ls.D2*dfdz_zbio_l2, ls.D2*dgdz_zbio_l2, ...
                        0, 0);
                end
            else
                n = bsd.ncl;
                % initialise to identity transformation
                ls.a = ones(1,n); ls.b = zeros(1,n); ls.c=zeros(1,n); ls.d = ones(1,n); ls.e = zeros(1,n); ls.f = zeros(1,n);
                
                % wholly within bioturbated layer
                ls.ltype(zL <= bsd.zbio) = 1;
                
                % wholly within non-bioturbated layer
                ls.ltype(zU >= bsd.zbio) = 2;
                
                % crossing boundary - sort out solution matching at zbio
                tl = zL > bsd.zbio & zU < bsd.zbio;
                ls.ltype(tl) = 3;
                
                [e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1] ...
                    = obj.calcfg_l1(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D1, ktemp);
                [e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2] ...
                    = obj.calcfg_l2(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D2, ktemp);
                
                % match solutions at zbio - continuous concentration and flux
                [a, b, c, d, e, f] = benthic_utils.matchsoln( ...
                    e_zbio_l1, f_zbio_l1, g_zbio_l1, ls.D1.*dedz_zbio_l1, ls.D1.*dfdz_zbio_l1, ls.D1.*dgdz_zbio_l1, ...
                    e_zbio_l2, f_zbio_l2, g_zbio_l2, ls.D2.*dedz_zbio_l2, ls.D2.*dfdz_zbio_l2, ls.D2.*dgdz_zbio_l2, ...
                    0, 0);
                
                ls.a(tl) = a(tl); ls.b(tl)=b(tl); ls.c(tl)=c(tl); ls.d(tl)=d(tl); ls.e(tl)=e(tl); ls.f(tl)=f(tl);
            end
            
        end
        
        
        function FReac = calcReac(obj, zU, zL, reac1, reac2, bsd, swi, res)
            % Integral of reacted organic matter from zU to zL,
            % multiplied by stoichiometric factors reac1, reac2 (for the two OC phases)
            
            % Vector-friendly way of handling 3 cases:
            % 1) wholly within bioturbated layer:    calcReac_l1(zU,zL)     + (0 =) calcReac_l2(bsd.zbio, bsd.zbio)
            % 2) wholly within non-bio     layer:  (0=) calcReac_l1(zbio, zbio) +   calcReac_l2(zU, zL)
            % 3) crossing zbio                       calcRead_l1(zU,zbio)   +       calcReac_l2(zbio, zL)
            
            FReac=obj.calcReac_l1(min(zU,bsd.zbio), min(zL,bsd.zbio), reac1, reac2, bsd, swi, res) ...
                + obj.calcReac_l2(max(zU,bsd.zbio), max(zL, bsd.zbio), reac1, reac2, bsd, swi, res);
            
            
            % TODO confirm (1-bsd.por).*  has been added (to k1 & k2 ?)
        end
        
        function FReac1 = calcReac_l1(obj, zU, zL, reac1, reac2, bsd, swi, res)
            % Integral of reacted organic matter, zU and zL within layer 1 (bioturbated)
            
            r = res.rTOC;
            
            reacf1=res.zTOC.k1.*reac1;
            reacf2=res.zTOC.k2.*reac2;
            
            FReac1 =-reacf1.*(r.A11.*(exp(r.a11.*zU).*r.b11 - exp(r.b11.*zU).*r.a11 - exp(r.a11.*zL).*r.b11 + exp(r.b11.*zL).*r.a11)...
                +swi.C01.*exp(r.b11.*zU).*r.a11 -swi.C01.*exp(r.b11.*zL).*r.a11)./(r.a11.*r.b11) ...
                -reacf2.*(r.A12.*(exp(r.a12.*zU).*r.b12 - exp(r.b12.*zU).*r.a12 - exp(r.a12.*zL).*r.b12 + exp(r.b12.*zL).*r.a12) ...
                +swi.C02.*exp(r.b12.*zU).*r.a12 -swi.C02.*exp(r.b12.*zL).*r.a12)./(r.a12.*r.b12);
        end
        
        function FReac2 = calcReac_l2(obj, zU, zL, reac1, reac2, bsd, swi, res)
            % Integral of reacted organic matter, zU and zL within layer 2 (non bioturbated)
            
            r = res.rTOC;
            
            reacf1=res.zTOC.k1.*reac1;
            reacf2=res.zTOC.k2.*reac2;
            
            FReac2= -reacf1.*r.A21.*(exp(r.a21.*zU) - exp(r.a21.*zL))./r.a21...
                -reacf2.*r.A22.*(exp(r.a22.*zU) - exp(r.a22.*zL))./r.a22;
        end
        
        
        function FOM = calcOM(obj, zU, zL, reac1, reac2, bsd, swi, res)
            % Integral of organic matter from zU to zL,
            % multiplied by stoichiometric factors reac1, reac2 (for the two OC phases)
            
            % Vector-friendly way of handling 3 cases:
            % 1) wholly within bioturbated layer:    calcReac_l1(zU,zL)     + (0 =) calcReac_l2(bsd.zbio, bsd.zbio)
            % 2) wholly within non-bio     layer:  (0=) calcReac_l1(zbio, zbio) +   calcReac_l2(zU, zL)
            % 3) crossing zbio                       calcRead_l1(zU,zbio)   +       calcReac_l2(zbio, zL)
            
            FOM=obj.calcOM_l1(min(zU,bsd.zbio), min(zL,bsd.zbio), reac1, reac2, bsd, swi, res) ...
                + obj.calcOM_l2(max(zU,bsd.zbio), max(zL, bsd.zbio), reac1, reac2, bsd, swi, res);
            
            
            % TODO confirm (1-bsd.por).*  has been added (to k1 & k2 ?)
        end
        
        function FOM1 = calcOM_l1(obj, zU, zL, reac1, reac2, bsd, swi, res)
            % Integral of reacted organic matter, zU and zL within layer 1 (bioturbated)
            
            r = res.rTOC;
            
            reacf1=reac1;
            reacf2=reac2;
            
            FOM1 =-reacf1.*(r.A11.*(exp(r.a11.*zU).*r.b11 - exp(r.b11.*zU).*r.a11 - exp(r.a11.*zL).*r.b11 + exp(r.b11.*zL).*r.a11)...
                +swi.C01.*exp(r.b11.*zU).*r.a11 -swi.C01.*exp(r.b11.*zL).*r.a11)./(r.a11.*r.b11) ...
                -reacf2.*(r.A12.*(exp(r.a12.*zU).*r.b12 - exp(r.b12.*zU).*r.a12 - exp(r.a12.*zL).*r.b12 + exp(r.b12.*zL).*r.a12) ...
                +swi.C02.*exp(r.b12.*zU).*r.a12 -swi.C02.*exp(r.b12.*zL).*r.a12)./(r.a12.*r.b12);
        end
        
        function FOM2 = calcOM_l2(obj, zU, zL, reac1, reac2, bsd, swi, res)
            % Integral of reacted organic matter, zU and zL within layer 2 (non bioturbated)
            
            r = res.rTOC;
            
            reacf1=reac1;
            reacf2=reac2;
            
            FOM2= -reacf1.*r.A21.*(exp(r.a21.*zU) - exp(r.a21.*zL))./r.a21...
                -reacf2.*r.A22.*(exp(r.a22.*zU) - exp(r.a22.*zL))./r.a22;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %    FCTs FOR PO4 & M
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P,...
                e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M] ...
                = calcfg_l12_PO4_M(obj, z, bsd, swi, res, reac1P, reac2P, ktempP, QtempP, alphaP, ls, ktempM, QtempM, alphaM)
            % calculate solution basis functions, for layer which may cross bioturbation boundary
            %
            % reac1, reac2        - mol/mol S released per organic carbon C
            %
            % General solution for solute S is given by
            %  S(z) = A .* e(z) + B .* f(z) + g(z)
            %
            % Where e,f,g are generated by matching solutions across bioturbation boundary (if necessary)
            % Solution properties (matching etc) are input in ls
            % On input, ls should contain fields generated by prepfg_l12
            
            if bsd.usescalarcode
                switch ls.ltype
                    case 1  % bioturbated
                        if(alphaP==0)  % was z<=res.zox oxic layer -> call PO4 first
                            [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, a1_P, b1_P, Phi1_P] ...
                                = obj.calcfg_l1_PO4( z, bsd, swi, res, reac1P, reac2P, ls.D1P, ktempP, QtempP, 0,0, alphaP);
                            [e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, a1_M, b1_M] ...
                                = obj.calcfg_l1_M( z, bsd, swi, res, ls.D1M, ktempM, QtempM, a1_P, b1_P, Phi1_P, alphaM);
                        else    % anoxic layer -> call M first
                            [e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, a1_M, b1_M] ...
                                = obj.calcfg_l1_M(z, bsd, swi, res, ls.D1M, ktempM, QtempM, 0, 0, 0, alphaM);
                            [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, a1_P, b1_P, Phi1_P] ...
                                = obj.calcfg_l1_PO4(z, bsd, swi, res, reac1P, reac2P, ls.D1P, ktempP, QtempP, a1_M, b1_M, alphaP);
                        end
                        
                    case 2 % not bioturbated
                        if(alphaP==0)  % was z<=res.zox oxic layer -> call PO4 first
                            [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, a2_P, b2_P, Phi2_P] ...
                                = obj.calcfg_l2_PO4(z, bsd, swi, res, reac1P, reac2P, ls.D2P, ktempP, QtempP, 0, alphaP);
                            [e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, a2_M] ...
                                = obj.calcfg_l2_M( z, bsd, swi, res, ktempM, QtempM, a2_P, b2_P, Phi2_P, alphaM);
                        else    % anoxic layer -> call M first
                            [e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, a2_M] ...
                                = obj.calcfg_l2_M( z, bsd, swi, res, ktempM, QtempM, 0, 0, 0, alphaM);
                            [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, a2_P, b2_P, Phi2_P] ...
                                = obj.calcfg_l2_PO4(z, bsd, swi, res, reac1P, reac2P, ls.D2P, ktempP, QtempP, a2_M, alphaP);
                        end
                        
                    case 3 % crossing boundary
                        if z > bsd.zbio  % not bioturbated region
                            if(alphaP==0) % was z<=res.zox oxic layer -> call PO4 first BUT DECIDE VIA ALPHA_M NOT WITH <+ ZOX!!! DOESN't WORK FOR BOUNDARY ZOX
                                %%% DOMINIK: CASE 2: LAYER 2: have 3 int. const.
                                [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, a2_P, b2_P, Phi2_P] ...
                                    = obj.calcfg_l2_PO4(z, bsd, swi, res, reac1P, reac2P, ls.D2P, ktempP, QtempP, 0, alphaP);
                                [e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, a2_M] ...
                                    = obj.calcfg_l2_M( z, bsd, swi, res, ktempM, QtempM, a2_P, b2_P, Phi2_P, alphaM);
                            else    % anoxic layer -> call M first
                                %%% DOMINIK: CASE 1&2: LAYER 3: have 3 int. const.
                                [e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, a2_M] ...
                                    = obj.calcfg_l2_M( z, bsd, swi, res, ktempM, QtempM, 0, 0, 0, alphaM);
                                [e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, a2_P, b2_P, Phi2_P] ...
                                    = obj.calcfg_l2_PO4(z, bsd, swi, res, reac1P, reac2P, ls.D2P, ktempP, QtempP, a2_M, alphaP);
                            end
                        else   % bioturbated region z <= zbio
                            if(alphaP==0)  % was z<=res.zox oxic layer -> call PO4 first
                                %%% DOMINIK: CASE 1 & 2: LAYER 1: have 4 int. const.
                                [e_P_1, dedz_P_1, f_P_1, dfdz_P_1, g_P_1, dgdz_P_1, p_P_1, dpdz_P_1, q_P_1, dqdz_P_1, a1_P, b1_P, Phi1_P] ...
                                    = obj.calcfg_l1_PO4( z, bsd, swi, res, reac1P, reac2P, ls.D1P, ktempP, QtempP, 0,0, alphaP);
                                [e_M_1, dedz_M_1, f_M_1, dfdz_M_1, g_M_1, dgdz_M_1, p_M_1, dpdz_M_1, q_M_1, dqdz_M_1, a1_M, b1_M] ...
                                    = obj.calcfg_l1_M( z, bsd, swi, res, ls.D1M, ktempM, QtempM, a1_P, b1_P, Phi1_P, alphaM);
                                %%% DOMINIK: FOR CASE 2: DON'T HAVE D FROM LAYER 2
                            else    % anoxic layer -> call M first
                                %%% DOMINIK: CASE 1: LAYER 2: have 4 int. const.
                                [e_M_1, dedz_M_1, f_M_1, dfdz_M_1, g_M_1, dgdz_M_1, p_M_1, dpdz_M_1, q_M_1, dqdz_M_1, a1_M, b1_M] ...
                                    = obj.calcfg_l1_M(z, bsd, swi, res, ls.D1M, ktempM, QtempM, 0, 0, 0, alphaM);
                                [e_P_1, dedz_P_1, f_P_1, dfdz_P_1, g_P_1, dgdz_P_1, p_P_1, dpdz_P_1, q_P_1, dqdz_P_1, a1_P, b1_P, Phi1_P] ...
                                    = obj.calcfg_l1_PO4(z, bsd, swi, res, reac1P, reac2P, ls.D1P, ktempP, QtempP, a1_M, b1_M, alphaP);
                            end
                            
                            % Now find 'transformed' basis functions such that in layer 1,
                            % O2 = A_2*et + B_2*ft + gt  (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
                            % DOMINIK: TODO: don't clculate B (or rather F) like this:
                            EFPQ_P = [e_P_1; f_P_1; p_P_1; q_P_1];
                            dEFPQdz_P = [dedz_P_1; dfdz_P_1; dpdz_P_1; dqdz_P_1];
                            EFPQ_M = [e_M_1; f_M_1; p_M_1; q_M_1];
                            dEFPQdz_M = [dedz_M_1; dfdz_M_1; dpdz_M_1; dqdz_M_1];
                            [EFPQ_P_t, g_P, dEFPQ_P_t, dgdz_P, EFPQ_M_t, g_M, dEFPQ_M_t, dgdz_M] ...
                                = benthic_utils.xformsoln_PO4_M(EFPQ_P, EFPQ_M, dEFPQdz_P, dEFPQdz_M, g_P_1, g_M_1,dgdz_P_1, dgdz_M_1, ls.C, ls.D);
                            
                            % CHECK/TODO: OR SHALL I RATHER LEAVE THE VALUES IN MATRICES???? OR
                            % SAVE EACH IN SINGLE VARIABLE??? ... CASES BEFORE HAVE THEM IN
                            % ORDINARY VARIABLES!!!
                            
                            %%% WHEN lTYPE=3
                            %%% DEAL WITH ONE VARIABLE SHORT FROM LAYER BELOW
                            % FOR CASE 1: no Q from layer below
                            % FOR CASE2: no Q from layer below
                            if(res.zox<=bsd.zbio)   % CASE 1: no F from layer below - no e_M & f_M as well!!!??? but 0 anyway at the moment
                                e_P = EFPQ_P_t(1);
                                f_P = EFPQ_P_t(2); % was = EFPQ_P_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                p_P = EFPQ_P_t(3);
                                q_P = q_P_1;
                                
                                dedz_P = dEFPQ_P_t(1);
                                dfdz_P = dEFPQ_P_t(2);  % was = dEFPQ_P_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                dpdz_P = dEFPQ_P_t(3);
                                dqdz_P = dqdz_P_1;
                                
                                e_M = EFPQ_M_t(1);
                                f_M = EFPQ_M_t(2);  % was = EFPQ_M_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                p_M = EFPQ_M_t(3);
                                q_M = q_M_1;
                                
                                dedz_M = dEFPQ_M_t(1);
                                dfdz_M = dEFPQ_M_t(2);  % was = dEFPQ_M_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                dpdz_M = dEFPQ_M_t(3);
                                dqdz_M = dqdz_M_1;
                            else    % CASE 2: no Q from layer below - no p_P as well!!!!!
                                e_P = EFPQ_P_t(1);
                                f_P = EFPQ_P_t(2); % was = EFPQ_P_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                p_P = p_P_1; % was = EFPQ_P_t(3);
                                q_P = q_P_1;
                                
                                dedz_P = dEFPQ_P_t(1);
                                dfdz_P = dEFPQ_P_t(2);  % was = dEFPQ_P_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                dpdz_P = dpdz_P_1; % was = dEFPQ_P_t(3);
                                dqdz_P = dqdz_P_1;
                                
                                e_M = EFPQ_M_t(1);
                                f_M = EFPQ_M_t(2);  % was = EFPQ_M_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                p_M = EFPQ_M_t(3);
                                q_M = q_M_1;
                                
                                dedz_M = dEFPQ_M_t(1);
                                dfdz_M = dEFPQ_M_t(2);  % was = dEFPQ_M_t(2) BUT IS ZERO IN LAYER BELOW, SO USE VALUE FROM THIS LAYER
                                dpdz_M = dEFPQ_M_t(3);
                                dqdz_M = dqdz_M_1;
                            end
                            
                            
                        end
                        
                        
                        
                        
                    otherwise
                        error('unrecognized ltype  %g\n',ls.ltype);
                end  % switch
            else % same logic as scalar code, written in vectorised form
                
            end
        end
        
        
        
        function ls = prepfg_l12_PO4_M(obj, bsd, swi, res, reac1, reac2, ktempP, QtempP, zU, zL, D1P, D2P, alphaP, ...
                ktempM, QtempM, D1M,   D2M, alphaM)
            % calculate solution matching (if necessary) for layer that may cross a bioturbation boundary
            ls.zU = zU;
            ls.zL = zL;
            ls.D1P = D1P;
            ls.D2P = D2P;
            ls.D1M = D1M;
            ls.D2M = D2M;
            
            if bsd.usescalarcode
                if zL <= bsd.zbio  % wholly within bioturbated layer
                    ls.ltype = 1;
                elseif zU >= bsd.zbio % wholly within non-bioturbated layer
                    ls.ltype = 2;
                else             % crossing boundary - sort out solution matching at zbio
                    ls.ltype = 3;
                    
                    if(zL<=res.zox)   % oxic layer -> call PO4 first  DOMINIK: TODO: MAYBE HERE ALSO ask alphaP~=0!?
                        [e_zbio_l1_P, dedz_zbio_l1_P, f_zbio_l1_P, dfdz_zbio_l1_P, g_zbio_l1_P, ...
                            dgdz_zbio_l1_P, p_zbio_l1_P, dpdz_zbio_l1_P, q_zbio_l1_P, dqdz_zbio_l1_P, a1_P, b1_P, Phi1_P] ...
                            = obj.calcfg_l1_PO4(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D1P, ktempP, QtempP, 0,0, alphaP);
                        
                        [e_zbio_l2_P, dedz_zbio_l2_P, f_zbio_l2_P, dfdz_zbio_l2_P, g_zbio_l2_P, ...
                            dgdz_zbio_l2_P, p_zbio_l2_P, dpdz_zbio_l2_P, q_zbio_l2_P, dqdz_zbio_l2_P, a2_P, b2_P, Phi2_P] ...
                            = obj.calcfg_l2_PO4(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D2P, ktempP, QtempP, 0, alphaP);
                        
                        [e_zbio_l1_M, dedz_zbio_l1_M, f_zbio_l1_M, dfdz_zbio_l1_M, g_zbio_l1_M, dgdz_zbio_l1_M, ...
                            p_zbio_l1_M, dpdz_zbio_l1_M, q_zbio_l1_M, dqdz_zbio_l1_M, a1_M, b1_M] ...
                            = obj.calcfg_l1_M(bsd.zbio, bsd, swi, res, ls.D1M, ktempM, QtempM, a1_P, b1_P, Phi1_P, alphaM);
                        
                        [e_zbio_l2_M, dedz_zbio_l2_M, f_zbio_l2_M, dfdz_zbio_l2_M, g_zbio_l2_M, dgdz_zbio_l2_M, ...
                            p_zbio_l2_M, dpdz_zbio_l2_M, q_zbio_l2_M, dqdz_zbio_l2_M, a2_M] ...
                            = obj.calcfg_l2_M(bsd.zbio, bsd, swi, res, ktempM, QtempM, a2_P, b2_P, Phi2_P, alphaM);
                        
                    else        % anoxic layer -> call M first
                        
                        [e_zbio_l1_M, dedz_zbio_l1_M, f_zbio_l1_M, dfdz_zbio_l1_M, g_zbio_l1_M, dgdz_zbio_l1_M, ...
                            p_zbio_l1_M, dpdz_zbio_l1_M, q_zbio_l1_M, dqdz_zbio_l1_M, a1_M, b1_M] ...
                            = obj.calcfg_l1_M(bsd.zbio, bsd, swi, res, ls.D1M, ktempM, QtempM, 0, 0, 0, alphaM);
                        
                        [e_zbio_l2_M, dedz_zbio_l2_M, f_zbio_l2_M, dfdz_zbio_l2_M, g_zbio_l2_M, dgdz_zbio_l2_M, ...
                            p_zbio_l2_M, dpdz_zbio_l2_M, q_zbio_l2_M, dqdz_zbio_l2_M, a2_M] ...
                            = obj.calcfg_l2_M(bsd.zbio, bsd, swi, res, ktempM, QtempM, 0, 0, 0, alphaM);
                        
                        [e_zbio_l1_P, dedz_zbio_l1_P, f_zbio_l1_P, dfdz_zbio_l1_P, g_zbio_l1_P, ...
                            dgdz_zbio_l1_P, p_zbio_l1_P, dpdz_zbio_l1_P, q_zbio_l1_P, dqdz_zbio_l1_P, a1_P, b1_P, Phi1_P] ...
                            = obj.calcfg_l1_PO4(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D1P, ktempP, QtempP, a1_M, b1_M, alphaP);
                        
                        [e_zbio_l2_P, dedz_zbio_l2_P, f_zbio_l2_P, dfdz_zbio_l2_P, g_zbio_l2_P, ...
                            dgdz_zbio_l2_P, p_zbio_l2_P, dpdz_zbio_l2_P, q_zbio_l2_P, dqdz_zbio_l2_P, a2_P, b2_P, Phi2_P] ...
                            = obj.calcfg_l2_PO4(bsd.zbio, bsd, swi, res, reac1, reac2, ls.D2P, ktempP, QtempP, a2_M, alphaP);
                        
                    end
                    
                    % match solutions at zbio - continuous concentration and flux
                    % organize the data in matrices and let matlab do the calculation
                    %  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity
                    %  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux
                    %  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity
                    %  |       x16|   | D_l |      |        y16| | D_r|    |z4|    SD always M _diffusive_ flux  = 0 (cf org C)
                    
                    % discontinuity constants
                    Vb = 0;
                    Fb = 0;
                    
                    new_version = false;      % Dominik 08.02.2016: Use old version as for match at zbio always 3 int. const. at lower layer
                    if (new_version == true)
                        % Dominik 28.01.2016: now with diff case depending on oxydation as
                        % in zPO4_M
                        if(res.zox <= bsd.zbio)   % 1. CASE: 4 int const. in each layer
                            % SD zox is bioturbated
                            X = [e_zbio_l1_P, f_zbio_l1_P, p_zbio_l1_P, q_zbio_l1_P; ...
                                ls.D1P*dedz_zbio_l1_P, ls.D1P*dfdz_zbio_l1_P, ls.D1P*dpdz_zbio_l1_P, ls.D1P*dqdz_zbio_l1_P; ...
                                e_zbio_l1_M, f_zbio_l1_M, p_zbio_l1_M, q_zbio_l1_M; ...
                                ls.D1M*dedz_zbio_l1_M, ls.D1M*dfdz_zbio_l1_M, ls.D1M*dpdz_zbio_l1_M, ls.D1M*dqdz_zbio_l1_M];
                            Y = [e_zbio_l2_P, f_zbio_l2_P, p_zbio_l2_P, q_zbio_l2_P; ...
                                ls.D2P*dedz_zbio_l2_P, ls.D2P*dfdz_zbio_l2_P, ls.D2P*dpdz_zbio_l2_P, ls.D2P*dqdz_zbio_l2_P; ...
                                e_zbio_l2_M, f_zbio_l2_M, p_zbio_l2_M, q_zbio_l2_M; ...
                                ls.D2M*dedz_zbio_l2_M, ls.D2M*dfdz_zbio_l2_M, ls.D2M*dpdz_zbio_l2_M, ls.D2M*dqdz_zbio_l2_M];
                            Z = [g_zbio_l2_P-g_zbio_l1_P + Vb; ...
                                ls.D2P*dgdz_zbio_l2_P - ls.D1P*dgdz_zbio_l1_P + Fb - bsd.w.*Vb; ...
                                g_zbio_l2_M-g_zbio_l1_M + Vb; ...
                                ls.D2M*dgdz_zbio_l2_M - ls.D1M*dgdz_zbio_l1_M + Fb - bsd.w.*Vb];
                            case_flag=1;
                        else    % 2. CASE: 3 int const. in each layer
                            % SD zox non-bioturbated
                            % SD this should generate 3x3 matrices as no M flux bc (and then 4x4 C with zeros)
                            X = [e_zbio_l1_P, f_zbio_l1_P, p_zbio_l1_P; ...
                                ls.D1P*dedz_zbio_l1_P, ls.D1P*dfdz_zbio_l1_P, ls.D1P*dpdz_zbio_l1_P; ...
                                e_zbio_l1_M, f_zbio_l1_M, p_zbio_l1_M];
                            Y = [e_zbio_l2_P, f_zbio_l2_P, p_zbio_l2_P; ...
                                ls.D2P*dedz_zbio_l2_P, ls.D2P*dfdz_zbio_l2_P, ls.D2P*dpdz_zbio_l2_P; ...
                                e_zbio_l2_M, f_zbio_l2_M, p_zbio_l2_M];
                            Z = [g_zbio_l2_P-g_zbio_l1_P + Vb; ...
                                ls.D2P*dgdz_zbio_l2_P - ls.D1P*dgdz_zbio_l1_P + Fb - bsd.w.*Vb; ...
                                g_zbio_l2_M-g_zbio_l1_M + Vb];
                            
                            case_flag=2;
                        end
                        
                        [ls.C, ls.D] = benthic_utils.matchsoln_PO4_M(X, Y, Z,case_flag);
                    else
                        
                        % Dominik 28.01.2016 that's how it was  : Dominik 08.02.2016: add here M-diffusive flux always zero!!!
                        X = [e_zbio_l1_P, f_zbio_l1_P, p_zbio_l1_P, q_zbio_l1_P; ...
                            ls.D1P*dedz_zbio_l1_P, ls.D1P*dfdz_zbio_l1_P, ls.D1P*dpdz_zbio_l1_P, ls.D1P*dqdz_zbio_l1_P; ...
                            e_zbio_l1_M, f_zbio_l1_M, p_zbio_l1_M, q_zbio_l1_M; ...
                            ls.D1M*dedz_zbio_l1_M, ls.D1M*dfdz_zbio_l1_M, ls.D1M*dpdz_zbio_l1_M, ls.D1M*dqdz_zbio_l1_M];
                        Y = [e_zbio_l2_P, f_zbio_l2_P, p_zbio_l2_P, q_zbio_l2_P; ...
                            ls.D2P*dedz_zbio_l2_P, ls.D2P*dfdz_zbio_l2_P, ls.D2P*dpdz_zbio_l2_P, ls.D2P*dqdz_zbio_l2_P; ...
                            e_zbio_l2_M, f_zbio_l2_M, p_zbio_l2_M, q_zbio_l2_M; ...
                            ls.D2M*dedz_zbio_l2_M, ls.D2M*dfdz_zbio_l2_M, ls.D2M*dpdz_zbio_l2_M, ls.D2M*dqdz_zbio_l2_M];
                        Z = [g_zbio_l2_P-g_zbio_l1_P + Vb; ...
                            ls.D2P*dgdz_zbio_l2_P - ls.D1P*dgdz_zbio_l1_P + Fb - bsd.w.*Vb; ...
                            g_zbio_l2_M-g_zbio_l1_M + Vb; ...
                            ls.D2M*dgdz_zbio_l2_M - ls.D1M*dgdz_zbio_l1_M + Fb - bsd.w.*Vb];
                        
                        [ls.C, ls.D] = benthic_utils.matchsoln_PO4_M(X, Y, Z, 1);
                    end
                    
                end  % end for bioturbation cases
                
            else
                
                %%% would be the vectorized version
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c1, d1] ...
                = calcfg_l1_M(obj, z, bsd, swi, res, Dtemp, ktemp, Qtemp, a1_P, b1_P, Phi1_P, alpha)
            % Basis functions for solutes, case z <= zbio
            %
            % reac1, reac2        - mol./mol S released per organic carbon C
            % depend1,   depend2 coming from other species
            %
            % General solution for solute S is given by
            %       S(z) = A .* e(z) + B .* f(z) + g(z)
            % and for dependent species
            %       S(z) = A .* e(z) + B .* f(z) + C .* p(z) +  D.* q(z) + g(z)
            %
            %  was: S(z) = C .* h(z) + D .* f(z) + A .* h(z) +  B.* i(z) + g(z)
            
            r = res.rTOC;
            
            % SD this does work for oxic case, as ktemp == 0 (and Qtemp == 0) ?
            c1=(bsd.w-sqrt(bsd.w.^2+4.*Dtemp.*ktemp))/(2.*Dtemp);
            % p = const as c1==0 for oxic case
            p=exp(z.*c1);           % was e = ones(1,bsd.ncl);
            dpdz = c1.*exp(z.*c1);       % was zeros(1,bsd.ncl);
            
            d1=(bsd.w+sqrt(bsd.w.^2+4.*Dtemp.*ktemp))/(2.*Dtemp);  % was bsd.w./Dtemp;
            q=exp(z.*d1);
            dqdz = d1.*exp(z.*d1);
            
            
            %pfac=1./bsd.por;   % assume org matter already .*(1-bsd.por)
            %            pfac = 1;          % in fact, already has (1-por)/por
            
            if(alpha ~= 0) % oxic layer: was z<=res.zox BUT problems at boundary. M is dependent on PO4
                c1=0;
                d1=0;
                % Dominik Dec 2015: Change 2 added -1 to existing alpha
                e = -alpha/(Dtemp.*a1_P.^2-bsd.w.*a1_P-ktemp).*exp(z.*a1_P);
                dedz = a1_P.*e;
                f = -alpha/(Dtemp.*b1_P.^2-bsd.w.*b1_P-ktemp).*exp(z.*b1_P);
                dfdz = b1_P.*f;
                
                % this is when M is dependent! % Dominik 28.01.2016: Think this
                % should not be in exp()!
                %    Dominik 28.01.2016 was before: and this was in g and dgdz instead of
                %    r.a11, r.b11, ...
                %                 ea11z = exp(r.a11.*z);
                %                 eb11z = exp(r.b11.*z);
                %                 ea12z = exp(r.a12.*z);
                %                 eb12z = exp(r.b12.*z);
                
                % SD - suspect this is missing a scaling with ksPO4  ?
                % (I'd expect \propto ksPO4 ?)
                % also looks wrong in the doc ....
                % Dominik 18.12.2015: Change 1 added -alpha here //
                % Dominik 20.01.2016 no just *-1
                g = -alpha*(Phi1_P.PhiI1/(Dtemp.*r.a11.^2-bsd.w.*r.a11-ktemp).*exp(z.*r.a11) + Phi1_P.PhiII1/(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp).*exp(z.*r.b11) + ...
                    Phi1_P.PhiIII1/(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp).*exp(z.*r.b11) + ...
                    Phi1_P.PhiI2/(Dtemp.*r.a12.^2-bsd.w.*r.a12-ktemp).*exp(z.*r.a12) + Phi1_P.PhiII2/(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp).*exp(z.*r.b12) + ...
                    Phi1_P.PhiIII2/(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp).*exp(z.*r.b12));
                dgdz = -alpha*(Phi1_P.PhiI1/(Dtemp.*r.a11.^2-bsd.w.*r.a11-ktemp).*exp(z.*r.a11).*r.a11 + Phi1_P.PhiII1/(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp).*exp(z.*r.b11).*r.b11 + ...
                    Phi1_P.PhiIII1/(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp).*exp(z.*r.b11).*r.b11 + ...
                    Phi1_P.PhiI2/(Dtemp.*r.a12.^2-bsd.w.*r.a12-ktemp).*exp(z.*r.a12).*r.a12 + Phi1_P.PhiII2/(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp).*exp(z.*r.b12).*r.b12 + ...
                    Phi1_P.PhiIII2/(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp).*exp(z.*r.b12).*r.b12);
                
            else        %anoxic layer: M is independent of PO4 (no value in alpha!)
                g = Qtemp/ktemp;
                dgdz = 0;
                e = 0;
                dedz = 0;
                f = 0;
                dfdz = 0;
            end
            
        end
        
        
        function [ e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c2] ...
                = calcfg_l2_M(obj, z, bsd, swi, res, ktemp, Qtemp, a2_P, b2_P, Phi2_P, alpha)
            % Basis functions for solutes, case z > zbio
            % reac1, reac2        - mol/mol S released per organic carbon C
            %
            % General solution for solute S is given by
            %  S(z) = A .* e(z) + B .* f(z) + g(z)
            
            r = res.rTOC;
            c2=0;
            
            if(alpha ~= 0)  % was z<=res.zox M is dependent of PO4
                
                % Dominik Dec 2015: Change 2
                e=alpha./(bsd.w.*a2_P).*exp(z.*a2_P);
                dedz = a2_P.*e;
                
                f=alpha./(bsd.w.*b2_P).*exp(z.*b2_P);
                dfdz = b2_P.*f;
                
                p = 1;  % CHECK/TODO: integration constant just C
                dpdz = 0;
                q=0;
                dqdz=0;
                
                %pfac=1./bsd.por;   % assume org matter already .*(1-bsd.por)
                %                pfac = 1;          % in fact, already has (1-por)/por
                
                % Dominik 18.12.2015: Change 1: added the minus here to the already existing alpha
                % Dominik 20.01.2016: no -1
                g = alpha./(bsd.w).*(Phi2_P.PhiI1./(r.a21).*exp(r.a21.*z) + ...
                    Phi2_P.PhiI2./(r.a22).*exp(r.a22.*z));
                dgdz = alpha./(bsd.w).*(Phi2_P.PhiI1.*exp(r.a21.*z) + ...
                    Phi2_P.PhiI2.*exp(r.a22.*z));
                
            else    % z > res.zox - M is independent of PO4
                c2=-ktemp/bsd.w;
                p=exp(c2.*z);           % was e = ones(1,bsd.ncl);
                dpdz = c2.*exp(c2.*z);       % was zeros(1,bsd.ncl);
                q=0;
                dqdz=0;
                g = Qtemp/ktemp;
                dgdz = 0;
                e=0;
                dedz=0;
                f=0;
                dfdz=0;
            end
            
        end
        
        
        
        
        function [e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a1, b1, Phi1] ...
                = calcfg_l1_PO4(obj, z, bsd, swi, res, reac1, reac2, Dtemp, ktemp, Qtemp, a1_M, b1_M, alpha)
            % Basis functions for solutes, case z <= zbio
            %
            % reac1, reac2        - mol./mol S released per organic carbon C
            % depend1,   depend2 coming from other species
            %
            % General solution for solute S is given by
            %       S(z) = A .* e(z) + B .* f(z) + g(z)
            % and for dependent species
            %       S(z) = A .* e(z) + B .* f(z) + C .* p(z) +  D.* q(z) + g(z)
            
            r = res.rTOC;
            
            a1=(bsd.w-sqrt(bsd.w.^2+4.*Dtemp.*ktemp))/(2.*Dtemp);
            e=exp(z.*a1);           % was e = ones(1,bsd.ncl);
            dedz = a1.*exp(z.*a1);       % was zeros(1,bsd.ncl);
            
            b1=(bsd.w+sqrt(bsd.w.^2+4.*Dtemp.*ktemp))/(2.*Dtemp);  %was bsd.w./Dtemp;
            f=exp(z.*b1);
            dfdz = b1.*exp(z.*b1);
            
            %pfac=1./bsd.por;   % assume org matter already .*(1-bsd.por)
            pfac = 1;          % in fact, already has (1-por)/por
            
            % NOW to OM reaction terms!
            % save all Phis in one variable to pass back
            Phi1.PhiI1   =-pfac.*obj.k1.*(reac1).*r.A11./(Dtemp.*r.a11.^2-bsd.w.*r.a11-ktemp);
            Phi1.PhiII1  = pfac.*obj.k1.*(reac1).*r.A11./(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp);
            Phi1.PhiIII1 =-pfac.*obj.k1.*(reac1).*swi.C01./(Dtemp.*r.b11.^2-bsd.w.*r.b11-ktemp);
            Phi1.PhiI2   =-pfac.*obj.k2.*(reac2).*r.A12./(Dtemp.*r.a12.^2-bsd.w.*r.a12-ktemp);
            Phi1.PhiII2  = pfac.*obj.k2.*(reac2).*r.A12./(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp);
            Phi1.PhiIII2 =-pfac.*obj.k2.*(reac2).*swi.C02./(Dtemp.*r.b12.^2-bsd.w.*r.b12-ktemp);
            
            ea11z = exp(r.a11.*z);
            eb11z = exp(r.b11.*z);
            ea12z = exp(r.a12.*z);
            eb12z = exp(r.b12.*z);
            
            if(ktemp==0) % CHECK: actually no need as ktemp always <> 0
                g =  Phi1.PhiI1.*ea11z + Phi1.PhiII1.*eb11z + Phi1.PhiIII1.*eb11z + ...
                    Phi1.PhiI2.*ea12z + Phi1.PhiII2.*eb12z + Phi1.PhiIII2.*eb12z;
            else
                g =  Phi1.PhiI1.*ea11z + Phi1.PhiII1.*eb11z + Phi1.PhiIII1.*eb11z + ...
                    Phi1.PhiI2.*ea12z + Phi1.PhiII2.*eb12z + Phi1.PhiIII2.*eb12z + Qtemp./ktemp;   % here problem if ktemp=0
            end
            dgdz = Phi1.PhiI1.*r.a11.*ea11z + Phi1.PhiII1.*r.b11.*eb11z + Phi1.PhiIII1.*r.b11.*eb11z + ...
                Phi1.PhiI2.*r.a12.*ea12z + Phi1.PhiII2.*r.b12.*eb12z + Phi1.PhiIII2.*r.b12.*eb12z;
            
            
            if(alpha == 0) % was z<=res.zox PO4 is independent of M (no info in alpha)
                p = 0;
                dpdz = 0;
                q = 0;
                dqdz = 0;
                
            else        % PO4 is dependent on M
                %%%%% dependent (if existent) otherwise ls_old should be zero
                % % % %                 p = e;
                % % % %                 dpdz = dedz;
                % % % %                 q = f;
                % % % %                 dqdz = dfdz;
                % Dominik 21.12.2015: Change 3: added *-1 here for p qnd q (as sign
                % changes for exponential part)
                p = -alpha/(Dtemp.*a1_M.^2-bsd.w.*a1_M-ktemp).*exp(z.*a1_M);
                dpdz = a1_M.*p;
                q = -alpha/(Dtemp.*b1_M.^2-bsd.w.*b1_M-ktemp).*exp(z.*b1_M);
                dqdz = b1_M.*q;
                
            end
            
        end
        
        
        function [ e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a2, b2, Phi2] ...
                = calcfg_l2_PO4(obj, z, bsd, swi, res, reac1, reac2, Dtemp, ktemp, Qtemp, a2_M, alpha)
            % Basis functions for solutes, case z > zbio
            % reac1, reac2        - mol/mol S released per organic carbon C
            %
            % General solution for solute S is given by
            %  S(z) = A .* e(z) + B .* f(z) + g(z)
            
            r = res.rTOC;
            
            a2=(bsd.w-sqrt(bsd.w.^2+4.*Dtemp.*ktemp))/(2.*Dtemp);
            e=exp(z.*a2);           % was e = ones(1,bsd.ncl);
            dedz = a2.*exp(z.*a2);       % was zeros(1,bsd.ncl);
            
            b2=(bsd.w+sqrt(bsd.w.^2+4.*Dtemp.*ktemp))/(2.*Dtemp);  %was bsd.w./Dtemp;
            f=exp(z.*b2);
            dfdz = b2.*exp(z.*b2);
            
            
            %pfac=1./bsd.por;   % assume org matter already .*(1-bsd.por)
            pfac = 1;          % in fact, already has (1-por)/por
            
            Phi2.PhiI1 = -pfac.*obj.k1.*(reac1).*r.A21./(Dtemp.*r.a21.^2-bsd.w.*r.a21-ktemp);
            Phi2.PhiI2 = -pfac.*obj.k2.*(reac2).*r.A22./(Dtemp.*r.a22.^2-bsd.w.*r.a22-ktemp);
            
            if(ktemp==0)    % CHECK: think no need for this as always ktemp <> 0
                g = Phi2.PhiI1.*exp(r.a21.*z) + ...
                    Phi2.PhiI2.*exp(r.a22.*z);
            else
                g = Phi2.PhiI1.*exp(r.a21.*z) + ...
                    Phi2.PhiI2.*exp(r.a22.*z) + Qtemp./ktemp;
            end
            dgdz = Phi2.PhiI1.*r.a21.*exp(r.a21.*z) + ...
                Phi2.PhiI2.*r.a22.*exp(r.a22.*z);
            
            if(alpha == 0) % was z<=res.zox PO4 is independent of M (no info in alpha)
                p = 0;
                dpdz = 0;
                q = 0;
                dqdz = 0;
            else            % PO4 is dependent on M
                % % % %                 p = e;
                % % % %                 dpdz = dedz;
                % % % %                 q = f;
                % % % %                 dqdz = dfdz;
                % Dominik 21.12.2015: Change 3: added *-1 here (as sign
                % changes for exponential part)
                p = -alpha/(Dtemp.*a2_M.^2-bsd.w.*a2_M-ktemp).*exp(a2_M.*z);
                dpdz = a2_M.*p;
                % subcase zox&zbio < z: here just one term in M
                q=0;
                dqdz=0;
                
                %f = alpha/(Dtemp.*ls_old.b1.^2-bsd.w.*ls_old.b1-ktemp).*exp(z.*ls_old.b1);
                %dfdz = ls_old.b1.*f;
                
            end
            
        end
        
    end
end

classdef benthic_zPO4_M
    % Solve PO4
    % check if PO4 < 0 not included here
    
    properties
        qdispPO4=112.90764;         % PO4 diffusion coefficient in water (cm2/yr)
        adispPO4=5.586252           % PO4 linear coefficient for temperature dependence (cm2/yr/oC)
        DPO41;                      % PO4 diffusion coefficient in bioturbated layer (cm2/yr)
        DPO42;                      % PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        KPO4_ox=200.0;              % Adsorption coefficient in oxic layer (-) (was 10-0)
        KPO4_anox=1.3;              % Adsorption coefficient in anoxic layer (-)
        ksPO4=3.65; %7.3; % 3.65 0.26*365;           	% Rate constant for kinetic P sorption (1/yr) (Palastanga: 3.65;)
        kmPO4=0.193;               	% Rate constant for Fe-bound P release upon Fe oxide reduction
        kaPO4=0.37; %0.365;             	%Rate constant for authigenic P formation (1/yr)
%        ksPO4=1e-15;
%        kmPO4= 1e-15 ;
%        kaPO4 = 0.0;
        PO4s=2.0e-9; %12.0e-9; % 10.0e-10;              % Equilibrium concentration for P sorption (mol/cm3)       was 1.5e-9; ; Slomp ea 1996
        PO4a= 3.7e-9; % 3.7e-9;                % Equilibrium concentration for authigenic P formation (mol/cm3) was 0.7e-9
        Minf=5.2e-9; % 0.0; 5.2e-9; %1.99e-10; %2.0e-9; %1.0e-10;       % asymptotic concentration for Fe-bound P (mol/cm3)  (was 5.2e-9;)
        %Minf = 0;
        
        % OM reactive terms in oxic layer
        reac1_ox;
%         reac2_ox;
        % OM reactive terms in anoxic layer
        reac1_anox;
%         reac2_anox;
    end
    
    methods
        function obj = benthic_zPO4_M(bsd, swi)
            obj.DPO41=((obj.qdispPO4+obj.adispPO4*swi.T).*bsd.dispFactor+bsd.Dbio); 	% PO4 diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DPO42=((obj.qdispPO4+obj.adispPO4*swi.T).*bsd.dispFactor);           	% PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
            
            %reactive terms: OM degradation
            obj.reac1_ox=1/(1+obj.KPO4_ox)*bsd.PC1;
%             obj.reac2_ox=1/(1+obj.KPO4_ox)*bsd.PC2;
            obj.reac1_anox=1/(1+obj.KPO4_anox)*bsd.PC1;
%             obj.reac2_anox=1/(1+obj.KPO4_anox)*bsd.PC2;
            
        end
        
        function r = calc(obj, bsd, swi, r)
            
%         if(r.zox == bsd.zinf)
% %            obj.Minf=1.0e-10; %1.99e-6;       % asymptotic concentration for Fe-bound P (mol/cm3)      TODO/CHECK: good value? is from Slomp et al. 1996 Dom was 1.99e-6
%             rPO4_M.ls1 = r.zTOC.prepfg_l12_PO4_M(bsd, swi, r, obj.reac1_ox, obj.reac2_ox, obj.ksPO4/(1+obj.KPO4_ox), obj.PO4s*obj.ksPO4/(1+obj.KPO4_ox),0, r.zox, obj.DPO41/(1+obj.KPO4_ox), obj.DPO42/(1+obj.KPO4_ox), 0, ...
%                 0, 0, bsd.Dbio, 0, (1/bsd.SD)*obj.ksPO4);
% 
%             % Solution at swi, top of layer 1
%             [ e1_0_P, dedz1_0_P, f1_0_P, dfdz1_0_P, g1_0_P, dgdz1_0_P, p1_0_P, dpdz1_0_P, q1_0_P, dqdz1_0_P, ...
%                 e1_0_M, dedz1_0_M, f1_0_M, dfdz1_0_M, g1_0_M, dgdz1_0_M, p1_0_M, dpdz1_0_M, q1_0_M, dqdz1_0_M]...
%                 = r.zTOC.calcfg_l12_PO4_M(0, bsd, swi, r, obj.reac1_ox, obj.reac2_ox, obj.ksPO4/(1+obj.KPO4_ox), obj.PO4s*obj.ksPO4/(1+obj.KPO4_ox), 0,  rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
% 
%             % Basis functions at bottom of layer 2 zinf
%             [ e2_zinf_P, dedz2_zinf_P, f2_zinf_P, dfdz2_zinf_P, g2_zinf_P, dgdz2_zinf_P, p2_zinf_P, dpdz2_zinf_P, q2_zinf_P, dqdz2_zinf_P, ...
%                 e2_zinf_M, dedz2_zinf_M, f2_zinf_M, dfdz2_zinf_M, g2_zinf_M, dgdz2_zinf_M, p2_zinf_M, dpdz2_zinf_M, q2_zinf_M, dqdz2_zinf_M] ...
%                 = r.zTOC.calcfg_l12_PO4_M(bsd.zinf, bsd, swi, r, obj.reac1_anox, obj.reac2_anox, obj.kaPO4/(1+obj.KPO4_anox), obj.PO4a*obj.kaPO4/(1+obj.KPO4_anox), bsd.SD*obj.kmPO4/(1+obj.KPO4_anox), rPO4_M.ls2, obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
% 
%         else
%            obj.Minf=1.0e-10;       % asymptotic concentration in anoxic conditions
%        end
            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            
            % layer 1: 0 < z < zox, OM degradation (-) Sorption to sediment Fe-oxides (ktemp)
            %       ls =         prepfg_l12_PO4M(bsd, swi, r,    reac1P,     ktempP                 ,     QP,                             zU, zL,            D1P,                 D2P ,                alphaP
            rPO4_M.ls1 = r.zTOC.prepfg_l12_PO4_M(bsd, swi, r, obj.reac1_ox, obj.ksPO4/(1+obj.KPO4_ox), obj.PO4s*obj.ksPO4/(1+obj.KPO4_ox),0, r.zox, obj.DPO41/(1+obj.KPO4_ox), obj.DPO42/(1+obj.KPO4_ox), 0, ...
                0, 0, bsd.Dbio, 0, (1/bsd.SD)*obj.ksPO4);
            % for M  ktempM, QM, D1M,   D2M,  alphaM
            
            % layer 2: zox < z < zinf,
            % OM degradation (-) authigenic P formation (ktemp) (+) P desorption due to Fe-bound P release upon Fe oxide reduction
            
            %             ls =  prepfg_l12_PO4M( bsd, swi, r,     reac1P,      ktempP                ,     QP,                                   zU,   zL,        D1P,       D2P,    alphaP
            rPO4_M.ls2 = r.zTOC.prepfg_l12_PO4_M(bsd, swi, r, obj.reac1_anox, obj.kaPO4/(1+obj.KPO4_anox), obj.PO4a*obj.kaPO4/(1+obj.KPO4_anox), r.zox, bsd.zinf, obj.DPO41/(1+obj.KPO4_anox), obj.DPO42/(1+obj.KPO4_anox), bsd.SD*obj.kmPO4/(1+obj.KPO4_anox), ...
                obj.kmPO4, obj.kmPO4.*obj.Minf, bsd.Dbio, 0, 0);
            % for M           ktempM,       QM,              D1M,  D2M, alphaM)
            
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 2 zinf
            [ e2_zinf_P, dedz2_zinf_P, f2_zinf_P, dfdz2_zinf_P, g2_zinf_P, dgdz2_zinf_P, p2_zinf_P, dpdz2_zinf_P, q2_zinf_P, dqdz2_zinf_P, ...
                e2_zinf_M, dedz2_zinf_M, f2_zinf_M, dfdz2_zinf_M, g2_zinf_M, dgdz2_zinf_M, p2_zinf_M, dpdz2_zinf_M, q2_zinf_M, dqdz2_zinf_M] ...
                = r.zTOC.calcfg_l12_PO4_M(bsd.zinf, bsd, swi, r, obj.reac1_anox, obj.kaPO4/(1+obj.KPO4_anox), obj.PO4a*obj.kaPO4/(1+obj.KPO4_anox), bsd.SD*obj.kmPO4/(1+obj.KPO4_anox), rPO4_M.ls2, obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
            % calcfg_l12_PO4_M(obj,            z,   bsd, swi, res,     reac1P,   ktempP,                        QtempP,                     alphaP,                                  ls,      ktempM,           QtempM,    alphaM)
            
            % Match at zox, layer 1 - layer 2 (continuity and flux)
            % basis functions at bottom of layer 1
            [ e1_zox_P, dedz1_zox_P, f1_zox_P, dfdz1_zox_P, g1_zox_P, dgdz1_zox_P, p1_zox_P, dpdz1_zox_P, q1_zox_P, dqdz1_zox_P, ...
                e1_zox_M, dedz1_zox_M, f1_zox_M, dfdz1_zox_M, g1_zox_M, dgdz1_zox_M, p1_zox_M, dpdz1_zox_M, q1_zox_M, dqdz1_zox_M]...
                = r.zTOC.calcfg_l12_PO4_M(r.zox, bsd, swi, r,  obj.reac1_ox, obj.ksPO4/(1+obj.KPO4_ox), obj.PO4s*obj.ksPO4/(1+obj.KPO4_ox), 0, rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
            % calcfg_l12_PO4_M(obj, z, bsd, swi, res,          reac1P,       ktempP,                     QtempP,                alphaP,       ls,  ktempM, QtempM, alphaM)
            
            % ... and top of layer 2
            [ e2_zox_P, dedz2_zox_P, f2_zox_P, dfdz2_zox_P, g2_zox_P, dgdz2_zox_P, p2_zox_P, dpdz2_zox_P, q2_zox_P, dqdz2_zox_P, ...
                e2_zox_M, dedz2_zox_M, f2_zox_M, dfdz2_zox_M, g2_zox_M, dgdz2_zox_M, p2_zox_M, dpdz2_zox_M, q2_zox_M, dqdz2_zox_M] ...
                = r.zTOC.calcfg_l12_PO4_M(r.zox, bsd, swi, r, obj.reac1_anox, obj.kaPO4/(1+obj.KPO4_anox), obj.PO4a*obj.kaPO4/(1+obj.KPO4_anox), bsd.SD*obj.kmPO4/(1+obj.KPO4_anox), rPO4_M.ls2, obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
            % calcfg_l12_PO4_M(obj, z, bsd, swi, res,         reac1P,         ktempP,                        QtempP,                          alphaP,                           ls,         ktempM,           QtempM,    alphaM)
            
            % match solutions at zox - continuous concentration and flux
            % organize the data in matrices and let matlab do the calculation
            %  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity
            %  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux
            %  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity
            %  |       x16|   | D_l |      |        y16| | D_r|    |z4|    SD M flux  only in bioturbated case, otherwise not an  independent constraint
            
            % discontinuity constants
            Vb = 0;
            Fb = 0;
            
            if(r.zox < bsd.zbio)   % 1. CASE: 4 int const. in each layer
                % SD zox is bioturbated
                X = [e1_zox_P, f1_zox_P, p1_zox_P, q1_zox_P; ...
                    dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P, dqdz1_zox_P; ...
                    e1_zox_M, f1_zox_M, p1_zox_M, q1_zox_M; ...
                    dedz1_zox_M, dfdz1_zox_M, dpdz1_zox_M, dqdz1_zox_M];
                Y = [e2_zox_P, f2_zox_P, p2_zox_P, q2_zox_P; ...
                    dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P, dqdz2_zox_P; ...
                    e2_zox_M, f2_zox_M, p2_zox_M, q2_zox_M; ...
                    dedz2_zox_M, dfdz2_zox_M, dpdz2_zox_M, dqdz2_zox_M];
                Z = [g2_zox_P-g1_zox_P + Vb; ...
                    dgdz2_zox_P - dgdz1_zox_P + Fb - bsd.w.*Vb; ...
                    g2_zox_M-g1_zox_M + Vb; ...
                    dgdz2_zox_M - dgdz1_zox_M + Fb - bsd.w.*Vb];
                case_flag=1;
            else    % 2. CASE: 3 int const. in each layer
                % SD zox non-bioturbated
                % SD this should generate 3x3 matrices as no M flux bc (and then 4x4 C with zeros)
                X = [e1_zox_P, f1_zox_P, p1_zox_P; ...
                    dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P; ...
                    e1_zox_M, f1_zox_M, p1_zox_M];
                Y = [e2_zox_P, f2_zox_P, p2_zox_P; ...
                    dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P; ...
                    e2_zox_M, f2_zox_M, p2_zox_M];
                Z = [g2_zox_P-g1_zox_P + Vb; ...
                    dgdz2_zox_P - dgdz1_zox_P + Fb - bsd.w.*Vb; ...
                    g2_zox_M-g1_zox_M + Vb];
                % SD old code - this _might_ work, but only because the two
                % M bc are degenerate anyway...
                %                 X = [e1_zox_P, f1_zox_P, p1_zox_P; ...
                %                      dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P; ...
                %                      e1_zox_M, f1_zox_M, p1_zox_M; ...
                %                      dedz1_zox_M, dfdz1_zox_M, dpdz1_zox_M];
                %                 Y = [e2_zox_P, f2_zox_P, p2_zox_P; ...
                %                      dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P; ...
                %                      e2_zox_M, f2_zox_M, p2_zox_M; ...
                %                      dedz2_zox_M, dfdz2_zox_M, dpdz2_zox_M];
                %                 Z = [g2_zox_P-g1_zox_P + Vb; ...
                %                      dgdz2_zox_P - dgdz1_zox_P + Fb - bsd.w.*Vb; ...
                %                      g2_zox_M-g1_zox_M + Vb; ...
                %                      dgdz2_zox_M - dgdz1_zox_M + Fb - bsd.w.*Vb];
                case_flag=2;
            end
            
            [zox.C, zox.D] = benthic_utils.matchsoln_PO4_M(X, Y, Z,case_flag);
            
            
            % Solution at swi, top of layer 1
            [ e1_0_P, dedz1_0_P, f1_0_P, dfdz1_0_P, g1_0_P, dgdz1_0_P, p1_0_P, dpdz1_0_P, q1_0_P, dqdz1_0_P, ...
                e1_0_M, dedz1_0_M, f1_0_M, dfdz1_0_M, g1_0_M, dgdz1_0_M, p1_0_M, dpdz1_0_M, q1_0_M, dqdz1_0_M]...
                = r.zTOC.calcfg_l12_PO4_M(0, bsd, swi, r, obj.reac1_ox, obj.ksPO4/(1+obj.KPO4_ox), obj.PO4s*obj.ksPO4/(1+obj.KPO4_ox), 0,  rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
            % calcfg_l12_PO4_M(obj, z, bsd, swi, res,      reac1P,      ktempP,                       QtempP,             alphaP,   ls,       tempM, QtempM, alphaM)
            
            % transform to use coeffs from l2
            % Now find 'transformed' basis functions such that in layer 1, O2 = A_2*et + B_2*ft + gt  (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
            
            EFPQ_P = [e1_0_P; f1_0_P; p1_0_P; q1_0_P];
            dEFPQdz_P = [dedz1_0_P; dfdz1_0_P; dpdz1_0_P; dqdz1_0_P];
            EFPQ_M = [e1_0_M; f1_0_M; p1_0_M; q1_0_M];
            dEFPQdz_M = [dedz1_0_M; dfdz1_0_M; dpdz1_0_M; dqdz1_0_M];
            [EFPQ_P, g1_0_P, dEFPQdz_P, dgdz1_0_P, EFPQ_M, g1_0_M, dEFPQdz_M, dgdz1_0_M] ...
                = benthic_utils.xformsoln_PO4_M(EFPQ_P, EFPQ_M, dEFPQdz_P, dEFPQdz_M, g1_0_P, g1_0_M,dgdz1_0_P, dgdz1_0_M, zox.C, zox.D);
            
            % % % %     VERSION WITH SINGLE VALUES (NO MATRICES)
            % % % %             [e1_0_P, f1_0_P, p1_0_P, q1_0_P, g1_0_P, dedz1_0_P, dfdz1_0_P, dpdz1_0_P, dqdz1_0_P, dgdz1_0_P, e1_0_M, f1_0_M, p1_0_M, q1_0_M, g1_0_M, dedz1_0_M, dfdz1_0_M, dpdz1_0_M, dqdz1_0_M, dgdz1_0_M] ...
            % % % %             = benthic_utils.xformsoln_PO4_M(e1_0_P, f1_0_P, p1_0_P, q1_0_P, e1_0_M, f1_0_M, p1_0_M, q1_0_M, dedz1_0_P, dfdz1_0_P, dpdz1_0_P, dqdz1_0_P, dedz1_0_M, dfdz1_0_M, dpdz1_0_M, dqdz1_0_M, g1_0_P, g1_0_M,dgdz1_0_P, dgdz1_0_M, zox.C, zox.D);
            
            % SD We solve for 3 unknowns (given DPO4M == 0) with 3 bc
            % Solve for APO4M, BPO4M, CPO4M given boundary conditions (expressed in terms of transformed basis fns, layer 2 A, B, C, D)
            % BC(i) for PO4: zero flux at z=oo
            % BC(ii) PO4 concentration equal at z=0
            % BC(iii) for M: Known flux at z=0  NB: advective flux is what is needed here (bioturbation flux == 0)
            % and NOT a BC for concentration = Minf at z=oo
            
            % SD TODO include porosity factors for advective fluxes eg bsd.w*(1-bsd.por) ?
            % APO4M*dedz2_zinf_P   +  BPO4M*dfdz2_zinf_P  + CPO4M*dpdz2_zinf_P   + dgdz2_zinf_P = 0;
            % APO4M*e1_0_P         +  BPO4M*f1_0_P        + CPO4M*p1_0_P         + g1_0_P       = swi.PO40;
            % bsd.w*(APO4M*e1_0_M      +  BPO4M*f1_0_M     + CPO4M*p1_0_M      + g1_0_M )       = swi.Mflux0;
            
            % | dedz2_zinf_P dfdz2_zinf_P dpdz2_zinf_P |  |APO4M|     | - dgdz2_zinf_P       |
            % |    e1_0_P      f1_0_P        p1_0_P    |  |BPO4M|     | swi.PO40 - g1_0_P    |
            % | bsd.w*e1_0_M   bsd.w*f1_0_M  bsd.w*p1_0_M||CPO4M|     |swi.Mflux0 - bsd.w*g1_0_M|
            % and set DPO4M = 0
            
            % SD assume D2 == 0 (as q, dqdz2_zinf = 0 ) and solve for 3 unknowns
            % Dominik 01.02.2016 Check: Do I need check for cases...?
            X = [dedz2_zinf_P, dfdz2_zinf_P, dpdz2_zinf_P; ...
                EFPQ_P(1), EFPQ_P(2), EFPQ_P(3); ...
                bsd.w*EFPQ_M(1), bsd.w*EFPQ_M(2), bsd.w*EFPQ_M(3)];
            Y = [-dgdz2_zinf_P; ...
                swi.PO40 - g1_0_P; ...
                swi.Mflux0 - bsd.w*g1_0_M];
            
            
            [ rPO4_M.A2, rPO4_M.B2, rPO4_M.C2]  = benthic_utils.solve2eqn_PO4_M(X,Y);
            rPO4_M.D2 = 0;
            
                        
            % calculate concentration at zinf
            r.conczinfPO4 = rPO4_M.A2.*e2_zinf_P+rPO4_M.B2.*f2_zinf_P + g2_zinf_P;
            
            
            % DH 2405: need to check if anoxic bc of adsorption coeff? flux PO4 at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            r.flxswi_P = bsd.por.*(obj.DPO41/(1+obj.KPO4_ox).*(rPO4_M.A2.*dEFPQdz_P(1)+rPO4_M.B2.*dEFPQdz_P(2) + rPO4_M.C2.*dEFPQdz_P(3)+rPO4_M.D2.*dEFPQdz_P(4) + dgdz1_0_P) - bsd.w.*(swi.PO40 - r.conczinfPO4));   % NB: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
            
            % flux M at swi - DO include por so this is per cm^2 water column area
            r.flxswi_M = bsd.por.*bsd.Dbio*(rPO4_M.A2.*dEFPQdz_M(1)+rPO4_M.B2.*dEFPQdz_M(2) + rPO4_M.C2.*dEFPQdz_M(3)+rPO4_M.D2.*dEFPQdz_M(4) + dgdz1_0_M);   % NB: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
            %DH2405 was r.flxswi_M = bsd.por.*obj.DPO41.*(rPO4_M.A2.*dEFPQdz_M(1)+rPO4_M.B2.*dEFPQdz_M(2) + rPO4_M.C2.*dEFPQdz_M(3)+rPO4_M.D2.*dEFPQdz_M(4) + dgdz1_0_M);   % NB: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
            
            
            % save coeffs for layer 1
            L1_P = zox.C*[ rPO4_M.A2; rPO4_M.B2; rPO4_M.C2; rPO4_M.D2]+zox.D;
            
            
            rPO4_M.A1 = L1_P(1);
            rPO4_M.B1 = L1_P(2);
            rPO4_M.C1 = L1_P(3);
            rPO4_M.D1 = L1_P(4);
            
            
            r.rPO4_M = rPO4_M;
%        end
        end
        
        
        
        function FPO4 = calcFPO4(obj, zPO4, bsd, swi, r)
            % Calculate PO4 consumption below zPO4, by organic matter and indirectly via methane oxidation
            % Dom 20.01.2016 not used in code....
            tmpreac1    = bsd.PO4C.*bsd.gammaCH4;
% %             tmpreac2    = bsd.PO4C.*bsd.gammaCH4
            
            FPO4 = r.zTOC.calcReac(zPO4, bsd.zinf, tmpreac1, bsd, swi, r);
            % TODO confirm (1-bsd.por)*  has been added (to k1 & k2 ?)
        end
        
        
        
        function [PO4, flxPO4, M, flxM, e_M, f_M, p_M, q_M, g_M, dedz_M, dfdz_M, dpdz_M, dqdz_M, dgdz_M] = calcPO4_M(obj, z, bsd, swi, r)
            % Calculate PO4 concentration and flux at depth z from solution
            
            if(r.zox == bsd.zinf)
                obj.Minf=1.99e-10; %1.99e-6;       % asymptotic concentration for Fe-bound P (mol/cm3)      TODO/CHECK: good value? is from Slomp et al. 1996 Dom was 1.99e-6
            else
                obj.Minf=1.99e-10;       % asymptotic concentration in anoxic conditions
            end
            
            
            rPO4_M = r.rPO4_M;
            if z <= bsd.zbio
                D_P = obj.DPO41;
                D_M = bsd.Dbio;
            else
                D_P = obj.DPO42;
                D_M = 0;
            end
            
            
            
            if (z <= r.zox)   % layer 1
                [ e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, ...
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M]...
                    = r.zTOC.calcfg_l12_PO4_M(z, bsd, swi, r,  obj.reac1_ox, obj.ksPO4/(1+obj.KPO4_ox), obj.PO4s*obj.ksPO4/(1+obj.KPO4_ox), 0, rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
                % calcfg_l12_PO4_M(obj, z, bsd, swi, res,       reac1P,      ktempP,                     QtempP,                alphaP, ls,   ktempM, QtempM, alphaM)
                
                
                PO4     = r.rPO4_M.A1.*e_P + r.rPO4_M.B1.*f_P + r.rPO4_M.C1.*p_P + r.rPO4_M.D1.*q_P + g_P;
                flxPO4  = D_P/(1+obj.KPO4_ox).*(r.rPO4_M.A1.*dedz_P+r.rPO4_M.B1.*dfdz_P + r.rPO4_M.C1.*dpdz_P+r.rPO4_M.D1.*dqdz_P + dgdz_P);
                M     = r.rPO4_M.A1.*e_M + r.rPO4_M.B1.*f_M + r.rPO4_M.C1.*p_M + r.rPO4_M.D1.*q_M + g_M;
                flxM  = D_M.*(r.rPO4_M.A1.*dedz_M+r.rPO4_M.B1.*dfdz_M + r.rPO4_M.C1.*dpdz_M+r.rPO4_M.D1.*dqdz_M + dgdz_M);
                
            else        % layer 2
                [ e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, ...
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M]...
                    = r.zTOC.calcfg_l12_PO4_M(z, bsd, swi, r, obj.reac1_anox, obj.kaPO4/(1+obj.KPO4_anox), obj.PO4a*obj.kaPO4/(1+obj.KPO4_anox), bsd.SD*obj.kmPO4/(1+obj.KPO4_anox), rPO4_M.ls2, obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
                % calcfg_l12_PO4_M(obj, z, bsd, swi, res,     reac1P,         reac2P,          ktempP,                        QtempP,                     alphaP,                    ls,      ktempM,           QtempM,    alphaM)
                
                
                PO4     = r.rPO4_M.A2.*e_P + r.rPO4_M.B2.*f_P + r.rPO4_M.C2.*p_P + r.rPO4_M.D2.*q_P + g_P;
                flxPO4  = D_P/(1+obj.KPO4_anox).*(r.rPO4_M.A2.*dedz_P+r.rPO4_M.B2.*dfdz_P + r.rPO4_M.C2.*dpdz_P+r.rPO4_M.D2.*dqdz_P + dgdz_P);
                M     = r.rPO4_M.A2.*e_M + r.rPO4_M.B2.*f_M + r.rPO4_M.C2.*p_M + r.rPO4_M.D2.*q_M + g_M;
                flxM  = D_M.*(r.rPO4_M.A2.*dedz_M+r.rPO4_M.B2.*dfdz_M + r.rPO4_M.C2.*dpdz_M+r.rPO4_M.D2.*dqdz_M + dgdz_M);
            end
            
        end
        
    end
    
end


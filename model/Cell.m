classdef Cell
    % These are the main properties of a cell object
    properties
        c_no
        ngh
        par
        var
        fx
    end
    methods
        % Constructor Method -> This gives me an object
        function c = Cell(cell_number,neighbours,parameters)
            c.c_no = cell_number;
            c.ngh  = neighbours;
            c.par  = Cell.pardecider(cell_number,neighbours,parameters);
        end
    end
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Fluxes of the model
        function [var,f] = fluxes(c_no,neigh,var_vec,par,vtx,vty)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            var.w  = var_vec(1);
            var.xi = var_vec(2);
            var.yi = var_vec(3);
            var.zi = var_vec(4);
            var.bi = var_vec(5);
            var.ca = var_vec(6);
            var.ip = var_vec(7);
            var.hh = var_vec(8);
            var.gg = var_vec(9);
            if length(var_vec)==15
                var.xl=[var_vec(10);var_vec(13)];
                var.yl=[var_vec(11);var_vec(14)];
                var.zl=[var_vec(12);var_vec(15)];
            elseif length(var_vec)==18
                var.xl=[var_vec(10);var_vec(13);var_vec(16)];
                var.yl=[var_vec(11);var_vec(14);var_vec(17)];
                var.zl=[var_vec(12);var_vec(15);var_vec(18)];
            elseif length(var_vec)==21
                var.xl=[var_vec(10);var_vec(13);var_vec(16);var_vec(19)];
                var.yl=[var_vec(11);var_vec(14);var_vec(17);var_vec(20)];
                var.zl=[var_vec(12);var_vec(15);var_vec(18);var_vec(21)];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ce = par.gamma * ( ( par.v0 / var.w )* par.ct - var.ca );
            
            phi_c = var.ca^4 / ( var.ca^4 + par.K_c^4 );
            
            phi_p = var.ip^2 / ( par.K_p^2 + var.ip^2);
            
            phi_p_down = par.K_p^2 / ( par.K_p^2 + var.ip^2);
            
            H_inf = par.K_h^4 / ( par.K_h^4 + var.ca^4 );
            
            f.TAU = par.tau_max * ( par.K_tau^4 ) ...
                                              / ( par.K_tau^4 + var.ca^4 );
            f.J_SERCA = par.V_p * ( var.ca^2 - par.K_bar * ce^2 ) ...
                                                / ( var.ca^2 + par.K_p^2 );
            beta = phi_p * phi_c * var.hh;
            alpha = phi_p_down * ( 1 - phi_c * H_inf );
            f.Po = beta / ( beta + par.k_beta * ( beta + alpha ) );
            f.J_RyR = par.V_RyR * ( var.ca^4 / ( var.ca^4 + ...
                                                   par.K_RyR^4 ) )* var.gg;
            f.g_inf = par.K_RyR2^2 / ( par.K_RyR2^2 + var.ca^2 );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Ie = (par.Nae+par.Ke+par.Cle+par.HCO3e);
            Ii = 2*(var.xi+var.yi+par.Hy)+par.CO20;
            % Basolateral Osmolarities
            f.qb = par.Lb*(Ii-Ie);
            % Basolateral Na/K ATPase pump (NaK-ATPase)
            f.jNaK=par.Sb*par.aNaK*(par.r*par.Ke^2*var.xi^3 ...
                /(par.Ke^2+par.alpha1*var.xi^3));
            % Basolateral cotransporter (Nkcc1)
            f.jNk1=par.aNkcc1 * par.Sb ...
                * ( par.a1 - par.a2 * var.xi * var.yi * var.zi^2 ) ...
                / ( par.a3 + par.a4 * var.xi * var.yi * var.zi^2 );
            % Anion Exchanger 4 (Ae4)
            f.j4 = par.Sb * par.G4 ...
                * ( (par.Cle / (par.Cle + par.KCl)) ...
                * (var.xi / (var.xi + par.KNa)) ...
                * (var.bi / (var.bi + par.KB))^2);
            % Basolateral Na/H+ Antiporter (NHye1)
            f.j1 = par.Sb * par.G1 ...
                * ( ( par.Nae / ( par.Nae + par.KNa ) ) ...
                * ( par.Hy / ( par.KH + par.Hy ) ) ...
                - ( var.xi / ( var.xi + par.KNa ) ) ...
                *(par.Hye / ( par.KH + par.Hye)));
            % Intracellular Carbonic Amylases
            f.jB = var.w * par.GB * (par.kp * par.CO20 ...
                - par.kn * var.bi * par.Hy);
            % K+ Nernst Potential
            f.vy = par.RTF * log( par.Ke / var.yi );
            % Basolateral K+ Channels Open Probability (Mesh-Dependent)
            f.Py=sum((1./(1+(par.KCaKC./var.ca).^par.eta2))...
                .*par.Sb_k)/par.Sb;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ApicalWater Fluxes, Tight Junctional and Cl channel NPs'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for m = 1:length(neigh)
                if c_no >= neigh(m)
                    % Apical PM Osmolarity
                    f.qa(m,1) = par.La*(2*(var.xl(m)+var.yl(m)...
                        -var.xi-var.yi-par.Hy)-par.CO20+par.Ul);
                    % Tight Junctional Osmolarity
                    f.qt(m,1) = par.Lt*(2*(var.xl(m)+var.yl(m))+par.Ul...
                        - (par.Nae+par.Ke+par.Cle+par.HCO3e));
                    % Tight Junctional Na+ Nernst Potential
                    f.vtx(m,1) = par.RTF * log(var.xl(m)/par.Nae);
                    % Tight Junctional K+ Nernst Potential
                    f.vty(m,1)  = par.RTF * log(var.yl(m)/par.Ke);
                    
                    % Apical Cl- Nernst Potential
                    f.vz(m,1) = par.RTF * log(var.zl(m)/var.zi);
                    % Apical Cl Channels Open Probability (Mesh-Dependent)
                    f.Pz(m,1) = sum((1./(1+(par.KCaCC./var.ca).^par.eta1))...
                        .*par.Sa_ck{c_no,neigh(m)})./par.Sa;
                    % Total FFR
                    f.qT(m,1) = f.qa(m,1) + f.qt(m,1);
                else
                    % Apical PM Osmolarity
                    f.qa(m,1) = par.La*(2*(var.xl(m)+var.yl(m)...
                        -var.xi-var.yi-par.Hy)-par.CO20+par.Ul);
                    % Apical Cl- Nernst Potential
                    f.vz(m,1) = par.RTF * log(var.zl(m)/var.zi);
                    % Apical Cl Channels Open Probability (Mesh-Dependent)
                    f.Pz(m,1) = sum((1./(1+(par.KCaCC./var.ca).^par.eta1))...
                        .*par.Sa_ck{c_no,neigh(m)})./par.Sa;
                    % Tight Junctional Na+ Nernst Potential
                    f.vtx(m,1) = vtx(m);
                    % Tight Junctional K+ Nernst Potential
                    f.vty(m,1) = vty(m);
                    % Total FFR
                    f.qT(m,1) = f.qa(m,1);
                end
            end
            f.Qa = sum(f.qa);
            %%%%%% Parameters to generate Membrane Potentials %%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            iNaK= par.F * f.jNaK;
            n = length(neigh);
            GtX = par.GtNa* par.St;
            GtY = par.GtK * par.St;
            GY  = par.GK  * f.Py;
            ND  = n*par.St*(par.GtNa+par.GtK);
            vtx = sum(f.vtx);
            vty = sum(f.vty);
            Gz  = par.GCl * sum(f.Pz);
            Gzv = par.GCl * dot(f.Pz,f.vz);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Membrane Potentials (Explicitly)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f.va=-(ND*iNaK-GY*(GtX*vtx+GtY*vty+ND*f.vy)...
                +Gzv*(GY+ND))/(ND*(GY+Gz)+GY*Gz);
            
            f.vb=-((ND+Gz)*(iNaK-GY*f.vy)+Gz*(GtX*vtx+GtY*vty)...
                +ND*Gzv)/(ND*GY+Gz*(GY+ND));
            % PM Potential-Dependent Fluxes
            for jj = 1:length(f.vz)
                vt = f.va - f.vb;
                % Apical Cl Channels Total Flux
                f.jz(jj,1) = par.GCl*f.Pz(jj)*(f.va+f.vz(jj))/par.F;
                % Tight Junctional Na+ Flux
                f.jtx(jj,1) = par.GtNa*par.St*(vt-f.vtx(jj))/par.F;
                % Tight Junctional K+ Flux
                f.jty(jj,1) = par.GtK*par.St*(vt-f.vty(jj))/par.F;
            end
            % Ca2+ Activated K+ Channel
            f.jy = par.GK*f.Py*(f.vb-f.vy)/par.F;
            % Volume Change
            f.jw = f.qb - f.Qa;
            
            f.jcalcium = ( ( par.k_f * f.Po + f.J_RyR ) ...
                          * ( ce - var.ca ) - f.J_SERCA ) ...
                              * ( par.v0 / var.w ) - f.jw * var.ca / var.w;
            
            f.jip3 = ( par.V_PLC * ( var.ca^2 ...
                    / ( var.ca^2 + par.K_PLC^2 ) ) ...
                        -(par.V_5K + par.V_3K * ( var.ca^2 / ...
                        ( var.ca^2 + par.K_3K^2 ) ) ) * var.ip ) ...
                              * ( par.v0 / var.w ) - f.jw * var.ip / var.w;
            
            f.jhh = ( H_inf - var.hh ) / f.TAU;
            f.jgg = ( f.g_inf - var.gg ) / par.tau_RyR;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Sorts out which parameter set I want to use
        function param = pardecider(i,neighbours,par)
            v0=[852.745711287631;
                709.618811198231;
                937.841588280235;
                857.349976197479;
                416.039162418089;
                916.352201001994;
                979.554983847748];
            
%             v0 = [770.745711287632; 639.618811198232; 850.841588280235; ...
%                 780.349976197479; 365.101779562562; 835.352201001994; 895.054983847748];
            % Mesh Dependent Parameters (Geometric Cellular Values)
            param.St 	 = par.St{i};
            param.w0 	 = par.w0{i};
            param.v0     = v0(i,1);
            param.Sb     = par.Sb{i};
            param.Sa     = par.Sa{i};
            param.Sa_k   = par.Sa_k{i};
            param.Sb_k   = par.Sb_k{i};
            param.sb_tri = par.sb_tri;
            param.com_tri_ap = par.com_tri_ap;
            for j = 1:length(neighbours)
                param.Sa_ck{i,neighbours(j)} = par.Sa_ck{i,neighbours(j)};
            end
            
            % Nkcc1
            param.aNkcc1 = 2.0712;
            param.a1 = 157.55;
            param.a2 = (2.0096*1e7)*(10^-12);
            param.a3 = 1.0306;
            param.a4 = 1.3852*1e6*(10^-12);
            % NaK-ATPases
            param.r = 1.305e6 * 10^(-9);
            param.alpha1 = 0.641 * 10^(-3);
            param.aNaK = 0.1606;
            
            % Tight Junctions
            param.GtNa = 2.0179e+06;
            param.GtK = 1.4600e+05;
            
            % Ca2+ Activated Cl Channels
            param.GCl = 1.8239e+09;
            param.KCaCC = 0.26;
            param.eta1  = 4.49;
            
            % Ca2+ Activated K Channels
            param.GK = 1.6722e+08;
            param.KCaKC = 0.3;
            param.eta2  = 2;
            
            % Nhe1
            param.G1 = 90.8356;
            param.KNa = 15;
            param.KH = 4.5e-4;
            
            % Ae4
            param.G4 = 1.3361e+03;
            param.KCl = 5.6;
            param.KB = 100;
            
            % Bicarbonate Buffers
            param.GB = 9.7680;
            param.kn = 2.64e4 * 0.012;
            param.kp = 11 * 0.012;
            
            % Constant Concentrations
            param.pHl = 6.81;
            param.pHi = 6.91;
            param.pHe = 7.41;
            param.HCO3l  = 0.000154881661891248;
            param.CO20   = 6.60012219974207;
            param.Ul     = 51.7310374376987;
            param.Cle    = 102.6;
            param.Nae    = 140.2;
            param.Ke     = 5.3;
            param.Hye    = 1e3*10^(-param.pHe);
            param.HCO3e  = 40;
            param.CO2e   = 1.6;
            param.Hl 	 = 1.548816618912483e-04;
            param.CO2l	 = 11.600000000000000;
            param.Hy 	 = 1.230268770812381e-04;
            
            % Thermodynamical Parameters
            param.RTF 	 = 8.314462100000000;
            param.T		 = 310;
            param.F		 = 9.645833650000000e+04;
            param.RTF	 = 26.721207772435513;
            
            % Membrane Permeabilities
            param.La	 = 4.433544803722930;
            param.Lb	 = 1.945808520616380;
            param.Lt	 = 0.059366380058720;
            
            % Calcium Parameters
%             param.k_f = 40; %45;
%             param.tau_max = 10;% 40
%             
%             %%%%%%%%%%%%% VPLC
%             param.V_PLC = 0.0013;
%             %%%%%%%%%%%%%
            param.k_f = 45; %45;
            param.tau_max = 100;% 40

            %%%%%%%%%%%%% VPLC
            param.V_PLC = 0.0013;
            %%%%%%%%%%%%%

            param.k_beta = 0.4;
            param.K_p = 0.2;
            param.K_c = 0.2;
            param.K_h = 0.08;
            param.V_p = 0.9;
            param.K_bar = 0.00001957;
            param.gamma = 5.5;
            param.K_PLC = 0.07;
           
            param.K_3K = 0.4;
            param.V_3K = 0.05;
            param.V_5K = 0.05;
            
            param.K_tau = 0.1;
            param.ct = 3;
            param.V_RyR = 0.04;
            param.K_RyR = 0.2;
            param.K_RyR2 = 0.15;
            param.tau_RyR = 1;
            
        end
        
    end
end












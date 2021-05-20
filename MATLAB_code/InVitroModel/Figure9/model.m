function [dXdt, J] = model(t, X, activity_array, solve_ode)

%% Unpack variables 
DPsi     = X(1); 
sumATP_x = X(2); 
sumADP_x = X(3); 
sumPi_x  = X(4); 
NADH_x   = X(5); 
QH2_x    = X(6); 
cred_i   = X(7); 
sumATP_c = X(8); 
sumADP_c = X(9); 
sumPi_c  = X(10); 

X_DH  = activity_array(1);
X_C1  = activity_array(2);
X_C3  = activity_array(3);
X_C4  = activity_array(4);
X_F   = activity_array(5);
E_ANT = activity_array(6);
E_PiC = activity_array(7);
X_H   = activity_array(8);
X_AtC = activity_array(9);

%%%%%% Constants defining metabolite pools %%%%%%
% Volume fractions and water space fractions
V_c = 1.0;          % buffer volume fraction        % L buffer (L cuvette)^(-1)
V_m = 0.0005;       % mitochondrial volume fraction % L mito (L cuvette)^(-1)
V_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cuvette)^(-1)
W_c = 1.0;          % buffer water space            % L buffer water (L buffer)^(-1)
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)
W_i = 0.1*W_m;     % intermembrane water space     % L IM water (L mito)^(-1)

% Total pool concentrations
NAD_tot = 2.97e-3;  % NAD+ and NADH conc            % mol (L matrix water)^(-1)
Q_tot   = 1.35e-3;  % Q and QH2 conc                % mol (L matrix water)^(-1)
c_tot   = 2.7e-3;   % cytochrome c ox and red conc  % mol (L IM water)^(-1)

% Membrane capacitance
Cm = 3.1e-3;

%%%%%% Set fixed pH, cation concentrations, and O2 partial pressure %%%%%%
% pH
pH_x = 7.40;
pH_c = 7.20;

% K+ concentrations
K_x  = 100e-3;      % mol (L matrix water)^(-1)
K_c  = 140e-3;      % mol (L cyto water)^(-1)

% Mg2+ concentrations
Mg_x = 1.0e-3;        % mol (L matrix water)^(-1)
Mg_c = 1.0e-3;        % mol (L cyto water)^(-1)

% Oxygen partial pressure
PO2 = 25; % mmHg

% Hydrogen ion concentration
H_x = 10^(-pH_x); % mol (L matrix water)^(-1)
H_c = 10^(-pH_c); % mol (L cuvette water)^(-1)

% Oxygen concentration
a_3  = 1.74e-6;   % oxygen solubility in cuvette   % mol (L matrix water * mmHg)^(-1)
O2_x = a_3*PO2;   % mol (L matrix water)^(-1)

% Thermochemical constants
R = 8.314;          % J (mol K)^(-1)
T = 37 + 273.15;    % K
F = 96485;          % C mol^(-1)

% Proton motive force parameters (dimensionless)
n_F  = 8/3;
n_C1 = 4;
n_C3 = 2;
n_C4 = 4;

% Dissociation constants
K_MgATP = 10^(-3.88);
K_HATP  = 10^(-6.33);
K_KATP  = 10^(-1.02);
K_MgADP = 10^(-3.00);
K_HADP  = 10^(-6.26);
K_KADP  = 10^(-0.89);
K_MgPi  = 10^(-1.66);
K_HPi   = 10^(-6.62);
K_KPi   = 10^(-0.42);

% Other concentrations computed from the state variables
NAD_x = NAD_tot - NADH_x;  %% mol (L matrix water)^(-1)
Q_x   = Q_tot - QH2_x;     %% mol (L matrix water)^(-1)
cox_i = c_tot - cred_i;    %% mol (L matrix water)^(-1)

%% Binding polynomials
% Matrix species % mol (L mito water)^(-1)
PATP_x = 1 + H_x/K_HATP + Mg_x/K_MgATP + K_x/K_KATP;
PADP_x = 1 + H_x/K_HADP + Mg_x/K_MgADP + K_x/K_KADP;
PPi_x  = 1 + H_x/K_HPi  + Mg_x/K_MgPi  + K_x/K_KPi;

% Cytosol species % mol (L cuvette water)^(-1)
PATP_c = 1 + H_c/K_HATP + Mg_c/K_MgATP + K_c/K_KATP;
PADP_c = 1 + H_c/K_HADP + Mg_c/K_MgADP + K_c/K_KADP;
PPi_c  = 1 + H_c/K_HPi  + Mg_c/K_MgPi  + K_c/K_KPi;

%% Unbound species
% Matrix species
ATP_x = sumATP_x / PATP_x; % [ATP4-]_x
ADP_x = sumADP_x / PADP_x; % [ADP3-]_x
Pi_x  = sumPi_x  / PPi_x;  % [HPO42-]_x

% Cytosolic species 
ATP_c = sumATP_c / PATP_c; % [ATP4-]_c
ADP_c = sumADP_c / PADP_c; % [ADP3-]_c
Pi_c  = sumPi_c  / PPi_c;  % [HPO42-]_c


%%%%%% NADH Dehydrogenase %%%%%%

% Constants
r      = 6.8385;
k_Pi1  = 4.659e-4;    % mol (L matrix water)^(-1)
k_Pi2  = 6.578e-4;    % mol (L matrix water)^(-1)

% Flux (mol (s * L mito)^(-1))
J_DH = X_DH * (r * NAD_x - NADH_x) * ((1 + sumPi_x / k_Pi1) / (1+sumPi_x / k_Pi2));

%%%%%% Complex I %%%%%%
% NADH_x + Q_x + 5H+_x <-> NAD+_x + QH2_x + 4H+_i + 4DPsi

% Gibbs energy (J mol^(-1))
DrGo_C1 = -109680;
DrGapp_C1 = DrGo_C1 - R * T * log(H_x); 

% Apparent equilibrium constant 
Kapp_C1   = exp( -(DrGapp_C1 + n_C1 * F * DPsi) / (R * T)) * ((H_x / H_c)^n_C1);

% Flux (mol (s * L mito)^(-1))
J_C1 = X_C1 * (Kapp_C1 * NADH_x * Q_x - NAD_x * QH2_x);

%%%%%% Complex III %%%%%%
% QH2_x + 2cuvetteC(ox)3+_i + 2H+_x <-> Q_x + 2cuvetteC(red)2+_i + 4H+_i + 2DPsi

% Gibbs energy (J mol^(-1))
DrGo_C3 = 46690;
DrGapp_C3 = DrGo_C3 + 2 * R * T * log(H_c);

% Apparent equilibrium constant 
Kapp_C3   = exp(-(DrGapp_C3 + n_C3 * F * DPsi) / (R * T)) * (H_x / H_c)^n_C3; 

% Flux (mol (s * L mito)^(-1))
J_C3 = X_C3 * (Kapp_C3 * cox_i^2 * QH2_x - cred_i^2 * Q_x);

%%%%%% Complex IV %%%%%%
% 2 cytoC(red)2+_i + 0.5O2_x + 4H+_x <-> cytoC(ox)3+_x + H2O_x + 2H+_i +2DPsi

% Constant
k_O2 = 1.2e-4;      % mol (L matrix water)^(-1)

% Gibbs energy (J mol^(-1))
DrGo_C4 = -202160;  % J mol^(-1)
DrGapp_C4 = DrGo_C4 - 2 * R * T * log(H_c);

% Apparent equilibrium constant 
Kapp_C4   = exp(-(DrGapp_C4 + n_C4 * F * DPsi) / (R * T)) * (H_x / H_c)^n_C4;

% Flux (mol (s * L mito)^(-1))
J_C4 = X_C4 *(Kapp_C4^0.5 * cred_i * O2_x^0.25 - cox_i) * (1 / (1 + k_O2 / O2_x));

%%%%%% F0F1-ATPase %%%%%%
% ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

% Gibbs energy (J mol^(-1))
DrGo_F = 4990; 
DrGapp_F = DrGo_F + R * T * log( H_x * PATP_x / (PADP_x * PPi_x));

% Apparent equilibrium constant
Kapp_F   = exp( (DrGapp_F + n_F * F * DPsi ) / (R * T)) * (H_c / H_x)^n_F;

% Flux (mol (s * L mito)^(-1))
J_F = X_F * (Kapp_F * sumADP_x * sumPi_x - sumATP_x);

%%%%%% ANT %%%%%%
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

% Constants 
del_D   = 0.0167;
del_T   = 0.0699;
k2o_ANT = 9.54/60;      % s^(-1)
k3o_ANT = 30.05/60;     % s^(-1)
K0o_D   = 38.89e-6;     % mol (L cuvette water)^(-1)
K0o_T   = 56.05e-6;     % mol (L cuvette water)^(-1)
A       = +0.2829;
B       = -0.2086;
C       = +0.2372;

phi = F * DPsi / (R * T);

% Reaction rates
k2_ANT = k2o_ANT * exp((A*(-3) + B*(-4) + C)*phi);
k3_ANT = k3o_ANT * exp((A*(-4) + B*(-3) + C)*phi);

% Dissociation constants
K0_D = K0o_D * exp(3*del_D*phi);
K0_T = K0o_T * exp(4*del_T*phi);

q     = k3_ANT * K0_D * exp(phi) / (k2_ANT * K0_T);
term1 = k2_ANT * ATP_x * ADP_c * q / K0_D;
term2 = k3_ANT * ADP_x * ATP_c / K0_T;
num   = term1 - term2;
den   = (1 + ATP_c/K0_T + ADP_c/K0_D) * (ADP_x + ATP_x * q);

% Flux (mol (s * L mito)^(-1))
J_ANT = E_ANT * num / den;

%%%%%% H+-PI2 cotransporter %%%%%%
% H2PO42-_x + H+_x = H2PO42-_c + H+_c

% Constant
k_PiC = 1.61e-3;  % mol (L cuvette)^(-1)

% H2P04- species
HPi_c = Pi_c * (H_c / K_HPi);
HPi_x = Pi_x * (H_x / K_HPi);

% Flux (mol (s * L mito)^(-1))
J_PiC = E_PiC * (H_c * HPi_c - H_x * HPi_x) / (k_PiC + HPi_c);

%%%%%% H+ leak %%%%%%

% Flux (mol (s * L mito)^(-1))
J_H = X_H * (H_c * exp(phi/2) - H_x * exp(-phi/2));

%%%%%% ATPase %%%%%%
% ATP4- + H2O = ADP3- + PI2- + H+

%Flux (mol (s * L cyto)^(-1))
J_AtC = X_AtC / V_c;

%%%%%% Differential equations (equation 23) %%%%%%
% Membrane potential
dDPsi = (n_C1 * J_C1 + n_C3 * J_C3 + n_C4 * J_C4 - n_F * J_F - J_ANT - J_H) / Cm;

% Matrix species
dATP_x  = (J_F  - J_ANT) / W_x;
dADP_x  = (-J_F + J_ANT) / W_x;
dPi_x   = (-J_F + J_PiC) / W_x;
dNADH_x = (J_DH  - J_C1) / W_x;
dQH2_x  = (J_C1  - J_C3) / W_x;

% IMS species
dcred_i = 2 * (J_C3 - J_C4) / W_i;

% Buffer species
dATP_c = ( V_m2c * J_ANT - J_AtC ) / W_c;
dADP_c = (-V_m2c * J_ANT + J_AtC ) / W_c;
dPi_c  = (-V_m2c * J_PiC + J_AtC)  / W_c;

dXdt = [dDPsi; dATP_x; dADP_x; dPi_x; 
    dNADH_x; dQH2_x; dcred_i; 
    dATP_c; dADP_c; dPi_c];

% Calculate state-dependent quantities after model is solved 
if solve_ode == 1
    dXdt = [dDPsi; dATP_x; dADP_x; dPi_x; 
        dNADH_x; dQH2_x; dcred_i; 
        dATP_c; dADP_c; dPi_c];
else
    J = [PATP_x; PADP_x; PPi_x; PATP_c; PADP_c; PPi_c;
        J_DH; J_C1; J_C3; J_C4; J_F; J_ANT; J_PiC];
end

end 
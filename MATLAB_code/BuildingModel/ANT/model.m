function dXdt = model(t, X, activity_array)

%% Unpack variables 
sumATP_x = X(1); 
sumADP_x = X(2); 
sumPi_x  = X(3);
sumATP_c = X(4); 
sumADP_c = X(5); 

X_F = activity_array(1); 
E_ANT = activity_array(2);  

%% Parameters and constants 

%%%%%% Constants defining metabolite pools %%%%%%
% Volume fractions and water space fractions
V_c = 0.6601;       % cytosol volume fraction       % L cyto (L cell)^(-1)
V_m = 0.2882;       % mitochondrial volume fraction % L mito (L cell)^(-1)
V_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cuvette)^(-1)
W_c = 0.8425;       % cytosol water space           % L cyto water (L cyto)^(-1)
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)

% Membrane potential 
DPsi = 175/1000;

%%%%%% Set fixed pH, cation concentrations, and O2 partial pressure %%%%%%
% pH
pH_x = 7.40;
pH_c = 7.20;

% K+ concentrations
K_x  = 100e-3;     % mol (L matrix water)^(-1)
K_c  = 140e-3;     % mol (L cyto water)^(-1)

% Mg2+ concentrations
Mg_x = 1.0e-3;     % mol (L matrix water)^(-1)
Mg_c = 1.0e-3;     % mol (L cyto water)^(-1)


% Hydrogen ion concentration
H_x = 10^(-pH_x); % mol (L matrix water)^(-1)
H_c = 10^(-pH_c); % mol (L cuvette water)^(-1)

% Thermochemical constants
R = 8.314;          % J (mol K)^(-1)
T = 37 + 273.15;    % K
F = 96485;          % C mol^(-1)

% Proton motive force parameters (dimensionless)
n_F  = 8/3;

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

%% Binding polynomials
% Matrix species % mol (L mito water)^(-1)
PATP_x = 1 + H_x/K_HATP + Mg_x/K_MgATP + K_x/K_KATP;
PADP_x = 1 + H_x/K_HADP + Mg_x/K_MgADP + K_x/K_KADP;
PPi_x  = 1 + H_x/K_HPi  + Mg_x/K_MgPi  + K_x/K_KPi;

% Cytosol species % mol (L cuvette water)^(-1)
PATP_c = 1 + H_c/K_HATP + Mg_c/K_MgATP + K_c/K_KATP;
PADP_c = 1 + H_c/K_HADP + Mg_c/K_MgADP + K_c/K_KADP;

%% Unbound species
% Matrix species
ATP_x = sumATP_x / PATP_x; % [ATP4-]_x
ADP_x = sumADP_x / PADP_x; % [ADP3-]_x

% Cytosol species 
ATP_c = sumATP_c / PATP_c; % [ATP4-]_c
ADP_c = sumADP_c / PADP_c; % [ADP3-]_c

%%%%%% F0F1-ATPase %%%%%%
% ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

% Gibbs energy (J mol^(-1))
DrGo_F   = 4990; 
DrGapp_F = DrGo_F + R * T * log( H_x * PATP_x / (PADP_x * PPi_x));

% Apparent equilibrium constant 
Kapp_F = exp( (DrGapp_F + n_F * F * DPsi ) / (R * T)) * (H_c / H_x)^n_F;

% Flux (mol (s * L mito)^(-1))
J_F = X_F * (Kapp_F * sumADP_x * sumPi_x - sumATP_x);

%%%%%% ANT %%%%%%
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

%Constants
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

% Reaction rates (s^(-1))
k2_ANT = k2o_ANT * exp((A*(-3) + B*(-4) + C)*phi);
k3_ANT = k3o_ANT * exp((A*(-4) + B*(-3) + C)*phi);

% Dissociation constants (M)
K0_D = K0o_D * exp(3*del_D*phi);
K0_T = K0o_T * exp(4*del_T*phi);

q     = k3_ANT * K0_D * exp(phi) / (k2_ANT * K0_T);
term1 = k2_ANT * ATP_x * ADP_c * q / K0_D;
term2 = k3_ANT * ADP_x * ATP_c / K0_T;
num   = term1 - term2;
den   = (1 + ATP_c/K0_T + ADP_c/K0_D) * (ADP_x + ATP_x * q);

% Flux (mol (s * L mito)^(-1))
J_ANT = E_ANT * num / den;

%%%%%% Differential equations (equation 14) %%%%%%
% Matrix species
dATP_x = (J_F  - J_ANT) / W_x;
dADP_x = (-J_F + J_ANT) / W_x;
dPi_x  = 0; %(-J_F ) / W_x

% Cytosol species
dATP_c = ( V_m2c * J_ANT) / W_c;
dADP_c = (-V_m2c * J_ANT) / W_c;

dXdt = [dATP_x; dADP_x; dPi_x; dATP_c; dADP_c];
end
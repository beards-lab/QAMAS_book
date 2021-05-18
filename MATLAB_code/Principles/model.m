function dXdt = model(t, X, DPsi, pH_c)

%% Unpack X state variable

sumATP = X(1); 
sumADP = X(2);
sumPi  = X(3); 

%% Constants and parameters

%Biophysical constants 
R   = 8.314;          % J (mol * K)^(-1)
T   = 310.15;         % K
F   = 96485;          % C mol^(-1)

% F0F1 constants 
n_F    = 8/3;
X_F    = 1000;        % mol (s * L mito)^(-1)
DrGo_F = 4990;        % (J mol^(-1))

% Dissociation constants
K_MgATP = 10^(-3.88);
K_MgADP = 10^(-3.00);
K_MgPi  = 10^(-1.66);
K_HATP  = 10^(-6.33);
K_HADP  = 10^(-6.26);
K_HPi   = 10^(-6.62);
K_KATP  = 10^(-1.02);
K_KADP  = 10^(-0.89);
K_KPi   = 10^(-0.42);

% Environment concentrations 
pH_x = 7.4;           % pH in matrix
H_x  = 10^(-pH_x);    % M 
H_c  = 10^(-pH_c);    % M 
K_x  = 150e-3;        % M 
Mg_x = 1e-3;          % M 

% Volume ratios
W_m = 0.7238;         % (L mito water) (L mito)^(-1)
W_x = 0.9 * W_m;      % (L matrix water) (L mito)^(-1)

%% Equations

% Binding polynomials
P_ATP = 1 + H_x/K_HATP + K_x/K_KATP + Mg_x/K_MgATP; % equation 6
P_ADP = 1 + H_x/K_HADP + K_x/K_KADP + Mg_x/K_MgADP; % equation 7 
P_Pi  = 1 + H_x/K_HPi  + K_x/K_KPi  + Mg_x/K_MgPi;  % equation 8 

% Gibbs energy (equation 9)
DrGapp_F = DrGo_F + R * T * log(H_x * P_ATP / (P_ADP * P_Pi));

% Apparent equilibrium constant 
Kapp_F = exp((DrGapp_F + n_F * F * DPsi)/ (R * T)) * (H_c / H_x) ^ n_F;

% Flux (mol (s * L mito)^(-1))  
J_F = X_F * (Kapp_F * sumADP * sumPi - sumATP);

%%%%%% Differential equations (equation 12) %%%%%%
dATP = J_F / W_x;
dADP = -J_F / W_x;
dPi  = -J_F / W_x;

dXdt = [dATP; dADP; dPi]; 
end 
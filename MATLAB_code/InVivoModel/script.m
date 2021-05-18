

% Volume fractions and water space fractions
V_c = 0.6601;       % cytosol volume fraction       % L cyto (L cell)^(-1)
V_m = 0.2882;       % mitochondrial volume fraction % L mito (L cell)^(-1)
V_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cyto)^(-1)
W_c = 0.8425;      % cytosol water space           % L cyto water (L cyto)^(-1)
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)
W_i = 0.1*W_m;      % intermembrane water space     % L IM water (L mito)^(-1)

% Total pool concentrations
NAD_tot = 2.97e-3;  % NAD+ and NADH conc            % mol (L matrix water)^(-1)
Q_tot   = 1.35e-3;  % Q and QH2 conc                % mol (L matrix water)^(-1)
c_tot   = 2.7e-3;   % cytochrome c ox and red conc  % mol (L IM water)^(-1)

%%%%%% Parameter vector %%%%%% 
X_DH  = 0.1732;
X_C1  = 1.0e4;
X_C3  = 1.0e6;
X_C4  = 0.0125;
X_F   = 1.0e3;
E_ANT = 0.325;
E_PiC = 5.0e6;
X_H   = 1.0e3;
X_CK  = 1e7;
X_AtC = 0.5e-3;

activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_CK, X_AtC];


%%%%%% Initial Conditions %%%%%%
% Membrane Potential
Psi_0 = 175/1000;      % Volts

% Matrix species
sumATP_x_0 = 0.5e-3;  % mol (L matrix water)^(-1)
sumADP_x_0 = 9.5e-3;  % mol (L matrix water)^(-1)
sumPi_x_0 = 0.3e-3;   % mol (L matrix water)^(-1)
NADH_x_0 = 2/3 * NAD_tot;  % mol (L matrix water)^(-1)
QH2_x_0 = 0.1 * Q_tot;   % mol (L matrix water)^(-1)

% IMS species
cred_i_0 = 0.1 * c_tot;

% Cytoplasmic species
%sumATP_c_0 = 9.95e-3  % mol (L cyto water)^(-1)
sumADP_c_0 = 0.05e-3;  % mol s(L cyto water)^(-1)
%sumPi_c_0 = 5.0e-3  % mol (L cyto water)^(-1)

%%%%%% Healthy normal case %%%%%%
TAN = 0.0076; %(M per liter cell)
TEP = 0.0275;  %(M per liter cell)
Cr_tot  = 0.040;    %(M per liter cell)
     
sumATP_c_0 = (TAN - V_m*W_x*(sumATP_x_0 + sumADP_x_0))/(V_c*W_c+V_m*W_i) - sumADP_c_0;
Cr_tot_c   = Cr_tot / (V_c * W_c); % convert to mol (L cyto water)^(-1)
CrP_c_0    = .3 * Cr_tot_c;  % mol (L cyto water)^(-1)
sumPi_c_0  = (TEP-V_m*W_x*(sumATP_x_0 + sumADP_x_0 + sumPi_x_0 ))/(V_c*W_c+V_m*W_i) - 2*sumATP_c_0 - sumADP_c_0 - CrP_c_0; 

X_0 = [Psi_0; sumATP_x_0; sumADP_x_0; sumPi_x_0;
    NADH_x_0; QH2_x_0; cred_i_0;
    sumATP_c_0; sumADP_c_0; sumPi_c_0; CrP_c_0]; 

% range of ATP consumption rates
X_AtC = linspace(0.4e-3,1.2e-3, 60);   % Increase max hydrolysis to find apparent Km.
steady_state = zeros(length(X_0), length(X_AtC));
JO2 = zeros(size(X_AtC));
tspan = [0,10];
options = odeset('MaxStep',0.1); 

% looping through different ATP consumptions states
for i = 1:length(JO2)
   activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_CK, X_AtC(i)];
   % run for long time to acheive steady-state
   steady_state_temp_results = ode15s(@model,tspan, X_0, options, activity_array,1,Cr_tot_c); 
   steady_state(:,i) = steady_state_temp_results.y(:,end);
end 

steady_state = steady_state'; 

DPsi     = steady_state(:,1);
sumATP_x = steady_state(:,2);
sumADP_x = steady_state(:,3);
sumPi_x  = steady_state(:,4);
NADH_x   = steady_state(:,5);
QH2_x    = steady_state(:,6);
cred_i   = steady_state(:,7);
sumATP_c = steady_state(:,8);
sumADP_c = steady_state(:,9);
sumPi_c  = steady_state(:,10);
CrP_c    = steady_state(:,11);

%% Plot normal case
% CrP/ATP ratio

CrP_tot_normal = CrP_c * (V_c * W_c); 
ATP_tot_normal = sumATP_x * V_m * W_x+ sumATP_c *(V_c * W_c + V_m * W_i); 

figure(1)
clf
hold on 
h1 = plot(X_AtC * 1000, CrP_tot_normal./ATP_tot_normal,'b');      % 

% Pi_c
figure(2)
clf
hold on 
plot(X_AtC * 1000, sumPi_c * 1000,'b')     % Pi_c 

%% Heart Failure (HF/TAC) case %%%%%%

%Mean TAC pools
TAN = 0.006976; %(M per liter cell)
TEP = 0.02411;  %(M per liter cell)
Cr_tot = 0.02303; %(M per liter cell)
    
sumATP_c_0 = (TAN - V_m*W_x*(sumATP_x_0 + sumADP_x_0))/(V_c*W_c+V_m*W_i) - sumADP_c_0;
Cr_tot_c   = Cr_tot / (V_c * W_c); % convert to mol (L cyto water)^(-1)
CrP_c_0    = .3 * Cr_tot_c;  % mol (L cyto water)^(-1)
sumPi_c_0  = (TEP-V_m*W_x*(sumATP_x_0 + sumADP_x_0 + sumPi_x_0 ))/(V_c*W_c+V_m*W_i) - 2*sumATP_c_0 - sumADP_c_0 - CrP_c_0; 

X_0 = [Psi_0; sumATP_x_0; sumADP_x_0; sumPi_x_0; 
    NADH_x_0; QH2_x_0; cred_i_0; 
    sumATP_c_0; sumADP_c_0; sumPi_c_0; CrP_c_0]; 

% range of ATP consumption rates
steady_stateHF = zeros(length(X_0), length(X_AtC));
% looping through different ATP consumptions states
for i = 1:length(X_AtC)
   activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_CK, X_AtC(i)];
   % run for long time to acheive steady-state
   steady_state_temp_resultsHF = ode15s(@model,tspan, X_0, options, activity_array,1,Cr_tot_c); 
   steady_stateHF(:,i) = steady_state_temp_resultsHF.y(:,end);
end 

steady_stateHF = steady_stateHF'; 

DPsi     = steady_stateHF(:,1);
sumATP_x = steady_stateHF(:,2);
sumADP_x = steady_stateHF(:,3);
sumPi_x  = steady_stateHF(:,4);
NADH_x   = steady_stateHF(:,5);
QH2_x    = steady_stateHF(:,6);
cred_i   = steady_stateHF(:,7);
sumATP_c = steady_stateHF(:,8);
sumADP_c = steady_stateHF(:,9);
sumPi_c  = steady_stateHF(:,10);
CrP_c    = steady_stateHF(:,11);

%% Plot figures

CrP_tot_HF = CrP_c * (V_c * W_c); 
ATP_tot_HF = sumATP_x * V_m * W_x+ sumATP_c *(V_c * W_c + V_m * W_i); 

% CrP/ATP ratio
figure(1)
h2 = plot(X_AtC * 1000, CrP_tot_HF./ATP_tot_HF,'r');       
ylabel('[CrP]_c/[ATP]_c')
xlabel('ATP consumption rate (mmol s^{-1} (L cell)^{-1})')
xlim([0,1.3])
ylim([0.0,2.5])
legend([h1, h2],'Normal','HF')
set(gca,'FontSize',20)

print -dpng Figure_10a.png
print -depsc2 Figure_10a.eps

% Pi_c
figure(2)
plot(X_AtC * 1000, sumPi_c * 1000, 'r')     
ylabel('[Pi]_c (mM)')
xlabel('ATP consumption rate (mmol s^{-1} (L cell)^{-1})')
xlim([0,1.3])
ylim([0,5])
set(gca,'FontSize',20)

print -dpng Figure_10b.png
print -depsc2 Figure_10b.eps


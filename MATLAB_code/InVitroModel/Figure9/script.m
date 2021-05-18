

%%%%%% Parameter vector %%%%%% 
X_DH  = 0.1732;
X_C1  = 1.0e4;
X_C3  = 1.0e6;
X_C4  = 0.0125;
X_F   = 1.0e3;
E_ANT = 0.325;
E_PiC = 5.0e6;
X_H   = 1.0e3;
X_AtC = 0;

activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC];


%%%%%% Run Low Pi experiments %%%%%%

% Membrane potential
DPsi_0 = 175*1e-3;

% Matrix species
ATP_x_0  = 0.5e-3;
ADP_x_0  = 9.5e-3;
Pi_x_0   = 0.3e-3;
NADH_x_0 = 0.1 * NAD_tot;
QH2_x_0  = 0.1 * Q_tot;

% IMS species
cred_i_0 = 0.1 * c_tot;

% Cytosol species
ATP_c_0 = 5.0e-3;
ADP_c_0 = 0.0e-3;
Pi_c_0  = 1.0e-3;

X_0 = [DPsi_0, ATP_x_0, ADP_x_0, Pi_x_0, ...
    NADH_x_0, QH2_x_0, cred_i_0,...
    ATP_c_0, ADP_c_0, Pi_c_0];

% range of ATP consumption rates
X_AtC = linspace(0,6e-6, 60);   % Increase max hydrolysis to find apparent Km.
steady_state = zeros(length(X_0), length(X_AtC));
JO2 = zeros(size(X_AtC));

% looping through different ATP consumptions states
for i = 1:length(X_AtC)
   activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC(i)]; 
   % run for long time to acheive steady-state
   steady_state_temp_results = ode15s(@model,[0, 3000], X_0, [], activity_array,1); 
   steady_state(:,i) = steady_state_temp_results.y(:,end);
   [~,J] = model(3000,steady_state(:,i), activity_array, 0);
   J_C4 = J(9); % oxygen flux in mol O / sec / (L mito)

   % convert to units of nmol / min / UCS
   % using the conversion factor 0.0012232 mL of mito per UCS
   JO2(i) = J_C4/2 * 60 * 1e9 * 0.0000012232;
   
end 

steady_state = steady_state'; 

DPsi_low_pi   = steady_state(:,1); 
ATP_x_low_pi  = steady_state(:,2); 
ADP_x_low_pi  = steady_state(:,3); 
Pi_x_low_pi   = steady_state(:,4); 
NADH_x_low_pi = steady_state(:,5); 
QH2_x_low_pi  = steady_state(:,6); 
cred_i_low_pi = steady_state(:,7); 
ATP_c_low_pi  = steady_state(:,8); 
ADP_c_low_pi  = steady_state(:,9); 
Pi_c_low_pi   = steady_state(:,10); 


% Low Pi Plotting  Plotting

figure(1)
clf
hold on 
plot(JO2, NADH_x_low_pi/NAD_tot, 'b')   % NADH

figure(2)
clf
hold on 
h1 = plot(JO2, DPsi_low_pi * 1000, 'b');      % Membrane Potential

figure(3)
clf
hold on 
plot(JO2, cred_i_low_pi/c_tot,'b')      % Cytochrome C

figure(4)
clf
hold on 
plot(JO2, ADP_c_low_pi * 1000, 'b')     % ADP


%%%%%% Running High Pi experiments %%%%%%
 
% Membrane potential
DPsi_0 = 175*1e-3;

% Matrix species
ATP_x_0  = 0.5e-3;
ADP_x_0  = 9.5e-3;
Pi_x_0   = 0.3e-3;
NADH_x_0 = 0.1 * NAD_tot;
QH2_x_0  = 0.1 * Q_tot;

% IM species
cred_i_0 = 0.1 * c_tot;

% Cytosol species
ATP_c_0 = 5.0e-3;
ADP_c_0 = 0.0e-3;
Pi_c_0  = 5.0e-3;

X_0 = [DPsi_0, ATP_x_0, ADP_x_0, Pi_x_0, ...
    NADH_x_0, QH2_x_0, cred_i_0, ...
    ATP_c_0, ADP_c_0, Pi_c_0 ];

% range of ATP consumption rates
X_AtC = linspace(0,6e-6, 60);   % Increase max hydrolysis to find apparent Km.
steady_state = zeros(length(X_0), length(X_AtC));
JO2 = zeros(size(X_AtC));

% looping through different ATP consumptions states
for i = 1:length(X_AtC)
   activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC(i)]; 
   % run for long time to acheive steady-state
   steady_state_temp_results = ode15s(@model,[0, 3000], X_0, [], activity_array,1); 
   steady_state(:,i) = steady_state_temp_results.y(:,end);
   [~,J] = model(3000,steady_state(:,i), activity_array, 0);
   J_C4 = J(9); % oxygen flux in mol O / sec / (L mito)

   % convert to units of nmol / min / UCS
   % using the conversion factor 0.0012232 mL of mito per UCS
   JO2(i) = J_C4/2 * 60 * 1e9 * 0.0000012232;
   
end 

steady_state = steady_state'; 

DPsi_hi_pi   = steady_state(:,1); 
ATP_x_hi_pi  = steady_state(:,2); 
ADP_x_hi_pi  = steady_state(:,3); 
Pi_x_hi_pi   = steady_state(:,4); 
NADH_x_hi_pi = steady_state(:,5); 
QH2_x_hi_pi  = steady_state(:,6); 
cred_i_hi_pi = steady_state(:,7); 
ATP_c_hi_pi  = steady_state(:,8); 
ADP_c_hi_pi  = steady_state(:,9); 
Pi_c_hi_pi   = steady_state(:,10); 


% High Pi Plotting  Plotting
figure(1)
plot(JO2, NADH_x_hi_pi/NAD_tot, 'r')   % NADH

figure(2)
h2 = plot(JO2, DPsi_hi_pi * 1000, 'r');      % Membrane Potential

figure(3)
plot(JO2, cred_i_hi_pi/c_tot,'r')      % Cytochrome C

figure(4)
plot(JO2, ADP_c_hi_pi * 1000,  'r')     % ADP
ylim([-0.01,0.8])


%% Plot experimental data

% NADH data
data_NADH = [10.494,      0.7284,                14.509, 0.6642;
       51.927,        0.4879,                53.625, 0.4821;
       71.967,        0.4241,                77.07,    0.3901;
       86.95,          0.4063,                108.57, 0.3688;
       105.62,        0.3964,                144.36, 0.3695];

figure(1)
plot(data_NADH(:,1), data_NADH(:,2), 'bo')
plot(data_NADH(:,3), data_NADH(:,4), 'ro')
ylabel('NADH (normalized)')
set(gca,'FontSize',20)

print -dpng Figure_9a.png
print -depsc2 Figure_9a.eps

% Membrane potential data
data_DPsi = [10.494, 183.2317,                         14.509, 183.4706;
        51.927,       167.5529,                         53.625, 167.9643;
        105.62,       162.3817,                         144.36, 156.2841];
    
figure(2)
h3 = plot(data_DPsi(:,1), data_DPsi(:,2), 'bo');
h4 = plot(data_DPsi(:,3), data_DPsi(:,4), 'ro');
legend([h1,h2,h3,h4], '1 mM Pi Model','5 mM Pi Model','1 mM Pi Exp','5 mM Pi Exp')
ylabel('Membrane Potential \Delta\Psi (mV)')
set(gca,'FontSize',20)

print -dpng Figure_9b.png
print -depsc2 Figure_9b.eps

% Cytochrome C
data_cred = [10.494, 0.1862,                14.509, 0.1931;
        51.927,       0.2658,                53.625, 0.2576;
        105.62,       0.3147,                144.36, 0.2959];

figure(3)
plot(data_cred(:,1), data_cred(:,2), 'bo')
plot(data_cred(:,3), data_cred(:,4), 'ro')
xlabel('OCR (nmol O_2 min^{-1} U CS^{-1})')
ylabel('Cyt c^{2+} (Normalized)')
set(gca,'FontSize',20)

print -dpng Figure_9c.png
print -depsc2 Figure_9c.eps

% Buffer ADP
data_ADP = [10.494, 0, 14.509, 0;
            51.927, 0.0561, 53.625, 0.0049;
            71.967, 0.1827, 77.07,    0.0578;
            86.95, 0.3589, 108.57, 0.1648;
            105.62, 0.6622, 144.36, 0.4716];

figure(4)
plot(data_ADP(:, 1), data_ADP(:, 2),  'bo')
plot(data_ADP(:, 3), data_ADP(:, 4),  'ro')
xlabel('OCR (nmol O_2 min^{-1} U CS^{-1})')
ylabel('Buffer ADP (Normalized)')
set(gca,'FontSize',20)

print -dpng Figure_9d.png
print -depsc2 Figure_9d.eps

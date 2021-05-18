

%% Figure 7 

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
X_AtC = 0;

activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC]';

%%%%%% Initial Conditions %%%%%%
% Membrane Potential
DPsi_0 = 175/1000;      % V

% Matrix species
sumATP_x_0 = 0.5e-3;         % mol (L matrix water)^(-1)
sumADP_x_0 = 9.5e-3;         % mol (L matrix water)^(-1)
sumPi_x_0  = 1.0e-3;         % mol (L matrix water)^(-1)
NADH_x_0   = 2/3 * NAD_tot;  % mol (L matrix water)^(-1)
QH2_x_0    = 0.1 * Q_tot;    % mol (L matrix water)^(-1)

% IMS species
cred_i_0 = 0.1 * c_tot; % mol (L IMS water)^(-1)

% Cytosolic species
sumATP_c_0 = 0;       % mol (L cyto water)^(-1)
sumADP_c_0 = 10e-3;   % mol (L cyto water)^(-1)
sumPi_c_0  = 10e-3;   % mol (L cyto water)^(-1)

X_0 = [DPsi_0, sumATP_x_0, sumADP_x_0, sumPi_x_0, ...
    NADH_x_0, QH2_x_0, cred_i_0, ...
    sumATP_c_0, sumADP_c_0, sumPi_c_0]; 
   
    
% Time vector 
t = linspace(0,5,100);
    
% Solve ODE
results = ode15s(@model, [0, 5], X_0, [], activity_array, 1); 
results = deval(t,results); 
results = results'; 

DPsi     = results(:,1);
sumATP_x = results(:,2);
sumADP_x = results(:,3);
sumPi_x  = results(:,4);
NADH_x   = results(:,5);
QH2_x    = results(:,6);
cred_i   = results(:,7);
sumATP_c = results(:,8);
sumADP_c = results(:,9);
sumPi_c  = results(:,10);

% Plot figures 
figure(1)
clf
hold on 
h1 = plot(t, sumATP_x*1000);
h2 = plot(t, sumADP_x*1000);
h3 = plot(t, sumPi_x*1000);
legend([h1 h2 h3], '[$\Sigma$ATP]$_x$', '[$\Sigma$ADP]$_x$', '[$\Sigma$Pi]$_x$', ... 
    'interpreter','latex','location','east')
xlabel('Time (s)') 
ylabel('Concentration (mM)')
ylim([-.5,10.5])
set(gca,'FontSize',20)

print -dpng Figure_7a.png
print -depsc2 Figure_7a.eps

figure(2)
clf
hold on
h4 = plot(t, sumATP_c*1000);
h5 = plot(t, sumADP_c*1000);
h6 = plot(t, sumPi_c*1000);
legend([h4 h5 h6], '[$\Sigma$ATP]$_c$', '[$\Sigma$ADP]$_c$', '[$\Sigma$Pi]$_c$', ... 
    'interpreter','latex','location','east')
xlabel('Time (s)')
ylabel('Concentration (mM)')
ylim([-.5,10.5])
set(gca,'FontSize',20)

print -dpng Figure_7b.png
print -depsc2 Figure_7b.eps


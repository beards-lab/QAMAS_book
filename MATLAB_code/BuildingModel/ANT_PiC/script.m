

%% Figure 5

% Membrane potential 
DPsi = 175/1000;

%%%%%% Parameter vector %%%%%% 
X_F   = 100;        % Synthase activity 
E_ANT = 0.325;      % Nucleotide transporter activity 
E_PiC = 5.0e6;      % Phosphate transporter activity 

activity_array = [X_F, E_ANT, E_PiC]';     

%%%%%% Initial Conditions %%%%%%
% Matrix species
sumATP_x_0 = 0.5e-3;  % mol (L matrix water)^(-1)
sumADP_x_0 = 9.5e-3;  % mol (L matrix water)^(-1)
sumPi_x_0  = 1e-3;  % mol (L matrix water)^(-1)

% Cytosolic species
sumATP_c_0 = 0;       % mol (L cyto water)^(-1)
sumADP_c_0 = 10e-3;   % mol (L cyto water)^(-1)
sumPi_c_0  = 10e-3;   % mol (L cyto water)^(-1)

X_0 = [sumATP_x_0, sumADP_x_0, sumPi_x_0, sumATP_c_0, sumADP_c_0, sumPi_c_0]';


% Solve ODE
t = linspace(0,2,100);
results = ode15s(@model, [0, 2], X_0, [], activity_array, DPsi); 
results = deval(t,results);
results = results'; 

sumATP_x = results(:,1);
sumADP_x = results(:,2); 
sumPi_x  = results(:,3); 
sumATP_c = results(:,4); 
sumADP_c = results(:,5); 
sumPi_c  = results(:,6);

% Plot figures 
figure(1)
clf
hold on 
h1 = plot(t, sumATP_x*1000);
h2 = plot(t, sumADP_x*1000);
h3 = plot(t, sumPi_x*1000);
legend([h1 h2 h3], '[$\Sigma$ATP]$_x$', '[$\Sigma$ADP]$_x$', '[$\Sigma$Pi]$_x$', ... 
    'interpreter','latex','location','east')
ylim([-.5,10.5])
xlim([0,2])
xticks([0,1,2])
set(gca,'xtick',xticks)
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',20)

print -dpng Figure_5a.png
print -depsc2 Figure_5a.eps

figure(2)
clf
hold on 
h4 = plot(t, sumATP_c*1000);
h5 = plot(t, sumADP_c*1000);
h6 = plot(t, sumPi_c*1000);
ylim([-0.5,10.5])
xlim([0,2])
xticks([0,1,2])
set(gca,'xtick',xticks)
xlabel('Time (s)')
legend([h4 h5 h6], '[$\Sigma$ATP]$_c$', '[$\Sigma$ADP]$_c$', '[$\Sigma$Pi]$_c$', ... 
    'interpreter','latex','location','east')
set(gca,'FontSize',20)

print -dpng Figure_5b.png
print -depsc2 Figure_5b.eps

%% Figure 6 

%%% Simulate over a range of Membrane potential from 100 mV to 250 mV %%%

% Define array to iterate over
membrane_potential = linspace(100,250);    % mV

% Define arrays to store steady state results 
ATP_x_steady = zeros(size(membrane_potential));
ADP_x_steady = zeros(size(membrane_potential));
Pi_x_steady  = zeros(size(membrane_potential));

ATP_c_steady = zeros(size(membrane_potential));
ADP_c_steady = zeros(size(membrane_potential));
Pi_c_steady  = zeros(size(membrane_potential));

% Iterate through range of membrane potentials 
for i = 1:length(membrane_potential)
    DPsi = membrane_potential(i) / 1000;      % convert to V
    temp_results = ode15s(@model, [0, 200], X_0, [], activity_array, DPsi); 
    ATP_x_steady(i) = temp_results.y(1,end) * 1000; 
    ADP_x_steady(i) = temp_results.y(2,end) * 1000;
    Pi_x_steady(i)  = temp_results.y(3,end) * 1000;
    ATP_c_steady(i) = temp_results.y(4,end) * 1000;
    ADP_c_steady(i) = temp_results.y(5,end) * 1000;
    Pi_c_steady(i)  = temp_results.y(6,end) * 1000;
end 

% Plot figures 
figure(3)
clf
hold on 
h1 = plot(membrane_potential, ATP_x_steady);
h2 = plot(membrane_potential, ADP_x_steady);
h3 = plot(membrane_potential, Pi_x_steady);
legend([h1 h2 h3], '[$\Sigma$ATP]$_x$', '[$\Sigma$ADP]$_x$', '[$\Sigma$Pi]$_x$', ... 
    'interpreter','latex','location','east')
xlabel('Membrane potential (mV)')
ylabel('Concentration (mM)')
xlim([100, 250])
ylim([-0.5,13])
set(gca,'FontSize',20)

print -dpng Figure_6a.png
print -depsc2 Figure_6a.eps

figure(4) 
clf
hold on
h4 = plot(membrane_potential, ATP_c_steady);
h5 = plot(membrane_potential, ADP_c_steady);
h6 = plot(membrane_potential, Pi_c_steady);
legend([h4 h5 h6], '[$\Sigma$ATP]$_c$', '[$\Sigma$ADP]$_c$', '[$\Sigma$Pi]$_c$', ... 
    'interpreter','latex','location','east')
xlabel('Membrane potential (mV)')
ylabel('Concentration (mM)')
xlim([100, 250])
ylim([-0.5,13])
set(gca,'FontSize',20)

print -dpng Figure_6b.png
print -depsc2 Figure_6b.eps



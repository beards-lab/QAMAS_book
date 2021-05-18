

%% Figure 2

% Simple steady state simulation at 175 mV membrane potential 

% Initial conditions (M)
sumATP_0 = 0.5e-3;
sumADP_0 = 9.5e-3;
sumPi_0  = 1e-3;

X_0 = [sumATP_0, sumADP_0, sumPi_0];

% Inputs  
DPsi = 175e-3; % Constant membrane potential (V)
pH_c = 7.2;    % IMS/buffer pH 

solutions = ode15s(@model, [0 1], X_0, [], DPsi, pH_c);
t = solutions.x;
results = solutions.y; 
results = results * 1000;

% Plot figure 
figure(1)
clf
hold on 
h1 = plot(t, results(1,:));
h2 = plot(t, results(2,:));
h3 = plot(t, results(3,:));
legend([h1 h2 h3],'[$\Sigma$ATP]$_x$','[$\Sigma$ATP]$_x$','[$\Sigma$Pi]$_x$',...
    'interpreter','latex','location','east')
xlabel('Time (s)')
ylabel('Concentration (mM)')
ylim([0, 10])
set(gca,'FontSize',20)

print -dpng Figure_2.png
print -depsc2 Figure_2.eps


%% Figure 3 

%%% Simulate over a range of Membrane potential from 100 mV to 250 mV %%%

% Define array to iterate over
membrane_potential = linspace(100,250);    % mV

% Constant external pH
pH_c = 7.2; % IMS/buffer pH

% Define arrays to store steady state results 
ATP_steady_DPsi = zeros(size(membrane_potential));
ADP_steady_DPsi = zeros(size(membrane_potential));
Pi_steady_DPsi  = zeros(size(membrane_potential));

% Iterate through range of membrane potentials 
for i = 1: length(membrane_potential)
    DPsi = membrane_potential(i) / 1000;      % convert to V
    temp_results = ode15s(@model, [0 5], X_0, [], DPsi, pH_c); 
    ATP_steady_DPsi(i) = temp_results.y(1,end) * 1000; 
    ADP_steady_DPsi(i) = temp_results.y(2,end) * 1000;
    Pi_steady_DPsi(i)  = temp_results.y(3,end) * 1000;
end 
    
% Concentration vs DPsi
figure(3)
clf 
hold on 
h1 = plot(membrane_potential, ATP_steady_DPsi); 
h2 = plot(membrane_potential, ADP_steady_DPsi);
h3 = plot(membrane_potential, Pi_steady_DPsi);
legend([h1 h2 h3],'[$\Sigma$ATP]$_x$','[$\Sigma$ATP]$_x$','[$\Sigma$Pi]$_x$',...
    'interpreter','latex','location','east')
xlabel('Membrane potential (mV)')
ylabel('Concentration (mM)')
xlim([100, 250])
ylim([0, 10])
set(gca,'FontSize',20)  

print -dpng Figure_3.png
print -depsc2 Figure_3.eps
    
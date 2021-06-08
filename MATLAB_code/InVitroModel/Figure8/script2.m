

%% Figure 8 

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

activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC];

%%%%%% Initial Conditions %%%%%%
% Membrane Potential
DPsi_0 = 175/1000;      % V

% Matrix species
sumATP_x_0 = 0.5e-3;        % mol (L matrix water)^(-1)
sumADP_x_0 = 9.5e-3;        % mol (L matrix water)^(-1)
sumPi_x_0  = 0.3e-3;        % mol (L matrix water)^(-1)
NADH_x_0   = 0;             % mol (L matrix water)^(-1)
QH2_x_0    = 0.1 * Q_tot;   % mol (L matrix water)^(-1)

% IMS species
cred_i_0 = 0.1 * c_tot; % mol (L IMS water)^(-1)

% Cytosolic species
sumATP_c_0 = 0;       % mol (L cyto water)^(-1)
sumADP_c_0 = 0;       % mol (L cyto water)^(-1)
sumPi_c_0  = 5.0e-3;  % mol (L cyto water)^(-1)

X_0 = [DPsi_0, sumATP_x_0, sumADP_x_0, sumPi_x_0, ...
    NADH_x_0, QH2_x_0, cred_i_0, ...
    sumATP_c_0, sumADP_c_0, sumPi_c_0]'; 

%%% Four State Model %%%
% State 1 - no substrates
t_1 = linspace(0,25,25*2);
X_AtC = 0;
X_DH  = 0.001; % Kept non-zero for solver stability
activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC]';

state_1_results = ode15s(@model,[0 25],X_0,[],activity_array,1); 
state_1_results = deval(t_1,state_1_results); 

% State 2 - Add substrate (i.e. turn on X_DH)
t_2 = linspace(25,75,(75 - 25)*2);
X_DH = 0.0866 * 2;
activity_array = [X_DH, X_C1, X_C3, X_C4, X_F, E_ANT, E_PiC, X_H, X_AtC]; 

state_2_results = ode15s(@model,[25,75],state_1_results(:,end),[],activity_array,1); 
state_2_results = deval(t_2,state_2_results); 

% State 3 and 4 - Add ADP
t_3 = linspace(75,200,(200 - 75)*2);
state_2_results(9, end) = 0.375e-3; % Molar
state_3_results = ode15s(@model,[75,200],state_2_results(:,end),[],activity_array,1); 
state_3_results = deval(t_3,state_3_results); 

% Concatenate Results
% Note: need to prevent duplicating points
all_results = [state_1_results(:,1:end-1), state_2_results(:,1:end-1), state_3_results]'; 
t = [t_1(1:end-1), t_2(1:end-1), t_3]'; 

DPsi     = all_results(:,1); 
sumATP_x = all_results(:,2); 
sumADP_x = all_results(:,3); 
sumPi_x  = all_results(:,4); 
NADH_x   = all_results(:,5); 
QH2_x    = all_results(:,6); 
cred_i   = all_results(:,7); 
sumATP_c = all_results(:,8); 
sumADP_c = all_results(:,9); 
sumPi_c  = all_results(:,10); 


% Calculate complex IV Flux
J_C4 = zeros(size(t));
J_F  = zeros(size(t)); 
J_H  = zeros(size(t)); 
for i = 1:length(J_C4)
    [~,J] = model(t(i), all_results(i,:), activity_array, 0);
    J_C4(i) = J(10);
    J_F(i)  = J(11); 
    J_H(i)  = J(12); 
end 

% Convert complex IV flux to oxygen flux in nmol O2 / U Citrate Synthase
JO2 = J_C4/2 * 60 * 1e9 * 0.0000012232; 

figure(10)
clf
hold on 
plot(t,J_F * 60 * 1e9 * 0.0000012232,'b')
plot(t,J_H * 60 * 1e9 * 0.0000012232,'r')
plot(t,J_C4/2 * 60 * 1e9 * 0.0000012232,'g')
ylim([-200, 1000])
legend('J_F','J_H','J_{C4}')
set(gca,'FontSize',20)

%% Plot figures 


figure(1)
clf 
hold on 
plot([-5,205],[0,0],'k')
plot([0,0],[-5,205],'k:', 'linewidth',1.75)
plot([25,25],[-5,205],'k:', 'linewidth',1.75)
plot([75,75],[-5,205],'k:', 'linewidth',1.75)
plot([150,150],[-5,205],'k:', 'linewidth',1.75)
plot(t, JO2)
ylim([-5,205])
xlim([-5,205])
text(2, 190, 'State 1','fontsize',15)
text(37.5, 190, 'State 2', 'fontsize',15)
text(100, 190, 'State 3','fontsize',15)
text(165, 190, 'State 4','fontsize',15)
xlabel('Times (s)')
ylabel('OCR (nmol O_2 min^{-1} U CS^{-1})')
set(gca,'FontSize',20)

print -dpng Figure_8a.png
print -depsc2 Figure_8a.eps

figure(2)
clf
hold on
plot([-5,1],[0,0],'k')
plot([0,0],[-5,1],'k:', 'linewidth',1.75)
plot([25,25],[-5,1],'k:', 'linewidth',1.75)
plot([75,75],[-5,1],'k:', 'linewidth',1.75)
plot([150,150],[-5,1],'k:', 'linewidth',1.75)
plot(t, NADH_x/NAD_tot,'r')
ylim([-0.05,1])
xlim([-5,205])
text(2, 0.95, 'State 1','fontsize',15)
text(37.5, 0.95, 'State 2', 'fontsize',15)
text(100, 0.95, 'State 3','fontsize',15)
text(165, 0.95, 'State 4','fontsize',15)
xlabel('Times (s)')
ylabel('[NADH]_x (mM)')
set(gca,'FontSize',20)

print -dpng Figure_8b.png
print -depsc2 Figure_8b.eps

figure(3)
clf 
hold on 
plot([-5,210],[0,0],'k')
plot([0,0],[-5,210],'k:', 'linewidth',1.75)
plot([25,25],[-5,210],'k:', 'linewidth',1.75)
plot([75,75],[-5,210],'k:', 'linewidth',1.75)
plot([150,150],[-5,210],'k:', 'linewidth',1.75)
plot(t, DPsi*1000,'g')
ylim([75,210])
xlim([-5,205])
text(2, 200, 'State 1','fontsize',15)
text(37.5, 200, 'State 2', 'fontsize',15)
text(100, 200, 'State 3','fontsize',15)
text(165, 200, 'State 4','fontsize',15)
xlabel('Times (s)')
ylabel('\Delta\Psi (mV)')
set(gca,'FontSize',20)

print -dpng Figure_8c.png
print -depsc2 Figure_8c.eps

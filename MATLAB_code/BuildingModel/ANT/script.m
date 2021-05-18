


%%%%%% Parameter vector %%%%%% 
X_F   = 1000;      % Synthase activity 
E_ANT = 0.325;     % Nucleotide transporter activity 
activity_array = [X_F, E_ANT];     % Note: This array will be larger in the future parts

%%%%%% Initial Conditions %%%%%%
% Matrix species
sumATP_x_0 = 0.5e-3;  % mol (L matrix water)^(-1)
sumADP_x_0 = 9.5e-3;  % mol (L matrix water)^(-1)
sumPi_x_0  = 1e-3;  % mol (L matrix water)^(-1)

% Cytoplasmic species
sumATP_c_0 = 0; %9.95e-3        % mol (L cyto water)^(-1)
sumADP_c_0 = 10e-3; %0.05e-3        % mol (L cyto water)^(-1)

X_0 = [sumATP_x_0, sumADP_x_0, sumPi_x_0, sumATP_c_0, sumADP_c_0]';


% Solve ODE
results = ode15s(@model,[0, 2],X_0,[],activity_array); 
t = results.x;
sumATP_x = results.y(1,:);
sumADP_x = results.y(2,:);
sumPi_x  = results.y(3,:);
sumATP_c = results.y(4,:);
sumADP_c = results.y(5,:);

% Plot figures 
figure(1)
clf
hold on 
h1 = plot(t, sumATP_x*1000);
h2 = plot(t, sumADP_x*1000);
h3 = plot(t, sumPi_x*1000);
legend([h1 h2 h3],'[$\Sigma$ATP]$_x$','[$\Sigma$ADP]$_x$','[$\Sigma$Pi]$_x$',...
    'interpreter','latex','location','east')
ylim([-.5,10.5])
xlabel('Time (s)')
xticks([0,1,2])
ylabel('Concentration (mM)')
set(gca,'FontSize',20)

print -dpng Figure_4a.png
print -depsc2 Figure_4a.eps

figure(2)
clf
hold on 
h4 = plot(t, sumATP_c*1000);
h5 = plot(t, sumADP_c*1000);
ylim([-0.5,10.5])
xticks([0,1,2])
set(gca,'xtick',xticks)
legend([h4 h5],'[$\Sigma$ATP]$_c$','[$\Sigma$ADP]$_c$',...
    'interpreter','latex','location','east')
xlabel('Time (s)')
set(gca,'FontSize',20)

print -dpng Figure_4b.png
print -depsc2 Figure_4b.eps

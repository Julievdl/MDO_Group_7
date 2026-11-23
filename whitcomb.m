clear; clc;

[file, folder] = uigetfile('*.dat', 'Select withcomb135.dat');
data = readmatrix(fullfile(folder, file));

xi = data(:,1);
zeta = data(:,2);

figure
plot(xi, zeta, 'o-')
axis equal
grid on
xlabel('x/c')
ylabel('z/c')
title('Whitcomb airfoil coordinates')

[~, iLE] = min(xi);                                                         % finds the position in the vector 'xi' with the smallest value, to identify where the LE starts

% Upper surface: reverse to run from leading to trailing edge
xi_u = flipud(xi(1:iLE));
zeta_u = flipud(zeta(1:iLE));

% Lower surface: already runs from leading to trailing edge
xi_l = xi(iLE:end);
zeta_l = zeta(iLE:end);

% Separate trailing-edge heights
zetaTE_u = zeta_u(end);
zetaTE_l = zeta_l(end);

%%
% Upper code
n = 5;
Cu = sqrt(xi_u).*(1-xi_u);
Bu = zeros(length(xi_u), n+1);

for i = 0:n
    Bu(:,i+1) = nchoosek(n,i).*xi_u.^i.*(1-xi_u).^(n-i);
end

surface_u = @(Au) Cu.*(Bu*Au) + xi_u.*zetaTE_u;                             % note that the second term comes from the fact that this airfoil has a finite TE thickness and is therefore not equal to 0,
residual_u = @(Au) zeta_u - surface_u(Au);                                              % because remember for epsilon = 1, C = 0, and zeta = 0 in case the extra term is not there

objective_u = @(Au) sum(residual_u(Au).^2); 

Au0 = 0.2*ones(6,1);

[Au_fit, Ju] = fmincon( ...
    objective_u, Au0, [], [], [], [], [], []);

Au_fit;
Ju;
max_error_u = max(abs(surface_u(Au_fit)-zeta_u));

%%
% Lower code
n = 5;
Cl = sqrt(xi_l).*(1-xi_l);
Bu = zeros(length(xi_l), n+1);

for i = 0:n
    Bl(:,i+1) = nchoosek(n,i).*xi_l.^i.*(1-xi_l).^(n-i);
end

surface_l = @(Al) Cl.*(Bl*Al) + xi_l.*zetaTE_l;
residual_l = @(Al) zeta_l - surface_l(Al);  

objective_l = @(Al) sum(residual_l(Al).^2); 

Al0 = -0.2*ones(6,1);

[Al_fit, Jl] = fmincon( ...
    objective_l, Al0, [], [], [], [], [], []);

Al_fit;
Jl;
max_error_l = max(abs(surface_l(Al_fit)-zeta_l));

%%
% Plot
figure
plot(xi_u,zeta_u,'go', xi_l,zeta_l,'bo')
hold on
plot(xi_u,surface_u(Au_fit),'r-')
plot(xi_l,surface_l(Al_fit),'m-')
axis equal
grid on
xlabel('x/c')
ylabel('z/c')
legend('Upper target','Lower target','Upper CST','Lower CST')

%%
% wing modification
xi_check = linspace(0,1,501)';
Ccheck = sqrt(xi_check).*(1-xi_check);

Bcheck = zeros(length(xi_check), n+1);
for i = 0:n
    Bcheck(:,i+1) = nchoosek(n,i) .* xi_check.^i ...
        .* (1-xi_check).^(n-i);
end

zu = Ccheck.*(Bcheck*Au_fit) + xi_check.*zetaTE_u;
zl = Ccheck.*(Bcheck*Al_fit) + xi_check.*zetaTE_l;

thickness = zu-zl;
[tc_old, imax] = max(thickness);
xi_at_max = xi_check(imax);                                                 % shows the position 'xi' where thickness is maximum

min_thickness = min(thickness);                                             % should be zero (because of TE and LE), if negative, geometry is wrong
k = 0.10/tc_old;

%%
% Camber code
Ac = (Au_fit + Al_fit)/2;                                                   % normalized camber
At = Au_fit - Al_fit;                                                       % normalized thickness

Au_new = Ac + k*At/2;                                                       % calculates the new 'A' coefficients after a thickness modification based on the variable 'k' for the upper surface
Al_new = Ac - k*At/2;                                                       % calculates the new 'A' coefficients after a thickness modification based on the variable 'k' for the lower surface

TE_c = (zetaTE_u + zetaTE_l)/2;                                             % Mean trailing-edge height
TE_t = zetaTE_u - zetaTE_l;                                                 % Trailing-edge thickness

zetaTE_u_new = TE_c + k*TE_t/2;
zetaTE_l_new = TE_c - k*TE_t/2;

zu_new = Ccheck.*(Bcheck*Au_new) + xi_check.*zetaTE_u_new;                  % new height coordinates after modification for upper surface
zl_new = Ccheck.*(Bcheck*Al_new) + xi_check.*zetaTE_l_new;                  % new height coordinates after modification for lower surface
thickness_new = zu_new - zl_new;

tc_new = max(zu_new-zl_new);                                                % new t/c ratio after modification, which for the assignment should be 0.1
camber_error = max(abs( ...                                                 % camber error, which should be approximately equal to 0, as modifying t/c doesn't change camber location                                
    (zu_new+zl_new)/2 - (zu+zl)/2));

%%
% Tank volume calculation
xi_tank = linspace(0.15,0.60,301)';
C_tank = sqrt(xi_tank).*(1-xi_tank);

B_tank = zeros(length(xi_tank),6);
for i = 0:5
    B_tank(:,i+1) = nchoosek(5,i).*xi_tank.^i ...
        .*(1-xi_tank).^(5-i);
end

zu_tank = C_tank.*(B_tank*Au_new) ...
    + xi_tank.*zetaTE_u_new;

zl_tank = C_tank.*(B_tank*Al_new) ...
    + xi_tank.*zetaTE_l_new;

I_tank = trapz(xi_tank, zu_tank-zl_tank);
%%
save('whitcomb_CST_baseline.mat', ...
    'Au_new','Al_new','zetaTE_u_new','zetaTE_l_new');

save('whitcomb_CST_baseline.mat','I_tank','-append');
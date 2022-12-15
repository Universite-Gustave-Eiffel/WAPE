% Author : Bill Kayser
% https://orcid.org/0000-0002-3403-2540
% https://www.researchgate.net/profile/Bill-Kayser

% /*
%  * Copyright 2020-2021 UMRAE
%  *
%  * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the 
%  * European Commission - subsequent versions of the EUPL (the "Licence");
%  * 
%  * You may not use this work except in compliance with the Licence.
%  * You may obtain a copy of the Licence at:
%  *
%  * https://joinup.ec.europa.eu/software/page/eupl5
%  *
%  * Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
%  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  * See the Licence for the specific language governing permissions and limitations under the Licence.
%  */

clear all
close all
clc

tic
%% Source parameters

% [variables]
freq = 250;              % source frequency (Hz)
hS = 35;                % source height (m)

%% Medium properties

% I. Atmospheric parameters [variables]
T = 10;                    % atmospheric temperature at the surface (°C)
Tlog = 0.2;                 % temperature coefficient for log temperature profil -0.5<Tlog<0.5. if 0 no profil

shear_exp = 0.15;    % wind shear exponents (scalar) for power law wind profil see [van den Berg 2008]
v_ref = 5;                 % wind speed (m/s) measured at z_ref height, for power-law wind profile
z_ref = 80;               % reference height for wind speed v_ref

theta = 0;                % propagation angle with respect to the source (0° : downwind, 180° : upwind)

f_turb_ind = true;    % logical 'true' or 'false' to account for turbulence or not
gamT = 0;                % turbulence strength used in "fast_turb" function. gamT = (CT/T0)^2 + 22/3 (Cv/c0)^2, e.g. [van Renterghem et al 2022]

% II.General constant [fixed values] %
Patm = 101300;                    % atmospheric pressure (Pa)
Rg = 286.6896552;               % perfect gas constant (Rg=R/M=8.314E3/29)
gam = 1.41;                          % heat ratio
hr = 0.8;                                % relative humidity of air 0<hr<1
temp0 = 273.15+T;               % temperature (K)

% III. Acoustic parameters [fixed values]
c0 = sqrt(gam*Rg*temp0);     % adiabatic sound speed (m.s-1)
lambda = c0/freq;                   % wavelength (m)
k0 = 2*pi*freq/c0;                    % adiabatic wavenumber (rad/m)

%% Spatial grid 2D(x,z)

% [variables]
dim_x = 3000;              % horizontal dimension x of the domain (m)
dim_z = 300;                % vertical dimension z of the domain (m)

% absorbing layer parameters [fixed values]
haut_a = 0.8; % where the absorbing layer starts at the top of the domain, according to z axis (0.8 ==> start at 80% of z axis)
coeff_a = 10; % damping rate

% spatial discretization [variables, don't go too low (<8;6)]
discrx = 10;               % discretization x = lambda/discrx
discrz = 10;               % discretization z = lambda/discrz

% numerical parameters [fixed values]
delta_x = lambda/discrx;       % spatial discretization step for calculation (m)
delta_z = lambda/discrz;       % spatial discretization step for calculation (m)

nb_x = floor(dim_x/delta_x)+1;   % number of point in the domain along x
nb_z = floor(dim_z/delta_z)+1;   % number of point in the domain along z

% z vertical axis (m) [fixed values]
z = (0:nb_z-1)*delta_z;
z = z';

% x horizontal axis (m) [fixed values]
x = 0:delta_x:dim_x;

stock_x = 0.5;                % spatial discretisation step along x, for saving data (m)
stock_z = 0.5;                % spatial discretisation step along z, for saving data (m)

z_r = 1.5;                       % receiver height for plotting (m)
iz_r = round(z_r/stock_z); % index of receiver height in z vector

%% Ground properties

% [variables]
hv = 0;                                  % vegetation height (m), it affects the shape of atmospheric profils
z0 = 0.13*hv + 0.00001;       % atmospheric roughness length (m), can't be null
d = 0.66*hv;                          % displacement height of flux profiles (m), it's directly linked to vegetation height

lc = 0 ;        %[0.05-1] (m)        % correlation length (ground roughness parameter), if 0 : no ground rugositiy 
sigmah = 0; %[0.01-0.05] (m)  % standard deviation of roughness height (ground roughness parameter)
sigma = 1000;                                % airflow resistivity of the ground
h = 0;                                       % thickness of ground surface layer (m), if h = 0 no layer

% [fixed values]
incidence = pi/2;             % wave incidence angle on the ground, pi/2 : incidence rasante

%% Ground admittance calculation

[Beta]= Miki(freq,k0,sigma,h);

% ground roughness effect on admittance
[Beta_rough] = roughness(lc,sigmah,k0,incidence);
Beta_eff = Beta + Beta_rough;

% reflexion coefficient of ground
Rp = (1-Beta_eff)/(1+Beta_eff);

%% Starting pressure field

U0 = sourceSalomons_order2(z,hS,nb_z,k0,Rp); % using Salomons source function

%% Wind profile

[Mx,epsilon0,epsilon1] = Mach(v_ref,z_ref,z,z0,d,shear_exp,c0,Tlog,temp0,Rg,gam,theta);

%% Numerical variables

gamma_x = 1./sqrt(1-Mx.^2);                    % Lorentz factor, see eq 32 from [Ostashev et al 2020]
tau_bar1 = (Mx.*gamma_x.^2).*(sqrt(1+epsilon1) - Mx); % operator used in following development, see eq 32 from [Ostashev et al 2020]

%% WAPE main calculation

% calcul of velocity potential
[Phi] = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,epsilon0,epsilon1,Beta,haut_a,coeff_a);

%% Pressure calculation

% the pressure is given from eq 84 of [Ostashev et al 2020]
% we need phi(n-1), phi(n) and phi(n+1) to calculate p(n)
n =  2 : size(Phi,2)-1;
Phi_n_i = Phi(:,n);
Phi_prv = Phi(:,n-1);
Phi_nxt = Phi(:,n+1);
x_pressure = x(n); % horizontal axis for pressure is shorter than x for phi (due to n-1, n and n+1 terms)

% 2D pressure field -3dB when x = *2
pp2D = zeros(nb_z,nb_x-1);
pp2D(:,n) = abs(exp(1i*k0*x_pressure).* ((1-Mx).*Phi_n_i + (1i.*Mx./(2.*k0.*delta_x)).*(Phi_nxt- Phi_prv))); % (Pa)

% 3D pressure field -6dB when x = *2
pp3D = pp2D./sqrt(x(1:end-1));

%% Interpolation (reduce size of matrices for saving)

% interpolation on a grid delat_x,delta_z
[L_phi2D,x_phi,z_phi] = interpolation(Phi,stock_x,stock_z,z,x);
[L_pp2D,x_pp,z_pp] = interpolation(pp2D,stock_x,stock_z,z,x(1:end-1));
[L_pp3D] = interpolation(pp3D,stock_x,stock_z,z,x(1:end-1));

%% Atmospheric absorption

% according to ISO-9613
L_pp2D = SetAtmos(L_pp2D,x_pp,z_pp,hS,freq,hr,T);
L_pp3D = SetAtmos(L_pp3D,x_pp,z_pp,hS,freq,hr,T);

%% Attenuation relative to free field level

[xi,zi] = meshgrid(x_pp,z_pp); % xi & zi are used to avoid loop and calculate in matrix formalism

R1 = sqrt(xi.^2 + (zi-hS).^2); % taking into account acoustic divergence
DeltaL2D = L_pp2D + 20*log10(R1);
DeltaL3D = DeltaL2D - 20*log10(sqrt(x_pp)); % correction 3D by sqrt(x)

%% Atmospheric turbulence

if f_turb_ind == true
    DeltaL2D = fast_turb(DeltaL2D,xi,zi,hS,freq,gamT);
    DeltaL3D = fast_turb(DeltaL3D,xi,zi,hS,freq,gamT);
end

%% Normalisation by maximum

L_phi2D = L_phi2D - max(max(L_phi2D));
L_pp2D = L_pp2D - max(max(L_pp2D));
L_pp3D = L_pp3D - max(max(L_pp3D));

toc

%% Saving results

save([num2str(freq) 'Hz_DeltaL.mat'],'DeltaL3D');
save([num2str(freq) 'Hz_normalized_SPL.mat'],'L_pp3D');

%% Plotting

figure 

subplot(2,1,1) % 2D(x,z) attenuation to free field map (dB)
surf(xi,zi,DeltaL3D)
view(2)
shading interp     
caxis([-30 20])
set(gca,'FontSize',12)
xlabel('Range x (m)')
ylabel('Height z (m)')
cb = colorbar;
cb.Label.String = '\Delta L (dB)';
title('Attenuation rel. to free field')

subplot(2,1,2) % 2D(x,z) SPL (dB)
surf(xi,zi,L_pp3D)
view(2)
shading interp     
caxis([-100 0])
set(gca,'FontSize',12)
xlabel('Range x (m)')
ylabel('Height z (m)')
cb = colorbar;
cb.Label.String = 'SPL rel. to max (dB)';
title('SPL')


figure 

subplot(2,1,1)% attenuation to free field with distance (dB)
hold all
plot(x_pp,DeltaL3D(iz_r,:),'LineWidth',1.5)
xlim([0 dim_x])
ylim([-80 20])
grid on
set(gca,'TickLabelInterpreter','latex')
xlabel('Range x (m)')
ylabel('\Delta L (dB)')
lgd = legend([num2str(freq) ' Hz']);
title(['z mic = ' num2str(z_r) ' m'])

subplot(2,1,2)% SPL with distance (dB)
hold all
plot(x_pp,L_pp3D(iz_r,:),'LineWidth',1.5)
xlim([0 dim_x])
ylim([-100 0])
grid on
set(gca,'TickLabelInterpreter','latex')
xlabel('Range x (m)')
ylabel('SPL rel. to max (dB)')
lgd = legend([num2str(freq) ' Hz']);
title(['z mic = ' num2str(z_r) ' m'])


%% Archive

% figure
% subplot(2,2,1)
% surf(x_pp,z_pp,L_pp3D)
% view(2)
% shading interp     
% caxis([-70 0])
% xlim([0 dim_h])
% ylim([0 dim_v])
% %set(gca,'FontSize',12)
% xlabel('Range x (m)')
% ylabel('Height z (m)')
% cb = colorbar;
% cb.Label.String = 'SPL (dB)';
% %title('3D pressure field')
% 
% subplot(2,2,2)
% surf(xi,zi,DeltaL3D)
% view(2)
% shading interp     
% caxis([-50 20])
% %set(gca,'FontSize',12)
% xlabel('Range x (m)')
% %ylabel('Height z (m)')
% cb = colorbar;
% cb.Label.String = '\Delta L (dB)';
% %title('Attenuation rel. to free field')
% 
% subplot(2,2,[3 4])
% hold all
% plot(x_pp,L_pp3D(icoupe_z,:),'color',[0 0.4470 0.7410],'LineWidth',1.5)
% xlim([0 dim_h])
% ylim([-100 0])
% grid on
% %set(gca,'TickLabelInterpreter','latex')
% xlabel('Range x (m)')
% ylabel('SPL rel. to max (dB)')
% %lgd = legend([num2str(freq) ' Hz']);
% title(['z mic = ' num2str(zcoupe) ' m'])
% 
% yyaxis right
% plot(x_pp,DeltaL3D(icoupe_z,:),'color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
% xlim([0 dim_h])
% ylim([-50 20])
% ylabel('\Delta L (dB)')
% 
% ax = gca;
% ax.YAxis(1).Color = [0 0.4470 0.7410];
% ax.YAxis(2).Color = [0.8500 0.3250 0.0980];
% 
% sgtitle(['f = ' num2str(freq) ' Hz, z source = ' num2str(hS) ' m'])

% figure
% subplot(2,3,1)
% surf(x_phi,z_phi,L_phi2D)
% view(2)
% shading interp     
% caxis([-60 0])
% xlim([0 dim_h])
% ylim([0 dim_v])
% set(gca,'FontSize',12)
% xlabel('Range x (m)')
% ylabel('Height z (m)')
% cb = colorbar;
% %cb.Label.String = 'Velocity potential \phi EWAPE (dB)';
% title('2D velocity potential field \phi')
% 
% subplot(2,3,2)
% surf(x_pp,z_pp,L_pp2D)
% view(2)
% shading interp     
% caxis([-60 0])
% xlim([0 dim_h])
% ylim([0 dim_v])
% set(gca,'FontSize',12)
% xlabel('Range x (m)')
% %ylabel('Height z (m)')
% cb = colorbar;
% %cb.Label.String = 'pressure p EWAPE (dB)';
% title('2D pressure field')
% 
% subplot(2,2,1)
% surf(x_pp,z_pp,L_pp3D)
% view(2)
% shading interp     
% caxis([-60 0])
% xlim([0 dim_h])
% ylim([0 dim_v])
% set(gca,'FontSize',12)
% xlabel('Range x (m)')
% %ylabel('Height z (m)')
% cb = colorbar;
% cb.Label.String = '(dB)';
% title('3D pressure field')
% % 
% % 
% % 
% subplot(2,2,[3 4])
% hold all
% % plot(x_phi,L_phi2D(icoupe_z,:));
% % plot(x_pp,L_pp2D(icoupe_z,:));
% plot(x_pp,L_pp3D(icoupe_z,:));
% xlim([0 dim_h])
% ylim([-60 10])
% grid on
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Range x (m)')
% ylabel('SPL rel. to max (dB)')
% lgd = legend('\phi EWAPE','p 2D EWAPE','p 3D EWAPE');
% lgd.Title.String = ['z = ' [num2str(zcoupe)] ' (m)'];
% 
% sgtitle(['z source = ' num2str(hS) ' (m)'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% subplot(2,2,1)
% surf(xi,zi,DeltaL2D)
% view(2)
% shading interp     
% caxis([-50 20])
% set(gca,'FontSize',12)
% xlabel('Range x (m)')
% ylabel('Height z (m)')
% cb = colorbar;
% cb.Label.String = '\Delta L 2D EWAPE (dB)';
% 
% subplot(2,2,2)
% surf(xi,zi,DeltaL3D)
% view(2)
% shading interp     
% caxis([-50 20])
% set(gca,'FontSize',12)
% xlabel('Range x (m)')
% %ylabel('Height z (m)')
% cb = colorbar;
% cb.Label.String = '\Delta L 3D EWAPE (dB)';
% 
% icoupe_z = round(zcoupe/stock_z);
% 
% subplot(2,2,[3 4])
% hold all
% plot(x_pp,DeltaL2D(icoupe_z,:));
% plot(x_pp,DeltaL3D(icoupe_z,:));
% %ylim([-60 10])
% grid on
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Range x (m)')
% ylabel('SPL rel. to free field (dB)')
% lgd = legend('\Delta L 2D','\Delta L 3D');
% lgd.Title.String = ['z = ' [num2str(zcoupe)] ' (m)'];
% 
% sgtitle(['z source = ' num2str(hS) ' (m)'])
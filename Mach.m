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

function [Mx,epsilon0,epsilon1] = Mach(v_ref,z_ref,z,z0,d,shear_exp,c0,Tlog,temp0,Rg,gam,theta)
% Allows calculation of wind (power law) & temperature (log) profiles,
% as well as mach number Mx and deviation of refractive index epsilon

V_x = (v_ref.*(z./z_ref).^shear_exp)*cos(pi*theta./180);       % power-law wind speed profile through x direction, for each z (m/s)

[~,indice_d] = min(abs(z-d));                   % find the z index that correspond to displacement length

temp1 = temp0+Tlog*(log((z-d)./z0 +1)); % temperature log profil along z (K/m)

cel = sqrt(gam*Rg*temp1);                     % sound celerity that depend on T(z)
cel(1:indice_d) = c0;                                 % sound celerity < displacement height = c0 (no influence of profiles)

epsilon = (c0.^2./cel.^2)-1;                      % deviation of refractive index, see eq 31 from [Ostashev et al 2020]. (scalar along z dimension)
 
Mx = V_x./cel;               % longitudinal mach number through x dimension, for each z (scalar)

epsilon0 = epsilon;        % for A mat
epsilon1 = epsilon;        % for B mat

% figure
% hold on
% plot(V_x,z,'LineWidth',2)
% set(gca,'FontSize',15)
% grid on
% 
% figure(12)
% hold on
% plot(temp1,z,'LineWidth',2)
% set(gca,'FontSize',15)
% grid on

end
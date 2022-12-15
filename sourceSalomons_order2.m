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

function U0 = sourceSalomons_order2(z,hS,nb_z,k0,Rp)
% source: PhD thesis Didier Dragna (2011) and Salomon's book (2001)
% difference of -1/(4*pi) between the 2 probably coming from Green's
% function definition

% Rp : reflexion coefficient
% z : vertical axisl
% hS : height of source
% k0 : adiabatic wavenumber

U01 = zeros(nb_z,1); % starter
U02 = zeros(nb_z,1); % starter
U0 = zeros(nb_z,1); % starter

% direct wave
k02z2 = k0^2*(z-hS).^2;
U01 = sqrt(k0*1j)*(1.3717-0.3701*k02z2).*exp(-k02z2/3);
% reflected wave
k02z2 = k0^2*(z+hS).^2;
U02 = sqrt(k0*1j)*(1.3717-0.3701*k02z2).*exp(-k02z2/3);
% starter
U0 = U01 + Rp*U02;

% figure
% hold on
% plot(abs(U0),z,'LineWidth',2)
% set(gca,'FontSize',15)
% grid on
end
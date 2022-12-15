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

function alpha=atmos(f,hr,T,pa)
% alpha=atmos(f,hr,T,pa)
% in dB/m
% atmospheric absorption coefficient calculation
% according to ISO 9613-1 (1993)
%
% f: Hz (default 1000Hz)
% hr : % (default 50%)
% T : deg C (defaut 20 deg C)
% pa : pressure atmos in kPa (ref=101.325)
%

if nargin==0
 f=1000;
 hr=50;
 T=20;
 pa=101.325;
end
if nargin==1
 hr=50;
 T=20;
 pa=101.325;
end
if nargin==2
 T=20;
 pa=101.325;
end
if nargin==3
 pa=101.325;
end

% reference values
pr=101.325;
To=293.15;
To1=273.16;

% T in deg K
T=T+273.15;

% h : fraction molaire de vapeur d'eau
C=-6.8346.*(To1./T).^1.261+4.6151;
h=hr.*10.^C.*pr./pa;

% Eq (3) (4) and (5) p3 of ISO 9613-1 standard
% Frequence de relaxation O2    
fro=pa./pr.*(24+4.04.*1e4.*h.*(0.02+h)./(0.391+h));
% Frequence de relaxation N
frn=pa./pr.*(T/To).^(-0.5).*(9+280.*h*exp(-4.170*((T/To).^(-1/3)-1)));

% Attenuation in dB/m
alpha=8.686*f.^2.*((1.84*1e-11*(pr./pa).*(T./To).^(0.5))+(T./To).^(-5/2)*(0.01275*(exp(-2239.1/T)).*(fro+(f.^2/fro)).^(-1) + 0.1068*(exp(-3352.0./T)).*(frn+(f.^2/frn)).^(-1)));


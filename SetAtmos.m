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

function ppdBiAtmos = SetAtmos(ppdBi,x,z,haut,freq,hr,T)
% Apply atmospheric absorption to pressure field
hr = hr*100;        % hr in percent

% calculation of spatial grid from x z vectors
[Xi,Zi] = meshgrid(x,z);

xs = 0;
zs = haut;
Ri = sqrt((Xi-xs).*(Xi-xs)+(Zi-zs).*(Zi-zs));
alpha = atmos(freq,hr,T);
ppdBiAtmos = ppdBi-alpha*Ri;
clear Ri



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

function ppdBiTurb = fast_turb(ppdBi,Xi,Zi,haut,freq,gamT)

% see eq.11 of [Van Renterghem et al 2022] Applied Acoustics
%
%  taken from Harmonoise project

xs = 0; % position x of source (m)
zs = haut; % position z of source (m)
d = sqrt((Xi-xs).*(Xi-xs)+(Zi-zs).*(Zi-zs)); % source receiver distance (m)

Turb = 25 + 10*log10(gamT) + 3*log10(freq/1000) + 10*log10(d./100); 

ppdBiTurb = 10*log10(10.^(ppdBi./10)+10.^(Turb./10)); % pressure field with correction
end
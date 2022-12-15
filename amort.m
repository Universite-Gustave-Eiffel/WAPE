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

function U11=amort(U0,haut_a,coeff_a)
% Damping at the top of the domain to avoid unrealistic reflexions

dim=length(U0); % length of z axis

za = haut_a*dim; % height where the damping start
idex = [ceil(za):dim-1]'; % index of damping
 
U0(idex)=U0(idex).*exp(  -(     (idex-za)./(coeff_a*(dim-idex))         ).^2); % damping
U0(dim)=0;
U11=U0;

end
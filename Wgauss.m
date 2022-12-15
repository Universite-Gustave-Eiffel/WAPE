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

function W = Wgauss(kappa, sigmah, lc )
% Rugosity spectrum calculation needed to calculate the effect of ground roughness
% on effective admittance (see impedance.m)
%
% kappa = k0*sin(angle_incidence)   (modified wave number)
% sigmah : standard deviation of roughness height
% lc : correlation length

%%

W = ((sigmah.^2)*(lc/(2*sqrt(pi))))*exp(-((kappa.*lc).^2)/4) ;

end


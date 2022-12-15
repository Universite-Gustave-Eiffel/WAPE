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

function [ppdBi,xi,zi] = interpolation(field,stock_x,stock_z,z,x)
% Interpolation of the results on a regular grid, the grid is larger than
% lambda/10 or 20, which allow to reduce matrices size. The interpolation is
% sufficient to have precise results (each 0.5m)
%

% the actual 2D vectors discretized in lambda/xx (usually 10 or 20)
X = x;
Z = z;

% 2D vector of interpolation
xi = 0:stock_x:X(end);
zi = 0:stock_z:Z(end);

% linear interpolation, calculation of dB levels
ppdBi = 20*log10(abs(interp2(X,Z,field,xi,zi','linear')));

end    
    
    


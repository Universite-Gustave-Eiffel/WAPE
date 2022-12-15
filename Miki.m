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

function [beta]=Miki(freq,k0,sigma,h)
% Miki impedance model see [Miki 1990] and [Kirby 2014]

    if sigma >= 100000 
        beta = zeros(size(freq));         % perfecty flat and rigid
    else
        Zc = 1+5.50*(freq./sigma).^(-0.632)+1i*8.43*(freq./sigma).^(-0.632); % characteristic impedance // Z0
        kc = k0.*(1+7.81*(freq./sigma).^(-0.618)+1i*11.41*(freq./sigma).^(-0.618)); % characteristic wavenumber
        
        Zs = thickness(Zc,h,kc);  % surface impedance

        beta = 1./Zs;
    end
    
end


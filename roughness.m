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

function [beta_rough] = roughness(lc,sigmah,k0,incidence)
% Calculates admittance correction of a rough ground
% see Olivier Faure 2014 thesis, Annex A
% see [Kayser et al 2019]

kappa = k0*sin(incidence);
s1 = 1 ; 
s2 = -1 ; % numerical parameter

 % test if roughness
    if lc == 0
        beta_rough = 0;
    else
        d_int_alpha=sqrt(k0) / 100 ;
        u_alpha = 0:d_int_alpha:sqrt(k0) ;

        integrande_alpha1 = real((1./(k0*sqrt(-u_alpha.^2+2*k0))).*((k0^2+s1*kappa*(k0-u_alpha.^2)).^2).*Wgauss(kappa+s1*(k0-u_alpha.^2), sigmah, lc)) ;
        integrande_alpha2 = real((1./(k0*sqrt(-u_alpha.^2+2*k0))).*((k0^2+s2*kappa*(k0-u_alpha.^2)).^2).*Wgauss(kappa+s2*(k0-u_alpha.^2), sigmah, lc)) ;

        alpha1 = trapz(u_alpha, integrande_alpha1) ;
        alpha2 = trapz(u_alpha, integrande_alpha2) ;
        alpha = alpha1 + alpha2 ;
    
        d_int_beta=(6/lc)/100 ;
        u_beta = 0:d_int_beta:(6/lc) ;

        integrande_beta1 = real((1./(k0*sqrt(k0^2+u_beta.^2))).*((k0^2+s1*kappa*sqrt(k0^2+u_beta.^2)).^2).*Wgauss(kappa+s1*sqrt(k0.^2+u_beta.^2), sigmah, lc)) ;
        integrande_beta2 = real((1./(k0*sqrt(k0^2+u_beta.^2))).*((k0^2+s2*kappa*sqrt(k0^2+u_beta.^2)).^2).*Wgauss(kappa+s2*sqrt(k0.^2+u_beta.^2), sigmah, lc)) ;
   
        beta1 = -trapz(u_beta, integrande_beta1) ;
        beta2 = -trapz(u_beta, integrande_beta2) ;
        beta = beta1 + beta2 ;
    
        beta_rough = alpha + 1i*beta ;
    end

end
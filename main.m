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

function [Phi] = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,epsilon0,epsilon1,Beta,haut_a,coeff_a)

Phi = zeros(nb_z,nb_x); % initialization of variable that stores the velocity potential

%% Amat

a_vect = (2 - 1i*k0*delta_x .*(2* gamma_x.^2 - tau_bar1)   )./(8*k0.^2* gamma_x.^2*delta_z.^2 );

b_vect = 1 +   epsilon1/4 ...
    -  1i*k0.*delta_x/2.*(gamma_x.^2 .* epsilon1./2  - (1+ epsilon1/4).* tau_bar1 )...
    -  (2 - 1i*k0*delta_x .*(2* gamma_x.^2 - tau_bar1)   )./(4*k0.^2* gamma_x.^2*delta_z.^2 ) ;

c_vect = a_vect;

% boundary conditions
b_vect(1) = b_vect(1) + a_vect(1)* (2*1i*delta_z*k0*Beta);
c_vect(1) = 2*c_vect(1);

A_b = sparse(1:nb_z,1:nb_z,1,nb_z,nb_z);  % diagonal term
A_a = sparse(2:nb_z,1:nb_z-1,1,nb_z,nb_z); % bottom
A_c = sparse(1:nb_z-1,2:nb_z,1,nb_z,nb_z); % top

Amat = b_vect.*A_b +  a_vect.*A_a+c_vect.*A_c;

%% Bmat
d_vect = (2 + 1i*k0*delta_x .*(2* gamma_x.^2 - tau_bar1)   )./(8*k0.^2* gamma_x.^2*delta_z.^2 );

e_vect = 1 +   epsilon0/4 ...
    +  1i*k0.*delta_x/2.*(gamma_x.^2 .* epsilon0./2  - (1+ epsilon0/4).* tau_bar1 )...
    - (2 + 1i*k0*delta_x .*(2* gamma_x.^2 - tau_bar1)   )./(4*k0.^2* gamma_x.^2*delta_z.^2 ) ;

f_vect = d_vect;

% boundary conditions
e_vect(1) = e_vect(1) + d_vect(1)* (2*1i*delta_z*k0*Beta);
f_vect(1) = 2*f_vect(1);

B_e = sparse(1:nb_z,1:nb_z,1,nb_z,nb_z);  % diagonal term
B_d = sparse(2:nb_z,1:nb_z-1,1,nb_z,nb_z); % bottom
B_f = sparse(1:nb_z-1,2:nb_z,1,nb_z,nb_z); % top

Bmat = e_vect.*B_e +  d_vect.*B_d+f_vect.*B_f;

%% Main calculation
clear a_vect b_vect c_vect d_vect e_vect f_vect

for ix = 1:nb_x
    
    % first matrix resolution
    U0=Bmat*U0;
    
    % second matrix resolution
    U0 = Amat\U0;
    
    % absorbing layer
    U0 = amort(U0,haut_a,coeff_a);

    % save velocity field after 1 iteration
    Phi(:,ix) = U0; % It is "Phi_hat" of the EWAPE given in eq 79 from [Ostashev et al 2020]
       
end
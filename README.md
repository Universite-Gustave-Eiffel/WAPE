# Wide Angle Parabolic Equation model with arbitrary Mach number

<img src="/image.png" width="700" />

The WAPE model proposed here is an implementation of the work of [Ostashev et al 2020]. The WAPE is developed using a Padé(1,1) series expansion. A Crank-Nicholson algorithm is used to reduce the equation into a matrix system. Then, the 2D acoustic pressure field is calculated using a second-order centered finite difference scheme on velocity potential.

This code is intended to be shared with the research community.

## Installation
This WAPE model is designed to run with [Matlab](https://www.mathworks.com/products/matlab.html), without any additonal toolbox.

## Run a calculation
To run a calculation you need to open "initialisation_WAPE.m", set parameters values, and run the script. It is indicated which parameter can be modified and which can not.

The input parameters are:

```Matlab
%% Source parameters
freq = 250;  % source frequency (Hz)
hS = 35;     % source height (m)

%% Medium properties
T = 10;      % atmospheric temperature at the surface (°C)
Tlog = 0.2;  % temperature coefficient for the vertical gradient

shear_exp = 0.15;    % wind shear exponents (scalar) for power law wind profil
v_ref = 5;          % wind speed (m/s) measured at z_ref height, for power-law wind profile
z_ref = 80;          % reference height for wind speed v_ref
theta = 0;         % propagation angle with respect to the source (0° : downwind, 180° : upwind)

f_turb_ind = true;   % logical 'true' or 'false' to account for turbulence or not
gamT = 0;            % turbulence strength

%% 2D Spatial domain (x,z)
dim_x = 3000; % horizontal dimension x of the domain (m)
dim_z = 300;  % vertical dimension z of the domain (m)

haut_a = 0.8; % where the absorbing layer starts at the top of the domain, according to z axis (0.8 ==> start at 80% of z axis)
coeff_a = 10; % damping rate

discrx = 10;  % discretization x = lambda/discrx
discrz = 10;  % discretization z = lambda/discrz

z_r = 1.5;    % receiver height for plotting (m)

%% Ground properties
hv = 0;       % vegetation height (m), it affects the shape of atmospheric profils
z0 = 0.13*hv + 0.00001;  % atmospheric roughness length (m), can't be null
d = 0.66*hv;   % displacement height of flux profiles (m), it's directly linked to vegetation height

lc = 0 ;       % correlation length (m) (ground roughness parameter), if 0 : no ground rugositiy, [0.05-1]
sigmah = 0;    % standard deviation of roughness height (m) (ground roughness parameter), [0.01-0.05]
sigma = 10000000;  % airflow resistivity of the ground (kN.s.m-4)
h = 0;         % thickness of ground surface layer (m), if h = 0 no layer
```

## Post-processing
When the simulation has completed there will be 'DeltaL.mat' and  'normalized_SPL.mat' files which correspond respectively to attenuation to free field (dB), and sound pressure level field (dB) which is normalized by the maximum amplitude. You can post-process these signals to your liking. Here is an exemple of output:

<img src="/SPL.eps" width="500" /> 

## Warning
Instabilities can occur with PE simulations if:
- the spatial steps are too large,
- the absorbent layer at the top of the domain is not thick enough.

Calculation time can be really long and take a lot of memory if: 
- the spatial steps are too small,
- the domain is too big,
- the frequency is too high.

## License
This code is released under the EUPL license.

## Credits
The development of this code took part into the french [PIBE project](https://www.anr-pibe.com/) (contract ANR-18-CE04-0011).

## Citing this work
If this WAPE script contributes to an academic publication, please cite it as:

      @misc{kayser2023wape,
        title = {WAPE model},
        author = {Bill Kayser},
        note = {https://github.com/bkayser13/WAPE/},
        year = {2023}
      }

## Some background references

This code largely implements algorithms that have already been published. A non-exhaustive list is proposed below.

Tappert, F. D. (1977). The parabolic approximation method. Wave propagation and underwater acoustics, 224-287.

Gilbert, K. E., & White, M. J. (1989). Application of the parabolic equation to sound propagation in a refracting atmosphere. The Journal of the Acoustical Society of America, 85(2), 630-637.

West, M., Gilbert, K., & Sack, R. A. (1992). A tutorial on the parabolic equation (PE) model used for long range sound propagation in the atmosphere. Applied Acoustics, 37(1), 31-49.

Collins, M. D. (1993). A split‐step Padé solution for the parabolic equation method. The Journal of the Acoustical Society of America, 93(4), 1736-1742.

Ostashev, V. E., Juvé, D., & Blanc-Benon, P. (1997). Derivation of a wide-angle parabolic equation for sound waves in inhomogeneous moving media. Acta Acustica united with Acustica, 83(3), 455-460.

Dallois, L., Blanc-Benon, P., & Juvé, D. (2001). A wide-angle parabolic equation for acoustic waves in inhomogeneous moving media: Applications to atmospheric sound propagation. Journal of Computational Acoustics, 9(02), 477-494.

Lihoreau, B., Gauvreau, B., Bérengier, M., Blanc-Benon, P., & Calmet, I. (2006). Outdoor sound propagation modeling in realistic environments: Application of coupled parabolic and atmospheric models. The Journal of the Acoustical Society of America, 120(1), 110-119.

Cheinet, S. (2012). A numerical approach to sound levels in near-surface refractive shadows. The Journal of the Acoustical Society of America, 131(3), 1946-1958.

Kayser, B., Gauvreau, B., & Ecotière, D. (2019). Sensitivity analysis of a parabolic equation model to ground impedance and surface roughness for wind turbine noise. The Journal of the Acoustical Society of America, 146(5), 3222-3231.

Ostashev, V. E., Muhlestein, M. B., & Wilson, D. K. (2019). Extra-wide-angle parabolic equations in motionless and moving media. The Journal of the Acoustical Society of America, 145(2), 1031-1047.

Ostashev, V. E., Wilson, D. K., & Muhlestein, M. B. (2020). Wave and extra-wide-angle parabolic equations for sound propagation in a moving atmosphere. The Journal of the Acoustical Society of America, 147(6), 3969-3984.

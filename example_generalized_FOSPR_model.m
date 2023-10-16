%% reset workspace and add necessary folders to path
close all
clear variables
addpath(strcat(pwd,'\permittivities'));
addpath(strcat(pwd,'\generalized_FOSPR_model'));

%% angular probability density function of rays from Monte Carlo simulation
angular_pdf = generate_angular_pdf( ...
    'RTMCS_complete', ...
    'bin_number', 64, ...
    'plot_bool', true);

%% inputs from optimization of the generalized FO-SPR model
dAu = 48.936004132662674;% gold thickness [nm]

% variables of Drude-Lorentz equation for relative permittivity of gold
epsilon_infinity = 1.147814238704828e-08;
omega_p = 2.029622338315185e+03;
gamma_D = 28.587727162907544;
DELTA_epsilon = 4.539019406602927;
OMEGA_L = 8.131149415666220e+02;
GAMMA_L = 1.031231302322861e+02;

%% FO-SPR probe parameters
lambda = (550:1:900)';% [nm]
d_layers = dAu;% [nm]
d_layers_reference = d_layers;
L = 2 * 10;% [mm] times two because of back reflection
D = 400;% [µm]

%% convert to SI-units
lambda_m = lambda .* 1e-9;% nm to m
d_layers_m = d_layers .* 1e-9;% nm to m
d_layers_reference_m = d_layers_reference .* 1e-9;% nm to m
L_m = L * 1e-3; % mm to m
D_m = D * 1e-6; % µm to m


%% relative permittivities of layers
% relative permittivity of gold
permittivity_Au = @(lambda_m_var)permittivity_Au_Drude_Lorentz(...
    lambda_m_var, epsilon_infinity, omega_p, gamma_D, DELTA_epsilon, ...
    OMEGA_L, GAMMA_L);

% relative permitivities of [glass gold sample(water)]
permittivities_layers = [ ...
    permittivity_SiO2(lambda_m), ...
    permittivity_Au(lambda_m), ...
    permittivity_H2O(lambda_m) ...
    ];
% relative permitivities of [glass gold air]
permittivities_layers_reference = [ ...
    permittivity_SiO2(lambda_m), ...
    permittivity_Au(lambda_m), ...
    permittivity_air(lambda_m) ...
    ];

%% reflectance of FO-SPR
R = generalized_FOSPR_model( lambda_m, ...
    permittivities_layers, ...
    permittivities_layers_reference, ...
    d_layers_m, ...
    d_layers_reference_m, ...
    L_m, D_m, ...
    angular_pdf);

figure
plot(lambda, R)
xlabel('wavelength (nm)');
ylabel('R_{rel}');
function [epsilon_air] = permittivity_air(lambda)
% Dispersion relation of air.
%
% epsilon_air : relative permittivity of air
%
% lambda : wavelength of light; [m]
%
% Rayleigh-scattering calculations for the terrestrial atmosphere
% Bucholtz, Anthony
%
%% Cauchy relation (230 till 4000 nm at 15 °C)
lambda = lambda .* 1e6;% m to µm

A0 = 0.05791817;
A1 = 238.0185;
A2 = 0.00167909;
A3 = 57.362;

epsilon_air = (1 + A0./(A1 - lambda.^-2) + A2./(A3 - lambda.^-2)).^0.5;
end
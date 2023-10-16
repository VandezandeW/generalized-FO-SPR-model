function [ epsilon_Au ] = permittivity_Au_Drude_Lorentz( lambda, ...
    epsilon_infinity, omega_p, gamma_D, DELTA_epsilon, OMEGA_L, GAMMA_L)
% Calculates the relative permittivity of gold with the Drude-Lorentz
% equation.
%
% epsilon_Au : relative permittivity of gold
%
% lambda : light wavelengths; [m]
% epsilon_infinity : background relative permittivity
% omega_p : angular plasma-frequency; [THz/(2*pi)]
% gamma_D : angular dampening-frequency; [THz/(2*pi)]
% DELTA_epsilon : Lorentz weight factor
% OMEGA_L : Lorentz angular oscillator-strength; [THz/(2*pi)]
% GAMMA_L : Lorentz angular dampening-frequency; [THz/(2*pi)]

% Based on:
% Improved analytical fit of gold dispersion: Application to the modeling
% of extinction spectra with a finite-difference time-domain method
% Alexandre Vial,* Anne-Sophie Grimault, Demetrio Macías, Dominique
% Barchiesi, and Marc Lamy de la Chapelle

%%
c = physconst('LightSpeed');%[m/s]
omega = 2 .* pi .* c ./ lambda;%[rad * m/s / m] = [rad/s] (eq. 13)

omega_p = omega_p * 2 * pi * 1e12;% [Hz]
gamma_D = gamma_D * 2 * pi * 1e12;% [Hz]
OMEGA_L = OMEGA_L * 2 * pi * 1e12;% [Hz]
GAMMA_L = GAMMA_L * 2 * pi * 1e12;% [Hz]

% Drude-Lorentz equation(eq. 11)
epsilon_Au = epsilon_infinity ...
    - (omega_p.^2 ./ (omega .* (omega + 1i .* gamma_D))) ...
    - (DELTA_epsilon .* OMEGA_L.^2 ./ (omega.^2 - OMEGA_L.^2 + 1i .* ...
    GAMMA_L .* omega));
end
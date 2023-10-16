function [ epsilon_H2O ] = permittivity_H2O( lambda )
% Dispersion relation of water.
%
% epsilon_H2O : relative permittivity of water
%
% lambda : wavelength of light; [m]
%
% Linear refractive index and absorption
% measurements of nonlinear optical
% liquids in the visible and near-infrared
% spectral region
% S. Kedenburg, M. Vieweg, T. Gissibl, and H. Giessen
%
%% Sellmeier formula (500 - 1600 nm @ 20°C)
% warning('only valid for T at 20°C')
lambda = lambda .* 1e6;% m to µm

A1 = 0.75831;
B1 = 0.01007;
A2 = 0.08495;
B2 = 8.91377;

epsilon_H2O = 1 + ...
    (A1 .* lambda.^2) ./ (lambda.^2 - B1) + ...
    (A2 .* lambda.^2) ./ (lambda.^2 - B2);
end


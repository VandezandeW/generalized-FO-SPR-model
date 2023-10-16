function [ epsilon_EtOH ] = permittivity_EtOH( lambda )
% Dispersion relation of ethanol.
%
% epsilon_EtOH : relative permittivity of ethanol
%
% lambda : wavelength of light; [m]
%
% Linear refractive index and absorption
% measurements of nonlinear optical
% liquids in the visible and near-infrared
% spectral region
% S. Kedenburg, M. Vieweg, T. Gissibl, and H. Giessen
%
%% Snellmeier formula (500 - 1600 nm @ 20°C)
% warning('only valid for T at 20°C')
lambda = lambda .* 1e6;% m to µm

A1 = 0.83189;
B1 = 0.00930;
A2 = -0.15582;
B2 = -49.452;

epsilon_EtOH = 1 + ...
    (A1 .* lambda.^2) ./ (lambda.^2 - B1) + ...
    (A2 .* lambda.^2) ./ (lambda.^2 - B2);
end


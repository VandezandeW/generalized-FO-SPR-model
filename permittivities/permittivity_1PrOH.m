function [ epsilon_1PrOH ] = permittivity_1PrOH( lambda )
% Dispersion relation of 1-propanol.
%
% epsilon_1PrOH : relative permittivity of 1-propanol
%
% lambda : wavelength of light; [m]
%
% Refractive, dispersive and thermo-optic properties of twelve organic
% solvents in the visible and near-infrared
% Konstantinos Moutzouris, Myrtia Papamichael, Sokratis C. Betsis,
% Ilias Stavrakas, George Hloupis, Dimos Triantis
%
%% extended-Cauchy formula (450 - 1551 nm @ 27°C)
% warning('only valid for T at 27°C')
lambda = lambda .* 1e6;% m to µm

A0 = 1.893400242;
A1 = -0.003349425;
A2 = 0.004418653;
A3 = 0.00108023;
A4 = -0.000067337;

epsilon_1PrOH = (A0 + ...
    A1 .* lambda.^2 + ...
    A2 ./ lambda.^2 + ...
    A3 ./ lambda.^4 + ...
    A4 ./ lambda.^6);
end


function [ epsilon_1BuOH ] = permittivity_1BuOH( lambda )
% Dispersion relation of 1-butanol.
%
% epsilon_1BuOH : relative permittivity of 1-butanol
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

A0 = 1.917816501;
A1 = -0.00115077;
A2 = 0.01373734;
A3 = -0.00194084;
A4 = 0.000254077;

epsilon_1BuOH = A0 + ...
    A1 .* lambda.^2 + ...
    A2 ./ lambda.^2 + ...
    A3 ./ lambda.^4 + ...
    A4 ./ lambda.^6;
end


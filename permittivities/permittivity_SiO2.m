function [ epsilon_SiO2 ] = permittivity_SiO2( lambda )
% Dispersion relation of water.
%
% epsilon_SiO2 : relative permittivity of glass
%
% lambda : wavelength of light; [m]
%
%% Cauchy relation optimized with thorlabs data (eq. 11)
%http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6845
% RI:
%   436 nm: 1.467287
%   589.3 nm: 1.458965
%   1020 nm: 1.450703

lambda = lambda .* 1e6; % m to µm

A0 = 2.103317754190267;
A1 = -0.007799439766211;
A2 = 0.009713150628039;

epsilon_SiO2 = A0 + A1 .* lambda.^2 + A2 .* lambda.^-2;

end


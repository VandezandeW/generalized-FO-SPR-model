function [R] = generalized_FOSPR_model( ...
    lambda, ...
    permittivities_layers, ...
    permittivities_layers_reference, ...
    d_layers, d_layers_reference, L, D, ...
    angular_pdf)
% This model assumes that the FO-SPR sensing surface is made out of
% isotropic layers where the first and last 'layers' are semi-infinite
% volumes. Each layer is indexed, starting from 1.
%
% R: relative reflectance
%
% lambda [m]: light wavelengths;
% permittivities_layers: permitivities of the different layers
% permittivities_layers_reference: permitivities of the different layers
%   for referencing. Usually the sample is replaced with air.
% d_layers [m]: layer thicknesses of second to second last layer;
% d_layers_reference [m]: thicknesses of the referencing layers of second to
%   second last layer. Usually the same as d_layers;
% L [m]: length of optical fiber where SPR takes place, backreflection
%   included, e.g., 6 mm probelength with backreflection: L = 12 mm;
% D [m]: diameter of optical fiber where SPR takes place;
% angular_pdf: gridded data interpolant object of the angular probability
%   density function

%% reflectance of FO-SPR sensor (eq. 32)
P = angular_power_integration(d_layers, ...
    permittivities_layers);% reflected light power of sample
P_reference = angular_power_integration(d_layers_reference, ...
    permittivities_layers_reference);% reflected light power of reference
R = P ./ P_reference;% relative reflectance

    function [ P ] = angular_power_integration(d_layers, ...
            permittivities_layers)
        N = size(permittivities_layers, 2);
        %% integration of the reflection power over all angles
        alpha = angular_pdf.GridVectors{1};% reflection angle [rad]
        gamma = angular_pdf.GridVectors{2};% skew angle [rad]

        P_trapz = reflectance_power(...
            lambda, alpha, gamma, permittivities_layers, L, D, N, ...
            d_layers, angular_pdf);
        P = trapz(alpha, P_trapz, 2);
        if size(P_trapz, 3) ~= 1
            P = trapz(gamma, P, 3);% reflected light power
        end
    end
end

function [P] = reflectance_power(lambda, ...
    alpha, gamma, epsilon, L, D, N, d, angular_pdf)
% P: power of the reflected light.
%
% alpha [rad]: reflection angle in the first, semi-infinite layer;
% gamma [rad]: skew angle in the first, semi-infinite layer;
% lambda [m]:  light wavelength;
% epsilon : permitivities of the different layers
% L [m]: total length of optical fiber where SPR takes place,
%   back reflection included, e.g., 6 mm probelength with back reflection:
%   L = 12 mm;
% D [m]: diameter of optical fiber where SPR takes place;
% N : number of layers;
% d [m]: layer thicknesses of from second to second last layer;
% angular_pdf: gridded data interpolant object of the angular probability
%   density function.

%% angular power distribution
[alpha_grid, gamma_grid] = ndgrid(alpha, gamma);
A_alpha_gamma = angular_pdf(alpha_grid, gamma_grid);% (eq. 16)

% remove unsubstantial probability densities to shorten calculation time
A_alpha_gamma(A_alpha_gamma < mean(A_alpha_gamma, 'all') ./ 1e4) = 0;

if ~any(A_alpha_gamma)
    % If there are no rays, there is no power.
    P = zeros(size(lambda, 1), size(A_alpha_gamma, 1), ...
        size(A_alpha_gamma, 2));
else
    %% Fresnel coefficients of reflection (eq. 5-8)
    [r_s, r_p] = reflection_coefficients_of_multilayer(lambda, alpha, ...
        epsilon, N, d);
    
    %% meridional angle theta (eq. 15)
    theta = real(asin(cos(alpha') ./ cos(gamma)));
    
    %% numbder of reflections (eq. 27)
    % (Note that L is the total length of the FO-SPR probe. Hence, the
    % length needs to be doubled before it is given to this function for a
    % back-reflecting FO-SPR sensor.)
    N_refl = (L  .* tan(theta)) ./ (D .* cos(gamma));
    
    % Fraction of rays that have one reflection more compared to the other 
    % rays dependening on the position of these rays at the fiber entrance.
    % This assumes a homogeneous distribution over the fiber entrance.
    N_refl_plus_1_fraction = N_refl - floor(N_refl);
    N_refl = floor(N_refl);
    
    % If there are no rays, there are no reflections.
    N_refl(A_alpha_gamma == 0) = 0;
    
    %% reorder dimensions to [lambda alpha gamma]
    N_refl_plus_1_fraction = permute(N_refl_plus_1_fraction, [3 1 2]);
    A_alpha_gamma = permute(A_alpha_gamma, [3 1 2]);
    
    %% average reflectance of rays with the same reflection angles
    [R_N, R_N_plus_1] = reflectance_per_alpha_gamma(alpha, gamma, ...
        r_s, r_p, N_refl);
    R = (1 - N_refl_plus_1_fraction) .* R_N + ...
        N_refl_plus_1_fraction .* R_N_plus_1;
    
    %% power of the reflected light
    P = R .* A_alpha_gamma;
end
end
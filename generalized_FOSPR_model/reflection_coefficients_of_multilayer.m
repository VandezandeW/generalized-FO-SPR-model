function [ r_s, r_p ] = reflection_coefficients_of_multilayer( ...
    lambda, alpha_1, epsilon, N, d)
% r_s: reflection intensity coefficient of s-polarized light;
% r_p: reflection intensity coefficient of p-polarized light.
%
% alpha_1 [rad]: reflection angle in the first, semi-infinite layer;
% lambda [m]:  light wavelength;
% epsilon: permitivities of the different layers;
% N: number of layers;
% d [m]: layer thicknesses of from second to second to last layer.

%% Characteristic matrices M from the second to the second last layer
M_s_11_k = nan(size(lambda, 1), size(alpha_1, 2));
M_s_12_k = M_s_11_k;
M_s_21_k = M_s_11_k;
M_s_22_k = M_s_11_k;
M_p_11_k = M_s_11_k;
M_p_12_k = M_s_11_k;
M_p_21_k = M_s_11_k;
M_p_22_k = M_s_11_k;
for k = 2:N - 1
    epsilon_k = epsilon(:, k);
    [M_s_11_k(:, :, k), M_s_12_k(:, :, k), ...
        M_s_21_k(:, :, k), M_s_22_k(:, :, k), ...
        M_p_11_k(:, :, k), M_p_12_k(:, :, k), ...
        M_p_21_k(:, :, k), M_p_22_k(:, :, k)] = ...
        charactaristric_matrix(epsilon_k, alpha_1, d(k - 1), ...
        epsilon(:, 1), lambda);
end

%% Inverse dynamical matrices D_inv of the first semi-infinite layer
%(derived from eq. 5 and 6)
D_1_inv_s_11 = 1/2;
D_1_inv_s_12 = sec(alpha_1) ./ (2 .* epsilon(:, 1).^0.5);
D_1_inv_s_21 = D_1_inv_s_11;
D_1_inv_s_22 = -D_1_inv_s_12;

D_1_inv_p_11 = sec(alpha_1)./2;
D_1_inv_p_12 = 1./(2 .* epsilon(:,1).^0.5);
D_1_inv_p_21 = D_1_inv_p_11;
D_1_inv_p_22 = -D_1_inv_p_12;

%% Dynamical matrices D of the last semi-infinite layer (eq. 5 and 6)
alpha_N = asin((epsilon(:, 1).^0.5) ./ ...
    (epsilon(:, N).^0.5) .* sin(alpha_1));

D_N_s_11 = 1;
D_N_s_12 = D_N_s_11;
D_N_s_21 = epsilon(:, N).^0.5 .* cos(alpha_N);
D_N_s_22 = -D_N_s_21;

D_N_p_11 = cos(alpha_N);
D_N_p_12 = D_N_p_11;
D_N_p_21 = epsilon(:, N).^0.5;
D_N_p_22 = -D_N_p_21;

%% Characteristic matrices M of the combined structure (eq. 7)
M_s_11 = D_1_inv_s_11;
M_s_12 = D_1_inv_s_12;
M_s_21 = D_1_inv_s_21;
M_s_22 = D_1_inv_s_22;
for j = 2:N-1
    [M_s_11, M_s_12, M_s_21, M_s_22] = matrix_product_2by2(M_s_11, ...
        M_s_12, M_s_21, M_s_22, M_s_11_k(:,:,j), M_s_12_k(:,:,j), ...
        M_s_21_k(:,:,j), M_s_22_k(:,:,j));
end

M_p_11 = D_1_inv_p_11;
M_p_12 = D_1_inv_p_12;
M_p_21 = D_1_inv_p_21;
M_p_22 = D_1_inv_p_22;
for j = 2:N-1
    [M_p_11, M_p_12, M_p_21, M_p_22] = matrix_product_2by2(M_p_11, ...
        M_p_12, M_p_21, M_p_22, M_p_11_k(:,:,j), M_p_12_k(:,:,j), ...
        M_p_21_k(:,:,j), M_p_22_k(:,:,j));
end

[M_s_11, ~, M_s_21, ~] = matrix_product_2by2(M_s_11, M_s_12, M_s_21, ...
    M_s_22, D_N_s_11, D_N_s_12, D_N_s_21, D_N_s_22);
[M_p_11, ~, M_p_21, ~] = matrix_product_2by2(M_p_11, M_p_12, M_p_21, ...
    M_p_22, D_N_p_11, D_N_p_12, D_N_p_21, D_N_p_22);

%% Reflection intensity coefficients R (eq. 8)
r_s = M_s_21 ./ M_s_11;
r_p = M_p_21 ./ M_p_11;

end

function [M_k_s_11, M_k_s_12, M_k_s_21, M_k_s_22, M_k_p_11, M_k_p_12, ...
    M_k_p_21, M_k_p_22] = charactaristric_matrix(epsilon_k, alpha_1, ...
    d_k, epsilon_1, lambda)
% M_k_s_ij : the characteristic matrix element in row i and column j for
% s-polarized light in layer k
% M_k_p_ij : the characteristic matrix element in row i and column j for
% p-polarized light in layer k
%
% epsilon_k : permittivity of layer k
% theta_1 : reflection angle in the first, semi-infinite layer; [rad]
% d_k : layer thicknesses of layer k; [m]
% epsilon : permitivities of the different layers
% lambda :  light wavelength; [m]

%% the characteristic matrix of the kth layer (eq. 1)
q_k_s = (epsilon_k - epsilon_1 .* sin(alpha_1).^2).^0.5;% (eq. 3)
q_k_p = (epsilon_k ./ q_k_s);% (eq. 4)
beta_k = (2 * pi * d_k) ./ lambda .* q_k_s;% (eq. 2)
M_k_12_i = -1i .* sin(beta_k);

M_k_s_11 = cos(beta_k);
M_k_s_22 = M_k_s_11;
M_k_s_12 = M_k_12_i ./ q_k_s;
M_k_s_21 = M_k_12_i .* q_k_s;

M_k_p_11 = M_k_s_11;
M_k_p_22 = M_k_p_11;
M_k_p_12 = M_k_12_i ./ q_k_p;
M_k_p_21 = M_k_12_i .* q_k_p;
end

function [M_11, M_12, M_21, M_22] = matrix_product_2by2(M_1_11, M_1_12, ...
    M_1_21, M_1_22, M_2_11, M_2_12, M_2_21, M_2_22)
% matrix product of two 2x2 arrays
% M_ij: the matrix element in row i and column j
%
% M_1_ij: the element in row i and column j of the left matrix
% M_2_ij: the element in row i and column j of the right matrix

%%
M_11 = M_1_11 .* M_2_11 + M_1_12 .* M_2_21;
M_12 = M_1_11 .* M_2_12 + M_1_12 .* M_2_22;
M_21 = M_1_21 .* M_2_11 + M_1_22 .* M_2_21;
M_22 = M_1_21 .* M_2_12 + M_1_22 .* M_2_22;
end
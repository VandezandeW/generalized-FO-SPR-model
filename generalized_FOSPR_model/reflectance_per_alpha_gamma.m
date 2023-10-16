function [R_N, R_N_plus_1] = reflectance_per_alpha_gamma(alpha, gamma, ...
    r_s, r_p, N)
% Computes reflectance after N (R_N) and N+1 (R_N_1) reflections of a ray
% with reflection angle alpha, skew angle gamma, and reflection
% coefficients r_s and r_p.

% R_N: the reflectance after N reflections
% R_N_plus_1: the reflectance after N+1 reflections

% alpha: the reflection angle of the reflecting rays; [rad]
% gamma: the skew angle of the reflecting rays; [rad]
% r_s: the reflection coefficient for the electric field perpendicular
%   (senkrecht) to the plane of incidence
% r_p: the reflection coefficient for the electric field parallel to the
% 	plane of incidence
% N: number of reflections

% All inputs are in order [alpha gamma {s,p,k} {s',p',k'} lambda] with
% {s,p,k} the basis of the incoming ray at the reflection and {s',p',k'}
% the basis of the outgoing ray at the reflection.

%% meridional angle theta (eq. 15)
[alpha, gamma] = ndgrid(alpha, gamma);
theta = real(asin(cos(alpha) ./ cos(gamma)));

%% reorder dimensions
% from [alpha gamma {s,p,k} {s',p',k'} lambda] to
% [{s,p,k} {s',p',k'} lambda alpha gamma]
theta = permute(theta, [3 4 5 1 2]);
gamma = permute(gamma, [3 4 5 1 2]);
r_s = permute(r_s, [4 5 1 2 3]);
r_p = permute(r_p, [4 5 1 2 3]);
N = permute(N, [3 4 5 1 2]);

%% reflectance R_N and R_N_1
% initial propagation vector of reflecting light (eq. S39)
k_0 = [cos(gamma) .* sin(theta); ...
    sin(gamma) .* sin(theta); ...
    cos(theta)];

% the right-hand rule rotation matrix around the axis of the optical fiber
% (z-axis) for the sequential propagation vector (eq. 17)
% example: k_1 = k_rotation_matrix * k_0
% Direction of rotation is arbitrary. Due to symmetry, both directions are
% equivalent.
zeros_size_gamma = zeros(size(gamma));
ones_size_gamma = ones(size(gamma));
k_rotation_matrix = [...
    [-cos(2 * gamma), -sin(2 * gamma), zeros_size_gamma]; ...
    [sin(2 * gamma), -cos(2 * gamma), zeros_size_gamma]; ...
    [zeros_size_gamma, zeros_size_gamma, ones_size_gamma] ...
    ];

% Jones matrix of reflection (eq. 23)
% J = [[r_s, 0, 0];[0, r_p, 0];[0, 0, 1]];
zeros_size_r = zeros(size(r_s));
ones_size_r = ones(size(r_s));
J = [[r_s, zeros_size_r, zeros_size_r]; ...
    [zeros_size_r, r_p, zeros_size_r]; ...
    [zeros_size_r, zeros_size_r, ones_size_r]];

% Initial electric field of two perpendicular polarized rays
% E_0_in_local_1 = [sin(polarization_angle); cos(polarization_angle); 0];
% (eq. 29)
E_0_in_local_1__0 = [1; 0; 0];% polarization_angle = 0
E_N_in_local_N_1__0 = electric_field_vector_N_reflections( ...
    E_0_in_local_1__0);

% (eq. S50)
% E_N_1_in_local_N_1__0(:,:,i,j,k) = J(:,:,i,j,k) * ...
%   E_N_in_local_N_1__0(:,:,i,j,k);
E_N_1_in_local_N_1__0 = pagemtimes(J, E_N_in_local_N_1__0);

E_0_in_local_1__pi_2 = [0; 1; 0];% polarization_angle = pi/2
E_N_in_local_N_1__pi_2 = electric_field_vector_N_reflections( ...
    E_0_in_local_1__pi_2);
% E_N_1_in_local_N_1__pi_2(:,:,i,j,k) = J(:,:,i,j,k) * ...
%   E_N_in_local_N_1__pi_2(:,:,i,j,k);% (eq. S50)
E_N_1_in_local_N_1__pi_2 = pagemtimes(J, E_N_in_local_N_1__pi_2);

% The reflectance of unpolarized light is the average of two
% perpendicularly polarized light beams. (eq. 30)
R_N = (sum(vecnorm(E_N_in_local_N_1__0).^2, 1) + ...
    sum(vecnorm(E_N_in_local_N_1__pi_2).^2, 1)) ./ 2;
R_N_plus_1 = (sum(vecnorm(E_N_1_in_local_N_1__0).^2, 1) + ...
    sum(vecnorm(E_N_1_in_local_N_1__pi_2).^2, 1)) ./ 2;

%% reorder dimensions to [lambda alpha gamma]
R_N = permute(R_N, [3 4 5 1 2]);
R_N_plus_1 = permute(R_N_plus_1, [3 4 5 1 2]);

%%
    function E_N_in_local_N_1 = electric_field_vector_N_reflections( ...
            E_0_in_local_1)
        % Computes the electric field vector after the Nth reflection in
        % N+1th reflection's local coordinates from the initial electric
        % field vector with polarization angle polarization_angle. (eq. 26)
        
        % E_0_in_local_1: initial electric field in local coordinates
        %   ([s; p; k])
        
        % E_N_in_local_N_1: the electric field vector after the Nth
        %   reflection in N+1th reflection's local coordinates
        
        %% The polarization ray-tracing matrix of 1 reflection is Q_1.
        % propagation vectors after first and second reflection (eq. 17)
        % k_1 = k_rotation_matrix * k_0
        k_1 = pagemtimes(k_rotation_matrix, k_0);
        % k_2 = k_rotation_matrix * k_1
        k_2 = pagemtimes(k_rotation_matrix, k_1);
        
        % orthogonal matrices for transformation of electric fields
        % (see functions for more information)
        O_out_1 = ...
            propagation_vectors_to_outgoing_orthogonal_matrix(k_0, ...
            k_1);
        O_in_2_inv = ...
            propagation_vectors_to_incoming_orthogonal_matrix(k_1, ...
            k_2);
        
        % Polarization ray-tracing matrix Q, reflects ray in local
        % coordinates and transforms ray to the local-coordinates of the
        % next reflection. (eq. 25)
        % Q_1 = O_in_2_inv * O_out_1 * J;
        Q_1 = pagemtimes(O_in_2_inv, pagemtimes(O_out_1, J));
        
        %% The polarization ray-tracing matrix of N reflections is Q_1^N.
        % (eq. 26)
        % Reshape and initialize arrays for parfor-loop.
        size_Q_1 = size(Q_1);
        Q_1_temp_length = prod(size_Q_1(3:end), 'all');
        Q_1_temp = reshape(Q_1, [size_Q_1(1), size_Q_1(2), ...
            Q_1_temp_length]);
        N_temp = repmat(N, [1, 1, size(Q_1, 3), 1]);
        N_temp = reshape(N_temp, [1, 1, Q_1_temp_length]);
        % if N == 0 than Q_1_N = 3x3 identity matrix
        Q_1_N_temp = repmat(eye(3, 3), [1, 1, Q_1_temp_length]);
        
        % Remove rays with no reflections.
        indices = 1:Q_1_temp_length;
        indices_non_zero = indices(N_temp ~= 0);
        N_temp_par = N_temp(:, :, indices_non_zero);
        % [[a, b, 0]; [c, d, 0];[0, 0, 1]^N =
        %   [[[a, b]; [c, d]]^N, [0; 0]]; [0, 0, 1]]
        Q_1_temp_par = Q_1_temp(1:2, 1:2, indices_non_zero);
        Q_1_N_temp_par = Q_1_N_temp(1:2, 1:2, indices_non_zero);
        
        % Q_1_N = Q_1^N;
        %         for ii = 1:length(indices_non_zero)% for debugging
        parfor ii = 1:length(indices_non_zero)
            Q_1_N_temp_par(:, :, ii) = ...
                Q_1_temp_par(:, :, ii)^N_temp_par(:, :, ii);
        end
        Q_1_N_temp(1:2, 1:2, indices_non_zero) = Q_1_N_temp_par;
        
        % Reshape back to original array size.
        Q_1_N = reshape(Q_1_N_temp, size_Q_1);
        
        %% Electric field after N reflections
        % in the local coordinates of the unperformed reflection N+1.
        % (eq. 26)
        % E_N_in_local_N_1 = Q_1_N * E_0_in_local_1;
        E_N_in_local_N_1 = pagemtimes(Q_1_N, E_0_in_local_1);
    end
end

function [O_out] = propagation_vectors_to_outgoing_orthogonal_matrix(...
    k_in, k_out)
% Function calculates the outgoing orthogonal matrix, O_out, from the
% incoming and outgoing propagation vectors, k_in and k_out, of the
% reflecting light.
%
% k_in: propagation vector of the incoming light
% k_out: propagation vector of the outgoing light
% O_out: the orthogonal matrix for the transformation from local
%   coordinates ([s; p; k]) to global coordinates ([x; y; z]) for the
%   electric field of the outgoing ray

%%
% basis state vector perpendicular to the plane of incidence and to both
% propagation vectors (eq. 18)
% (note: s_out = s_in; eq. 20)
s_out = cross(k_in, k_out, 1);
% s_out = s_out ./ norm(s_out);
s_out = s_out ./ sqrt(sum(s_out.^2, 1));

% basis state vector parallel to the plane of incidence and perpendicular
% to the outgoing propagation vector (eq. 21)
p_out = cross(k_out, s_out, 1);

% the outgoing orthogonal matrix (eq. 24)
O_out = [s_out, p_out, k_out];
end

function [O_in_inv] = propagation_vectors_to_incoming_orthogonal_matrix(...
    k_in, k_out)
% Function calculates the incoming orthogonal matrices, O_in_inv, from the
% incoming and outgoing propagation vectors, k_in and k_out, of the
% reflecting light
%
% k_in: propagation vector of the incoming light
% k_out: propagation vector of the outgoing light
% O_in_inv: the orthogonal matrix for the transformation from global
%   coordinates ([x; y; z]) to local coordinates ([s; p; k]) for the
%   electric field of the incoming ray

%%
% basis state vector perpendicular (senkrecht) to the plane of incidence
% and to both propagation vectors (eq. 18)
s_in = cross(k_in, k_out, 1);
s_in = s_in ./ sqrt(sum(s_in.^2, 1));

% basis state vector parallel to the plane of incidence and perpendicular
% to the incoming propagation vector (eq. 19)
p_in = cross(k_in, s_in, 1);

% the orthogonal matrices from transposed basis state vectors
% O_in_inv = [s_in'; p_in'; k_in']; (eq. 24)
transposed_order = [2 1 3:ndims(s_in)];
s_in_T = permute(s_in, transposed_order);
p_in_T = permute(p_in, transposed_order);
k_in_T = permute(k_in, transposed_order);
O_in_inv = [s_in_T; p_in_T; k_in_T];
end
function angular_pdf = generate_angular_pdf(RTMCS_name, varargin)
% angular_pdf: gridded data interpolant object of the angular probability
%   density function
%
% rtmcs: name of desired ray tracing Monte Carlo simulation
%
% optional:
% square_bin_number: number of bins desired in angle_pdf_fun for refleciton
%   angle alpha and skew angle gamma
% plot_bool: boolean, if a plot of the Monte Carlo angular distribution
%%
default_bin_number = 64;% with a RMSD of 0.00074328
%   compared to 512 for lambda = 550:900 and water. (32 has a RMSD of 
%   0.0020954.)
default_plot_bool = false;
p = inputParser;
addOptional(p, 'bin_number', default_bin_number, ...
    @(v)isscalar(v));
addOptional(p, 'plot_bool', default_plot_bool, @(v)(v == false || ...
    v == true));
parse(p, varargin{:});
bin_number = p.Results.bin_number;
plot_bool = p.Results.plot_bool;

%% load monte carlo simulation ray data
% output_structure = load('output_lens_bifurcated_20211013_1023.mat');
location_string = strcat(pwd,'\ray_tracing_Monte_Carlo_simulations\');
output_structure = load(strcat(location_string, RTMCS_name, '.mat'));
gamma = output_structure.gamma;
alpha = output_structure.alpha;

% Assumes that the rays after backreflection can be seen as a different
% ray.
gamma = [gamma(:, 1); gamma(:, 2)];
alpha = [alpha(:, 1); alpha(:, 2)];

%% determine angle boundaries for distribution
alpha_max = max(alpha);
alpha_min = min(alpha);
delta_alpha_range = alpha_max - alpha_min;
alpha_max = alpha_max + delta_alpha_range / bin_number;
alpha_min = alpha_min - delta_alpha_range / bin_number;
if alpha_max > pi/2
    alpha_max = pi/2;
end
if alpha_min < 0
    alpha_min = 0;
end

gamma_max = max(gamma);
gamma_min = min(gamma);
delta_gamma_range = gamma_max - gamma_min;
gamma_max = gamma_max + delta_gamma_range / bin_number;
gamma_min = gamma_min - delta_gamma_range / bin_number;
if gamma_max > pi/2
    gamma_max = pi/2;
end
if gamma_min < 0
    gamma_min = 0;
end

%% generate numerical angular probability density function
square_bin_step = 1/bin_number;

alpha_edges = alpha_min + (0:square_bin_step:1) .* ...
    (alpha_max - alpha_min);
gamma_edges = gamma_min + (0:square_bin_step:1) .* ...
    (gamma_max - gamma_min);

ray_pdf = histcounts2(alpha, gamma, ...
    alpha_edges, gamma_edges, ...
    'Normalization', 'pdf');
ray_pdf = ray_pdf ./ sum(ray_pdf, 'all');

%% generate interpolant object of the angular probability density function
% pad edges with zeros for proper extrapolation
ray_pdf_temp = zeros(size(ray_pdf) + [1, 1]);
ray_pdf_temp(2:end, 1:end-1) = ray_pdf;
ray_pdf = ray_pdf_temp;

alpha_centers = alpha_edges(1:end - 1) + diff(alpha_edges) ./ 2;
gamma_centers = gamma_edges(1:end - 1) + diff(gamma_edges) ./ 2;
alpha_centers = [2 * alpha_centers(1) - alpha_centers(2), alpha_centers];
gamma_centers = [gamma_centers, 2*gamma_centers(end) - ...
    gamma_centers(end - 1)];
[alpha_grid, gamma_grid] = ndgrid(alpha_centers, gamma_centers);

angular_pdf = griddedInterpolant(alpha_grid, gamma_grid, ray_pdf, ...
    'nearest');

%% plot density function
if plot_bool
    figure
    surf(alpha_grid, gamma_grid, ray_pdf, 'EdgeColor','none')
    xlabel('\alpha (rad)');
    ylabel('\gamma (rad)');
    title('Monte Carlo angular distribution')
    view(2)
    colorbar
    drawnow
end
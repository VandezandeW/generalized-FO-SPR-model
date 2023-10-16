function [ epsilon_clad ] = permittivity_clad( lambda )
% input lambda is wavelength in m

%% base on thorlabs data
%http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6845
% 436 nm: 1.467287
% 589.3 nm: 1.458965
% 1020 nm: 1.450703

lambda = lambda .* 1e6; % wavelength: m to µm
% 
A0 = 1.927393801542563;
A1 = 0.002076629715773;
A2 = 0.009323722026072;

epsilon_clad = A0 + A1 .* lambda.^2 + A2 .* lambda.^-2;
end
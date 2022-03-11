function [particle, fluid_u] = advanced_PTM(Flow, Particle, Model)
%ADVANCED_PTM - Model of particle motion using Lab or (SPINS) model input where:
%   Force = 1/2 rho_0 U^2 A C_d - where C_d is drag coefficient, U is
%   relative velocity, A is object area
%   du_dt = F/m , so A/m = 1/rho_f
% Inputs:
%    Flow - Structure containing flow data. An (x by t) matrix of observed fluid flow
%    (U_flow), Timestep, X_flow (a vector of x coordinates). Optionally
%    rho_0 (upper layer density)
%
%    Particle - Structure containing particle data. A double or vector of x
%    start locations (StartLoc), C_d (drag coefficient), rho_f (particle
%    density) and optionally start velocities for particles (StartU)
%
%    Model - Either 'advanced' or 'basic'    
%
% Outputs:
%    particle - Structure containing the particle's calculated x location, u and du/dt
%    at each timestep
%    fluid_u - Local fluid velocity at the particle location
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 19-Jan-2022; Last revision: 19-Jan-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
if isfield(Flow, 'rho_0')
    rho_1 = rho_0;
else
    rho_1 = 1029;
end
X = Flow.x;

x_start_ind = nearest_index(X, Particle.StartX);
timestep = Flow.timestep;
times = 0:timestep:size(Flow.U_flow, 2)-1;
u = Flow.U_flow;
% Initialise all the variables
n_particles = length(x_start_ind);
particle_u = zeros(length(times), n_particles);
particle_x = zeros(length(times), n_particles);
particle_dudt = zeros(length(times), n_particles);
fluid_u = zeros(length(times),n_particles);

% Set starting conditions
particle_x(1, :) = X(x_start_ind, 1);
for ii = 1:length(times)-1
    t = times(ii);
    fluid_u(ii, :) = interp1(X, u(:, ii), particle_x(ii, :));
    switch Model
        case 'basic'
            % Do a basic particle tracking (velocity = fluid velocity)
            particle_u = fluid_u;
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii, :)*timestep;
            % Advanced particle tracking (velocity = past velocity + acceleration)
        case 'advanced'
            if ii == 1
                if isfield('Particle', 'StartU')
                    particle_u(ii, :) = fluid_u(ii, :);% Initialise particles at local fluid velocity
                else
                    particle_u(ii, :) = u_0;
                end
                particle_coef = (0.5*C_d*rho_1/rho_f);
            end
            
            particle_dudt(ii+1, :) = particle_coef.*((fluid_u(ii, :)-particle_u(ii, :)).^2).*sign(fluid_u(ii, :)); % Update particle velocity
            particle_u(ii+1, :) = particle_u(ii, :)+particle_dudt(ii+1, :)*timestep; % Update Particle Velocity (u_t = t_(t-1) + du/dt)
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii+1, :)*timestep; % Update Particle time
    end
end

particle.u = particle_u;
particle.x = particle_x;
particle.dudt = particle_dudt;
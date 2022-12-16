function [Particle, fluid_u] = FloatMotionModel(Flow, Particle, Model)
%FLOATMOTIONMODEL - Model of particle motion using input flow data where:
%   Force = 1/2 rho_1 U^2 A C_d - where C_d is drag coefficient, U is
%   relative velocity, A is object area
%   du_dt = F/m , so A/m = 1/rho_f
%
% Inputs:
%    Flow - Structure containing flow data. A (x by t) matrix of observed fluid flow
%    (u_flow), timestep, x (a vector of x coordinates). Optionally
%    rho_1 (upper layer density)
%
%    Particle - Structure containing particle data. A double or vector of x
%    start locations (StartLoc), C_d (drag coefficient), rho_f (particle
%    density) and optionally start velocities for particles (StartU)
%
%    Model - Switch - either 'advanced' or 'basic'
%
% Outputs:
%    Particle - Structure containing the particle's calculated x location, u and du/dt
%    at each timestep
%    fluid_u - Local fluid velocity at the particle location
%
% Other m-files required: nearest_index
% Subfunctions: circle_weights
% MAT-files required: none
%
% See also: run_lab_ptm, run_model_ptm, find_cd
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 19-Jan-2022; Last revision: 28-Nov-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%% Parse variables
if isfield(Flow, 'rho_1')
    rho_1 = Flow.rho_1;
else
    rho_1 = 1029;
end

% Spatial
X = Flow.x;
dx = abs(X(2)-X(1));
x_start_ind = nearest_index(X, Particle.StartLoc);
if strcmpi(Particle.Shape, 'circle') % Calculate the weighting for each x
    weights = circle_weights(Particle.r, Flow.x);
else
    dx = abs(Flow.x(2)-Flow.x(1));
    circle_x = -(Particle.r+dx):dx:(Particle.r+dx);
    weights = ones(1, length(circle_x)-1);
end

% Temporal
dt = Flow.timestep;
times = (0:size(Flow.u_flow, 2)-1)*dt;

% Velocity
u = Flow.u_flow;

% Initialise all the variables
n_particles = length(x_start_ind);
particle_u = zeros(length(times), n_particles);
particle_x = zeros(length(times), n_particles);
particle_dudt = zeros(length(times), n_particles);
fluid_u = zeros(length(times),n_particles);
particle_x(1, :) = X(x_start_ind, 1); % Set starting condition

%% Run the model
for ii = 1:length(times)-1
    tmp_particle_x = particle_x(ii,:)' + (-(Particle.r+dx/2):dx:(Particle.r+dx/2));
    fluid_us = interp1(X, u(:, ii), tmp_particle_x);
    fluid_u(ii, :) = sum(fluid_us.*weights, 'omitnan')/sum(weights); % Calculate mean (weighted) flow velocity under float
    switch Model
        case 'basic'
            %% basic particle tracking (velocity = fluid velocity)
            particle_u = fluid_u;
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii, :)*dt;

        case 'advanced'
            %% Advanced particle tracking (velocity = past velocity + acceleration)
            if ii == 1
                if ~isfield('Particle', 'StartU')
                    particle_u(ii, :) = fluid_u(ii, :);% Initialise particles at local fluid velocity
                else
                    particle_u(ii, :) = Particle.StartU;
                end
                particle_coef = (0.5*Particle.C_d*rho_1/Particle.rho_f);
            end
            tmp_rel_us = (fluid_us-particle_u(ii, :)).*abs(fluid_us-particle_u(ii, :)); % Difference in velocity for each grid point
            tmp_force = particle_coef.*sum(tmp_rel_us.*weights, 'omitnan')/sum(weights); % Force for each grid point
            particle_dudt(ii+1, :) = tmp_force; % Update particle velocity
            particle_u(ii+1, :) = particle_u(ii, :)+particle_dudt(ii+1, :)*dt; % Update Particle Velocity (u_t = t_(t-1) + du/dt)
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii+1, :)*dt; % Update Particle time
    end

end % End of model

%% Produce output structure
Particle.u = particle_u;
Particle.x = particle_x;
if strcmp('Model', 'advanced')
    Particle.dudt = particle_dudt;
end
end % End of main function

function [weights] = circle_weights(r, x)
%CIRCLE_WEIGHTS - Calculates the weighting for a circle at points x/r
dx = abs(x(2)-x(1));
circarea = @(x, r) 2*r.*sin(acos(x/r));
circle_x = -(r+dx):dx:(r+dx);
weights = NaN(1, length(circle_x)-1);
for i = 1:length(circle_x)-1
    weights(i) = integral(@(x) circarea(x, r+dx), circle_x(i), circle_x(i+1));
end
end

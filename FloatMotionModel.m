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
% Velocity
u = Flow.u_flow;

% Spatial Domains
x = Flow.x;
dx = abs(x(2)-x(1));

if (x(2)-x(1))<0
    x = flip(x);
    u = flip(u);
end

x_start_ind = nearest_index(x, Particle.StartLoc);
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


% Initialise all the variables
n_particles = length(x_start_ind);
particle_u = zeros(length(times), n_particles);
particle_x = zeros(length(times), n_particles);
particle_dudt = zeros(length(times), n_particles);
fluid_u = zeros(length(times),n_particles);
particle_x(1, :) = x(x_start_ind, 1); % Set starting condition

[XX, TT] = ndgrid(x, times);
F = griddedInterpolant(XX, TT, u); % Set up a gridded interpolant, it's quicker for the repeated queries needed in Runga Kutta

%% Run the model
for ii = 1:length(times)-1
    tmp_particle_x = particle_x(ii,:)' + (-(Particle.r+dx/2):dx:(Particle.r+dx/2));
    fluid_us = interp1(x, u(:, ii), tmp_particle_x);
    fluid_u(ii, :) = sum(fluid_us.*weights, 'omitnan')/sum(weights); % Calculate mean (weighted) flow velocity under float

    switch Model
        case 'basic'
            %% basic particle tracking (velocity = fluid velocity)
            k1_particle_x = particle_x(ii,:)' + (-(Particle.r+dx/2):dx:(Particle.r+dx/2));
            [K1_X, ttimes] = ndgrid(k1_particle_x, times(ii));

            k1_fluids = F(K1_X, ttimes)';
            k1 = sum(k1_fluids.*weights, 'omitnan')/sum(weights);

            k2_particle_x = k1_particle_x + (k1*dt/2);
            [K2_X, ttimes] = ndgrid(k2_particle_x, times(ii));
            k2_fluids = F(K2_X, ttimes+dt/2)';
            k2 = sum(k2_fluids.*weights, 'omitnan')/sum(weights);
    
            k3_particle_x = k1_particle_x + (k2*dt/2);
            [K3_X, ttimes] = ndgrid(k3_particle_x, times(ii));
            k3_fluids = F(K3_X, ttimes+dt/2)';
            k3 = sum(k3_fluids.*weights, 'omitnan')/sum(weights);

            k4_particle_x = k1_particle_x + (k3*dt/2);
            [K4_X, ttimes] = ndgrid(k4_particle_x, times(ii));
            k4_fluids = F(K4_X, ttimes+dt/2)';
            k4 = sum(k4_fluids.*weights, 'omitnan')/sum(weights);

            particle_u(ii) = (k1+2*k2+2*k3+k4)/6;
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii, :)*dt;
            if ii == length(times)-1
                particle_u(ii+1) = k4;
            end
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

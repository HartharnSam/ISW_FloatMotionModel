function [particle, fluid_u] = FloatMotionModel(Flow, Particle, Model)
%FLOATMOTIONMODEL - Model of particle motion using input flow data where:
%   Force = 1/2 rho_0 U^2 A C_d - where C_d is drag coefficient, U is
%   relative velocity, A is object area
%   du_dt = F/m , so A/m = 1/rho_f
% Inputs:
%    Flow - Structure containing flow data. An (x by t) matrix of observed fluid flow
%    (U_flow), timestep, x (a vector of x coordinates). Optionally
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
% Other m-files required: nearest_index
% Subfunctions: none
% MAT-files required: none
%
% See also: run_lab_ptm, run_model_ptm
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 19-Jan-2022; Last revision: Nov-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
if isfield(Flow, 'rho_0')
    rho_1 = Flow.rho_0;
else
    rho_1 = 1029;
end
X = Flow.x;

x_start_ind = nearest_index(X, Particle.StartLoc);
timestep = Flow.timestep;
times = (0:size(Flow.U_flow, 2)-1)*timestep;
u = Flow.U_flow;
% Initialise all the variables
n_particles = length(x_start_ind);
particle_u = zeros(length(times), n_particles);
particle_x = zeros(length(times), n_particles);
particle_dudt = zeros(length(times), n_particles);
fluid_u = zeros(length(times),n_particles);
if strcmpi(Particle.Shape, 'circle')
    weights = circle_weights(Particle.r, Flow.x);
else
    dx = abs(Flow.x(2)-Flow.x(1));
    circle_x = -(Particle.r+dx):dx:(Particle.r+dx);
    weights = ones(1, length(circle_x)-1);
end
% Set starting conditions
particle_x(1, :) = X(x_start_ind, 1);
dx = abs(X(2)-X(1));
for ii = 1:length(times)-1
    %t = times(ii);
    particle_xs = particle_x(ii,:)' + (-(Particle.r+dx/2):dx:(Particle.r+dx/2));

    fluid_us = interp1(X, u(:, ii), particle_xs);
    fluid_u(ii, :) = sum(fluid_us.*weights, 'omitnan')/sum(weights);
    switch Model
        case 'basic'
            % Do a basic particle tracking (velocity = fluid velocity)
            particle_u = fluid_u;
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii, :)*timestep;
            % Advanced particle tracking (velocity = past velocity + acceleration)
        case 'advanced'
            if ii == 1
                if ~isfield('Particle', 'StartU')
                    particle_u(ii, :) = fluid_u(ii, :);% Initialise particles at local fluid velocity
                else
                    particle_u(ii, :) = Particle.StartU;
                end
                particle_coef = (0.5*Particle.C_d*rho_1/Particle.rho_f);
            end

            particle_dudt(ii+1, :) = particle_coef.*((fluid_u(ii, :)-particle_u(ii, :))*abs(fluid_u(ii, :)-particle_u(ii, :))); % Update particle velocity
            particle_u(ii+1, :) = particle_u(ii, :)+particle_dudt(ii+1, :)*timestep; % Update Particle Velocity (u_t = t_(t-1) + du/dt)
            particle_x(ii+1, :) = particle_x(ii, :)+particle_u(ii+1, :)*timestep; % Update Particle time
    end
    % Do a CFL timestepping check:
    %CFL = particle_u(ii+1)*timestep/dx;
    %if CFL>1
        %disp(['CFL = ', num2str(CFL), ' >1 : Sort our a better timestep?'])
    %end
end

particle.u = particle_u;
particle.x = particle_x;
if strcmp('Model', 'advanced')
    particle.dudt = particle_dudt;
end


end

function [weights] = circle_weights(r, x)
dx = abs(x(2)-x(1));
circarea = @(x, r) 2*r.*sin(acos(x/r));
circle_x = -(r+dx):dx:(r+dx);
weights = NaN(1, length(circle_x)-1);
for i = 1:length(circle_x)-1
    weights(i) = integral(@(x) circarea(x, r+dx), circle_x(i), circle_x(i+1));
    %circle_xs(i) = (circle_x(i)+circle_x(i+1))/ 2;
end
end

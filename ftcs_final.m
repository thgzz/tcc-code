clc; clear; close all
% Parameters
Lx = 1;         % Length of the domain in the x-direction
Ly = 1;         % Length of the domain in the y-direction
Nx = 50;        % Number of grid points in the x-direction
Ny = 50;        % Number of grid points in the y-direction
Nt = 500;
T = 1;        % Total simulation time
alpha = 0.01;   % Thermal diffusivity
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
% dt = T / Nt;
dt = 0.01;      % Time step size

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initialize temperature field
u = 37 * ones(Nx, Ny);
u_new = u;

% Set initial condition (e.g., a Gaussian distribution)
% u(:, :) = exp(-((X - Lx/2).^2 + (Y - Ly/2).^2) / (2*0.1^2));
current_time = 0;

% Define a source term function (e.g., a Gaussian source)
source_rate = 0.1;
source_amplitude = 10;
source_center_x = Lx / 2;
source_center_y = Ly / 2;
source_radius = 0.2;   % Radius of influence
source = @(x, y, t) source_amplitude + source_rate * current_time * ...
                      (1 - exp(-((x - source_center_x).^2 + (y - source_center_y).^2) / (2*source_radius^2)));




% Time-stepping loop with boundary conditions
for t = 0:dt:T
    current_time = t;
    % Apply boundary conditions
    u_new(1, :) = 37;  % Dirichlet boundary condition on the left
    u_new(:, 1) = 37;  % Dirichlet boundary condition on the bottom
    u_new(end, :) = 37;  % Dirichlet boundary condition on the right
    
    % Neumann boundary condition on the top (zero gradient)
    u_new(:, end) = u_new(:, end-1); 

    % Update interior points using the FTCS method with the source term
    for i = 2:Nx-1
        for j = 2:Ny-1
            % FTCS scheme for 2D heat equation with source term
            u_new(i, j) = u(i, j) + alpha * dt * ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2) + dt * source(x(i), y(j), t);
        end
    end
    
    % Swap u and u_new for the next time step
    u = u_new;
    
    % Visualization (optional)
    if mod(t, 0.01) == 0  % Plot every 0.01 seconds
        contourf(X, Y, u, 20, 'EdgeColor', 'none');
        colorbar;
        axis square;
        title(['Time: ' num2str(t)]);
        xlabel('X');
        ylabel('Y');
        drawnow;
    end
end

% Final temperature distribution
contourf(X, Y, u, 20, 'EdgeColor', 'none');
colorbar;
axis square;
title('Final Temperature Distribution');
xlabel('X');
ylabel('Y');

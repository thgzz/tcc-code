clc; clear; close all;
% Parameters
Lx = 1;         % Length of the domain in x-direction
Ly = 1;         % Length of the domain in y-direction
Nx = 50;        % Number of grid points in x-direction
Ny = 50;        % Number of grid points in y-direction
T = 0.1;        % Total simulation time
alpha = 0.01;   % Thermal diffusivity

% Discretization
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dt = 0.001;  % You can adjust the time step size for stability

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initialize temperature field
u = 37*ones(Nx, Ny);
u_new = u;

% Set initial condition (e.g., a Gaussian distribution)
% cx = Nx/2;
% cy = Ny/2;
% u(cx,cy) = 37;
u(:, :) = exp(-((X - Lx/2).^2 + (Y - Ly/2).^2) / (2*0.1^2));

% Time-stepping loop
for t = 0:dt:T
    % Apply boundary conditions (e.g., Dirichlet or Neumann)
    u_new(1, :) = 37;  % Example: Dirichlet boundary condition on left
    u_new(end, :) = 37;  % Example: Dirichlet boundary condition on right
    u_new(:, 1) = 37;  % Example: Dirichlet boundary condition on bottom
    u_new(:, end) = 37;  % Example: Dirichlet boundary condition on top
    
    % Update interior points using the finite difference scheme
    for i = 2:Nx-1
        for j = 2:Ny-1
            u_new(i, j) = u(i, j) + alpha * dt * ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2);
        end
    end
    
    % Swap u and u_new for the next time step
    u = u_new;
    
    % Visualization (optional)
%     if mod(t, 0.01) == 0  % Plot every 0.01 seconds
%         contourf(X, Y, u, 20, 'EdgeColor', 'none');
%         colorbar;
%         axis square;
%         title(['Time: ' num2str(t)]);
%         xlabel('X');
%         ylabel('Y');
%         drawnow;
%     end
end

% Final temperature distribution
contourf(X, Y, u, 20, 'EdgeColor', 'none');
colorbar;
axis square;
title('Final Temperature Distribution');
xlabel('X');
ylabel('Y');

clc; clear; close all
% Parameters
Lx = 0.15;         % Length of the domain in the x-direction
Ly = 0.15;         % Length of the domain in the y-direction
Nx = 500;        % Number of grid points in the x-direction
Ny = 500;        % Number of grid points in the y-direction
T = 60;        % Total simulation time
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);

k = 0.54;       % Thermal conductivity
rho = 1066;
c = 3763;
alpha = k / ( rho * c);
w_blood = 0.002; % Blood perfusion rate
C_blood = 3617; % Specific heat of blood
rho_b = 1050;

% Stability
% st = (alpha*dt)/dx^2;
m = 0.25;
D = (w_blood * rho_b * C_blood) / (rho * c);
% dt = m / (alpha/dx^2+D/4)
dt = (m*dx^2)/alpha;

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initialize temperature field
u = 37 * ones(Nx, Ny);
u_new = u;

% Define a source term function (e.g., a Gaussian source)
r = Lx / 2;
r0 = 0.01;
A = 1.3e6;
% y moves in x and x moves in y
source_center_x = Lx / 2;
source_center_y = Ly / 2;

% source = @(x, y, t) A * exp(-((x - source_center_x).^2 + (y - source_center_y).^2) / (r0^2));
% source = @(x, y, t) A * exp(-((x - source_center_x).^2 + (y - source_center_y).^2) / (r0^2)) * sin(2*pi*t);
source = @(x, y, t) A * exp( -( (x - r).^2 + (y - source_center_y).^2 ) / r0^2);

    % Apply boundary conditions
    u_new(1, :) = 37;  % Dirichlet boundary condition on the left (top)
    u_new(:, 1) = 37;  % Dirichlet boundary condition on the bottom (left)
    u_new(end, :) = 37;  % Dirichlet boundary condition on the right (bottom)
    
    % Neumann boundary condition on the top (zero gradient) (right)
    u_new(:, end) = u_new(:, end-1); 


% Time-stepping loop with boundary conditions
for t = 0:dt:T


    % Update interior points using the FTCS method with the source term
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Source
            source_term = source(x(i), y(j), t);
            % Pennes' equation with source term
            u_new(i, j) = u(i, j) + alpha * dt * ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2) + dt * ((source_term - w_blood * (u(i, j) - 37) * C_blood * rho_b) / ( c * rho )); 
        end
    end
    
    % Swap u and u_new for the next time step
    u = u_new;
    
    % Visualization (optional)
    % if mod(t, 0.1) == 0  % Plot every 0.01 seconds
    %     contourf(X, Y, u, 20, 'EdgeColor', 'none');
    %     colorbar;
    %     axis square;
    %     title(['Time: ' num2str(t)]);
    %     xlabel('X');
    %     ylabel('Y');
    %     drawnow;
    % end
end

% Teste
% surf(X, Y, u, 'EdgeColor', 'none');
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
% zlabel('Temperature', 'Interpreter', 'latex', 'FontSize', 14);
% title('2D Heat Equation Solution (FTCS Method)', 'Interpreter', 'latex', 'FontSize', 16);
% colormap('jet');
% colorbar;
% view(3);


% Final temperature distribution
contourf(X, Y, u, 20, 'EdgeColor', 'none');
colorbar;
axis square;
title('Final Temperature Distribution');
xlabel('X');
ylabel('Y');

% Define the circle parameters
% circle_radius = 0.05;
% circle_center_x = Lx / 2;
% circle_center_y = Ly / 2;

% % Draw the circle
% rectangle('Position', [circle_center_x - circle_radius, circle_center_y - circle_radius, 2*circle_radius, 2*circle_radius], ...
%     'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);

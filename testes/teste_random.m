% Parameters
Lx = 1;         % Length of the domain in the x-direction
Ly = 1;         % Length of the domain in the y-direction
Nx = 50;        % Number of grid points in the x-direction
Ny = 50;        % Number of grid points in the y-direction
T = 60;         % Total simulation time in seconds
alpha = 0.01;   % Thermal diffusivity
k = 0.01;       % Thermal conductivity
w_blood = 4000; % Blood perfusion rate
C_blood = 4000; % Specific heat of blood
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dt = 0.01;      % Time step size

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initialize temperature field
u = 37 * ones(Nx, Ny);
u_new = u;

% Define the source center and parameters
source_center_x = Lx / 2;
source_center_y = Ly / 2;
A = 1.3e6;
r0 = 3.1e-3;

% Function to compute the source term
source = @(x, y, t) A * exp((x - source_center_x).^2 + (y - source_center_y).^2) / (r0^2);

% Time-stepping loop
for t = 0:dt:T
    for i = 2:Nx-1
        for j = 2:Ny-1
            r_squared = (x(i) - source_center_x)^2 + (y(j) - source_center_y)^2;
            % Pennes' equation with the source term
            u_new(i, j) = u(i, j) + alpha * dt * ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2) + dt * (source(x(i), y(j), t) - w_blood * (u(i, j) - 37) / (C_blood * k));
        end
    end
    
    % Swap u and u_new for the next time step
    u = u_new;
end

% Final temperature distribution
contourf(X, Y, u, 20, 'EdgeColor', 'none');
colorbar;
axis square;
title('Final Temperature Distribution');
xlabel('X');
ylabel('Y');

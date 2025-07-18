clc; clear; close all
% Parameters
Lx = 0.15;         % Length of the domain in the x-direction
Ly = 0.15;         % Length of the domain in the y-direction
Nx = 500;        % Number of grid points in the x-direction
Ny = 500;        % Number of grid points in the y-direction
T = 180;        % Total simulation time
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
sigma = 0.34;
E = 3e3;
% y moves in x and x moves in y
source_center_x = Lx / 2;
source_center_y = Ly / 2;
source_width = 0.01;

source = @(x, y, t) sigma*(E^2/2) * exp(-((x - source_center_x).^2 + (y - source_center_y).^2) / (2*source_width^2));

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
                (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2) + (dt/(c*rho)) * (source_term - w_blood * (37 - u(i, j)) * C_blood * rho_b); 
        end
    end
    
    % Swap u and u_new for the next time step
    u = u_new;
    
end

save dados.mat

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
% contourf(X, Y, u, 20, 'EdgeColor', 'none');
% colorbar;
% axis square;
% title('Final Temperature Distribution');
% xlabel('X');
% ylabel('Y');

hfig2 = figure;
contourf(X, Y, u, 20, 'EdgeColor', 'none');
a = colorbar;
axis square;
a.Label.String = 'Temperatura (\circC)';
a.Ticks = 0:10:100;
clim([37 100]);
xlabel('Comprimento (m)')
ylabel('Largura (m)')
fname2 = 'myfigure2';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig2,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document

set(findall(hfig2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

% Define the circle parameters
circle_radius = 0.05;
circle_radius_2 = 0.02;
circle_center_x = Lx / 2;
circle_center_y = Ly / 2;

% Draw the circle
rectangle('Position', [circle_center_x - circle_radius, circle_center_y - circle_radius, 2*circle_radius, 2*circle_radius], ...
    'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);

rectangle('Position', [circle_center_x - circle_radius_2, circle_center_y - circle_radius_2, 2*circle_radius_2, 2*circle_radius_2], ...
    'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);


print(hfig2,fname2,'-dpdf','-vector','-bestfit')


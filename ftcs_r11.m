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
m = 0.25;
dt = (m*dx^2)/alpha;

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initialize temperature field
u2 = 37 * ones(Nx, Ny);
u_new = u2;

% Define a source term function (e.g., a Gaussian source)
r = Lx / 2;
r0 = 11e-3; % nano spread r0 = 10mm
A = 1.3e6;
% y moves in x and x moves in y
source_center_x = Lx / 2;
source_center_y = Ly / 2;

source = @(x, y, t) A * exp( -( (x - r).^2 + (y - source_center_y).^2 ) / r0^2);

    % Apply boundary conditions
    u_new(1, :) = 37;  % Dirichlet boundary condition on the left (top)
    u_new(:, end) = 37; % (right)
    u_new(end, :) = 37;  % Dirichlet boundary condition on the right (bottom)
    
    % Neumann boundary condition on the top (zero gradient) (right)
    u_new(:, 1) = u_new(:,2); % (left)

tic
% Time-stepping loop with boundary conditions
for t = 0:dt:T


    % Update interior points using the FTCS method with the source term
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Source
            source_term = source(x(i), y(j), t);
            % Pennes' equation with source term
            u_new(i, j) = u2(i, j) + alpha * dt * ((u2(i+1, j) - 2*u2(i, j) + u2(i-1, j)) / dx^2 + ...
                (u2(i, j+1) - 2*u2(i, j) + u2(i, j-1)) / dy^2) + (dt/(c*rho)) * (source_term - w_blood * (37 - u2(i, j)) * C_blood * rho_b); 
        end
    end
    
    % Swap u and u_new for the next time step
    u2 = u_new;

end
%%%% Salvar dados %%%%
save dados_nanor11.mat
toc
%%%% Plots %%%%

hfig1 = figure;
contourf(X, Y, u2, 20, 'EdgeColor', 'none');
a = colorbar;
axis square;
a.Label.String = 'Temperatura (\circC)';
a.Ticks = 0:10:100;
clim([37 100]);
xlabel('Comprimento (m)')
ylabel('Largura (m)')
fname = 'myfigure3';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig1,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document

set(findall(hfig1,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig1,'Position');
set(hfig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

circle_radius = 0.05;
circle_radius_2 = 0.02;
circle_center_x = Lx / 2;
circle_center_y = Ly / 2;

% Draw the circle
rectangle('Position', [circle_center_x - circle_radius, circle_center_y - circle_radius, 2*circle_radius, 2*circle_radius], ...
    'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);

rectangle('Position', [circle_center_x - circle_radius_2, circle_center_y - circle_radius_2, 2*circle_radius_2, 2*circle_radius_2], ...
    'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);


% print(hfig1,fname,'-dpdf','-painters')
print(hfig1,fname,'-dpdf','-painters','-bestfit')



clc;
clear;
tic;
% Steger-Warming flux vector splitting scheme applied with finite volume methods in a 2D supersonic
% inlet with geometry defined by 'g321x033uf.dat' and inviscid flow

% Read in data, create primary grid
% The following line reads in the ".dat" file as a 'table' variable
table = readtable('g321x033uf.dat',"VariableNamingRule","preserve");
fileID = fopen('g321x033uf.dat','r','n');

% Format structure of the first line in '.dat' file 
format_spec = '%*8c%d %*4c%d%*\n';
[grid_sizing] = fscanf(fileID, format_spec);
nx = grid_sizing(1); 
ny =grid_sizing(2);
fclose(fileID);

% The table2array function converts the "table" data to an "array" type that can be operated on
cell_loc = table2array(table);
[numrows,~] = size(cell_loc);
clear cells_table;
clear table;
x_hold = zeros(numrows,1);
y_hold = zeros(numrows,1);

% Assigning the nx and ny values from the ".dat" file to x and y column vectors
for i = 1:numrows
    x_hold(i,1) = 0.015*cell_loc(i,1);
    y_hold(i,1) = 0.015*cell_loc(i,2);
end

% Using the reshape function to create the [nx,ny] primary grids for both x and y values. A Z matrix is defined for plotting purposes with "mesh()" 
x_primary = reshape(x_hold,[nx,ny]);
y_primary = reshape(y_hold,[nx,ny]);
z = zeros(nx,ny);
clear x_hold;
clear y_hold;
clear numrows;

% Define physical properties and reference values
r_gas = 287; % [J/(kg K)]
c_p = 1005; % [J/(kg K)]
gamma = 1.400; % Ratio of specific heats
p_ref = 11664; % [Pa]
T_ref = 216.7; % [K]
M_inlet = 3.000; % Mach number at inlet
c_ref = 295.0; % Speed of sound at inlet conditions [m/s]
u_ref = M_inlet*c_ref; % Initial axial convective velocity
v_ref = 0; % Assumes flow is perfectly parallel to inlet orientation
rho_ref = p_ref/(r_gas*T_ref); % Free stream reference density
et_ref = p_ref/(rho_ref*(gamma-1))+0.5*((u_ref)^2+(v_ref)^2); % Reference total energy
Q_ref = [rho_ref; rho_ref*u_ref; rho_ref*u_ref; rho_ref*et_ref]; % Reference state vector
CFL = 0.5; % Courant-Friedrichs-Lewy number set for the domain
tolerance = 1E-4; % Covergence criteria for the L_inf norm
epsilon = 0.04; % Added value for tagged eigenvalues
step_max = 8000; % Maximum number of psuedo-time steps performed

% Develop for ghost cell grid
% Develop the ghost cells that cover the periphery of the 2D grid while maintaining spacing
% Since the ghost cells begin before the 2D grid and end after, the arrays are length nx+2, ny+2
x_ghost = zeros(size(x_primary,1)+2,size(x_primary,2)+2);
y_ghost = zeros(size(x_ghost));
z_ghost = zeros(size(x_ghost));

% Develop the ghost grid x components
for j = 1:ny
    for i = 1:nx
        if i == 1
            x_ghost(i,j) = -x_primary(2,j);
            x_ghost(i,j+1) = -x_primary(2,j);
            x_ghost(i,j+2) = -x_primary(2,j);
        elseif i == nx
            x_ghost(i,j) = x_primary(i-1,j);
            x_ghost(i,j+1) = x_primary(i-1,j);
            x_ghost(i,j+2) = x_primary(i-1,j);
            x_ghost(i+1,j) = x_primary(i,j);
            x_ghost(i+1,j+1) = x_primary(i,j);
            x_ghost(i+1,j+2) = x_primary(i,j);
            x_ghost(i+2,j) = x_primary(nx,j) + (x_primary(nx,j) - x_primary(nx-1,j));
            x_ghost(i+2,j+1) = x_primary(nx,j) + (x_primary(nx,j) - x_primary(nx-1,j));
            x_ghost(i+2,j+2) = x_primary(nx,j) + (x_primary(nx,j) - x_primary(nx-1,j));
        else
            x_ghost(i,j) = x_primary(i-1,j);
            x_ghost(i,j+1) = x_primary(i-1,j);
            x_ghost(i,j+2) = x_primary(i-1,j);
        end
    end
end

% Develop the ghost grid for the y components
for j = 1:ny
    for i = 1:nx
        if j == 1
           y_ghost(i,j) = y_primary(i,1) - (y_primary(i,2)-y_primary(i,1));
           y_ghost(i+1,j) = y_primary(i,1) - (y_primary(i,2)-y_primary(i,1));
           y_ghost(i+2,j) = y_primary(i,1) - (y_primary(i,2)-y_primary(i,1));
        elseif j == ny
           y_ghost(i,j) = y_primary(i,j-1);
           y_ghost(i+1,j) = y_primary(i,j-1);
           y_ghost(i+2,j) = y_primary(i,j-1);
           y_ghost(i,j+1) = y_primary(i,j);
           y_ghost(i+1,j+1) = y_primary(i,j);
           y_ghost(i+2,j+1) = y_primary(i,j);
           y_ghost(i,j+2) = y_primary(i,j) + (y_primary(i,j)-y_primary(i,j-1));
           y_ghost(i+1,j+2) = y_primary(i,j) + (y_primary(i,j)-y_primary(i,j-1));
           y_ghost(i+2,j+2) = y_primary(i,j) + (y_primary(i,j)-y_primary(i,j-1));
        else
           y_ghost(i,j) = y_primary(i,j-1);
           y_ghost(i+1,j) = y_primary(i,j-1);
           y_ghost(i+2,j) = y_primary(i,j-1);
        end
    end
end


% Calculate cell face areas and cell volumes
% In the Xi direction, the projected cell face area matrix is expected to
% be of size [nx+2,ny+1]. This accounts for the ghost grid being defined at
% the cell vertices and the cell areas not including the topmost row of the
% grid.
S_xi_x = zeros(nx+2,ny+1);
S_xi_y = zeros(nx+2,ny+1);

% Calculate the cell surface projection for fluxes passing through the Xi
% direction normal to y_eta
for i = 1:nx+2
    for j = 1:ny+1
        S_xi_x(i,j) = y_ghost(i,j+1) - y_ghost(i,j);
    end
end

% Calculate the cell surface projection for fluxes passing through the Xi
% direction from contributions normal to x_eta
for j = 1:ny+1
    for i = 1:nx+2
        S_xi_y(i,j) = -(x_ghost(i,j+1) - x_ghost(i,j));
    end
end

% The same is applied to calculated the cell face areas in the eta
% direction. This time the cell face area matrix is expected to be size
% [nx+1,ny+2]
S_eta_x = zeros(nx+1,ny+2);
S_eta_y = zeros(nx+1,ny+2);

% Calculate the cell face projection from contributions in the Eta direction
% based on the change in x
for j = 1:ny+2
    for i = 1:nx+1
        S_eta_y(i,j) = x_ghost(i+1,j) - x_ghost(i,j);
    end
end

% Calculate the cell face projection from contributions in the Eta direction
% based on the change in y
for j = 1:ny+2
    for i = 1:nx+1
        S_eta_x(i,j) = -(y_ghost(i+1,j) - y_ghost(i,j));
    end
end

% Build the projected cell face area in the Xi and Eta direction. This is
% calculating the hypotenuse of the right triangle created by the projected cell face
% components
S_xi = ((S_xi_x).^2 + (S_xi_y).^2).^(1/2);
S_eta = ((S_eta_x).^2 + (S_eta_y).^2).^(1/2);

% To calculate the cell volumes, we calculated the cell area in 2D. This
% assumes unit depth in the Zeta direction. To get the volumes, we
% calculate a determinant using the 4 cell vertices of each cell
V_cell = zeros(nx+1,ny+1);
for j = 1:ny+1
    for i = 1:nx+1
        V_cell(i,j) = 0.5.*((x_ghost(i+1,j+1)-x_ghost(i,j)).*(y_ghost(i,j+1)-y_ghost(i+1,j))-(x_ghost(i,j+1)-x_ghost(i+1,j)).*(y_ghost(i+1,j+1)-y_ghost(i,j)));
    end
end


% Develop grid for cell centered plotting
% For plotting the contours of the projected cell face areas and cell
% volumes, a cell centered grid is created by interpolating halfpoints on
% the ghost grid. These halfpoints give the cell face and cell centroids
% where the area and volume data can be stored
x_cell_centered = zeros(nx+1,ny+1);
y_cell_centered = zeros(nx+1,ny+1);
z_cell_centered = zeros(nx+1,ny+1);

% Develop x grid half points
for j = 1:ny+1
    for i = 1:nx+1
        x_cell_centered(i,j) = (x_ghost(i+1,j)-x_ghost(i,j))/2 + x_ghost(i,j);
    end
end

% Develop y grid half points
for i = 1:nx+1
    for j = 1:ny+1
        y_cell_centered(i,j) = (y_ghost(i,j+1)-y_ghost(i,j))/2 + y_ghost(i,j);
    end 
end


% Initialize Q state vector and apply initial conditions
Q(2:nx,2:ny,1) = rho_ref;
Q(2:nx,2:ny,2) = rho_ref*u_ref;
Q(2:nx,2:ny,3) = rho_ref*v_ref;
Q(2:nx,2:ny,4) = rho_ref*et_ref;
Q(1,1:ny+1,1) = rho_ref;
Q(1,1:ny+1,2) = rho_ref*u_ref;
Q(1,1:ny+1,3) = rho_ref*v_ref;
Q(1,1:ny+1,4) = rho_ref*et_ref;
Q(nx+1,1:ny+1,1) =  rho_ref;
Q(nx+1,1:ny+1,2) =  rho_ref*u_ref;
Q(nx+1,1:ny+1,3) =  rho_ref*v_ref;
Q(nx+1,1:ny+1,4) =  rho_ref*et_ref;

% Initialize Q_v primitive vectors with initial Q state vectors
Q_v(2:nx,2:ny,1) = (gamma-1).*Q(2:nx,2:ny,4)-0.5.*(gamma-1).*Q(2:nx,2:ny,1).*((Q(2:nx,2:ny,2)./Q(2:nx,2:ny,1)).^2+(Q(2:nx,2:ny,3)./Q(2:nx,2:ny,1)).^2);
Q_v(2:nx,2:ny,2) = Q(2:nx,2:ny,2)./(Q(2:nx,2:ny,1));
Q_v(2:nx,2:ny,3) = Q(2:nx,2:ny,3)./Q(2:nx,2:ny,1);
Q_v(2:nx,2:ny,4) = (1./(r_gas.*Q(2:nx,2:ny,1))).*((gamma-1).*Q(2:nx,2:ny,4)-0.5.*(gamma-1).*Q(2:nx,2:ny,1).*((Q(2:nx,2:ny,2)./Q(2:nx,2:ny,1)).^2+(Q(2:nx,2:ny,3)./Q(2:nx,2:ny,1)).^2));
Q_v(1,1:ny+1,1) = (gamma-1).*Q(1,1:ny+1,4)-0.5.*(gamma-1).*Q(1,1:ny+1,1).*((Q(1,1:ny+1,2)./Q(1,1:ny+1,1)).^2+(Q(1,1:ny+1,3)./Q(1,1:ny+1,1)).^2);
Q_v(1,1:ny+1,2) = Q(1,1:ny+1,2)./(Q(1,1:ny+1,1));
Q_v(1,1:ny+1,3) = Q(1,1:ny+1,3)./(Q(1,1:ny+1,1));
Q_v(1,1:ny+1,4) = (1./(r_gas.*Q(1,1:ny+1,1))).*((gamma-1).*Q(1,ny+1,4)-0.5.*(gamma-1).*Q(1,1:ny+1,1).*((Q(1,ny+1,2)./Q(1,1:ny+1,1)).^2+(Q(1,1:ny+1,3)./Q(1,1:ny+1,1)).^2));
Q_v(nx+1,1:ny+1,1) = (gamma-1).*Q(nx+1,1:ny+1,4)-0.5.*(gamma-1).*Q(nx+1,1:ny+1,1).*((Q(nx+1,1:ny+1,2)./Q(nx+1,1:ny+1,1)).^2+(Q(nx+1,1:ny+1,3)./Q(nx+1,1:ny+1,1)).^2);
Q_v(nx+1,1:ny+1,2) = Q(nx+1,1:ny+1,2)./(Q(nx+1,1:ny+1,1));
Q_v(nx+1,1:ny+1,3) = Q(nx+1,1:ny+1,3)./(Q(nx+1,1:ny+1,1));
Q_v(nx+1,1:ny+1,4) = (1./(r_gas.*Q(nx+1,1:ny+1,1))).*((gamma-1).*Q(nx+1,ny+1,4)-0.5.*(gamma-1).*Q(nx+1,1:ny+1,1).*((Q(nx+1,ny+1,2)./Q(nx+1,1:ny+1,1)).^2+(Q(nx+1,1:ny+1,3)./Q(nx+1,1:ny+1,1)).^2));

% Calculate Q_v primitive vectors in the boundary halo cells (2:nx,1)
% (2:nx,ny+1)
for i = 2:nx
    % Bottom wall BC
    Q_v(i,1,1) = (gamma-1).*Q(i,2,4)-0.5.*(gamma-1).*Q(i,2,1).*((Q(i,2,2)./Q(i,2,1)).^2+(Q(i,2,3)./Q(i,2,1)).^2);
    Q_v(i,1,2) = (((S_eta_y(i,2))^2+(S_eta_x(i,2))^2)^(-1))*(Q(i,2,2)*((S_eta_y(i,2))^2-(S_eta_x(i,2))^2)+Q(i,2,3)*(-2*S_eta_y(i,2)*S_eta_x(i,2)))/Q(i,2,1);
    Q_v(i,1,3) = (((S_eta_y(i,2))^2+(S_eta_x(i,2))^2)^(-1))*(Q(i,2,2)*(-2*S_eta_y(i,2)*S_eta_x(i,2))+Q(i,2,3)*(-(S_eta_y(i,2))^2+(S_eta_x(i,2))^2))/Q(i,2,1);
    Q_v(i,1,4) = (1./(r_gas.*Q(i,2,1))).*((gamma-1).*Q(i,2,4)-0.5.*(gamma-1).*Q(i,2,1).*((Q(i,2,2)./Q(i,2,1)).^2+(Q(i,2,3)./Q(i,2,1)).^2));

    % Top wall BC
    Q_v(i,ny+1,1) = (gamma-1).*Q(i,ny,4)-0.5.*(gamma-1).*Q(i,ny,1).*((Q(i,ny,2)./Q(i,ny,1)).^2+(Q(i,ny,3)./Q(i,ny,1)).^2);
    Q_v(i,ny+1,2) = (((S_eta_y(i,ny))^2+(S_eta_x(i,ny))^2)^(-1))*(Q(i,ny,2)*((S_eta_y(i,ny))^2-(S_eta_x(i,ny))^2)+Q(i,ny,3)*(-2*S_eta_y(i,ny)*S_eta_x(i,ny)))/Q(i,ny,1);
    Q_v(i,ny+1,3) = (((S_eta_y(i,ny))^2+(S_eta_x(i,ny))^2)^(-1))*(Q(i,ny,2)*(-2*S_eta_y(i,ny)*S_eta_x(i,ny))+Q(i,ny,3)*(-(S_eta_y(i,ny))^2+(S_eta_x(i,ny))^2))/Q(i,ny,1);
    Q_v(i,ny+1,4) = (1./(r_gas.*Q(i,ny,1))).*((gamma-1).*Q(i,ny,4)-0.5.*(gamma-1).*Q(i,ny,1).*((Q(i,ny,2)./Q(i,ny,1)).^2+(Q(i,ny,3)./Q(i,ny,1)).^2));
end

function [Q_v] = Q_v_new(nx,ny,gamma,Q_new,r_gas,S_eta_x,S_eta_y)
    Q_v = zeros(nx+1,ny+1,4);
    % Q_v domain internal values and inlet/exit
    Q_v(2:nx,2:ny,1) = (gamma-1).*Q_new(2:nx,2:ny,4)-0.5.*(gamma-1).*Q_new(2:nx,2:ny,1).*((Q_new(2:nx,2:ny,2)./Q_new(2:nx,2:ny,1)).^2+(Q_new(2:nx,2:ny,3)./Q_new(2:nx,2:ny,1)).^2);
    Q_v(2:nx,2:ny,2) = Q_new(2:nx,2:ny,2)./(Q_new(2:nx,2:ny,1));
    Q_v(2:nx,2:ny,3) = Q_new(2:nx,2:ny,3)./Q_new(2:nx,2:ny,1);
    Q_v(2:nx,2:ny,4) = (1./(r_gas.*Q_new(2:nx,2:ny,1))).*((gamma-1).*Q_new(2:nx,2:ny,4)-0.5.*(gamma-1).*Q_new(2:nx,2:ny,1).*((Q_new(2:nx,2:ny,2)./Q_new(2:nx,2:ny,1)).^2+(Q_new(2:nx,2:ny,3)./Q_new(2:nx,2:ny,1)).^2));
    Q_v(1,1:ny+1,1) = (gamma-1).*Q_new(1,1:ny+1,4)-0.5.*(gamma-1).*Q_new(1,1:ny+1,1).*((Q_new(1,1:ny+1,2)./Q_new(1,1:ny+1,1)).^2+(Q_new(1,1:ny+1,3)./Q_new(1,1:ny+1,1)).^2);
    Q_v(1,1:ny+1,2) = Q_new(1,1:ny+1,2)./(Q_new(1,1:ny+1,1));
    Q_v(1,1:ny+1,3) = Q_new(1,1:ny+1,3)./(Q_new(1,1:ny+1,1));
    Q_v(1,1:ny+1,4) = (1./(r_gas.*Q_new(1,1:ny+1,1))).*((gamma-1).*Q_new(1,ny+1,4)-0.5.*(gamma-1).*Q_new(1,1:ny+1,1).*((Q_new(1,ny+1,2)./Q_new(1,1:ny+1,1)).^2+(Q_new(1,1:ny+1,3)./Q_new(1,1:ny+1,1)).^2));
    Q_v(nx+1,1:ny+1,1) = (gamma-1).*Q_new(nx+1,1:ny+1,4)-0.5.*(gamma-1).*Q_new(nx+1,1:ny+1,1).*((Q_new(nx+1,1:ny+1,2)./Q_new(nx+1,1:ny+1,1)).^2+(Q_new(nx+1,1:ny+1,3)./Q_new(nx+1,1:ny+1,1)).^2);
    Q_v(nx+1,1:ny+1,2) = Q_new(nx+1,1:ny+1,2)./(Q_new(nx+1,1:ny+1,1));
    Q_v(nx+1,1:ny+1,3) = Q_new(nx+1,1:ny+1,3)./(Q_new(nx+1,1:ny+1,1));
    Q_v(nx+1,1:ny+1,4) = (1./(r_gas.*Q_new(nx+1,1:ny+1,1))).*((gamma-1).*Q_new(nx+1,ny+1,4)-0.5.*(gamma-1).*Q_new(nx+1,1:ny+1,1).*((Q_new(nx+1,ny+1,2)./Q_new(nx+1,1:ny+1,1)).^2+(Q_new(nx+1,1:ny+1,3)./Q_new(nx+1,1:ny+1,1)).^2));
    
    % Halo cell values
    for i = 2:nx
        % Bottom wall BC
        Q_v(i,1,1) = (gamma-1).*Q_new(i,2,4)-0.5.*(gamma-1).*Q_new(i,2,1).*((Q_new(i,2,2)./Q_new(i,2,1)).^2+(Q_new(i,2,3)./Q_new(i,2,1)).^2);
        Q_v(i,1,2) = (((S_eta_y(i,2))^2+(S_eta_x(i,2))^2)^(-1))*(Q_new(i,2,2)*((S_eta_y(i,2))^2-(S_eta_x(i,2))^2)+Q_new(i,2,3)*(-2*S_eta_y(i,2)*S_eta_x(i,2)))/Q_new(i,2,1);
        Q_v(i,1,3) = (((S_eta_y(i,2))^2+(S_eta_x(i,2))^2)^(-1))*(Q_new(i,2,2)*(-2*S_eta_y(i,2)*S_eta_x(i,2))+Q_new(i,2,3)*(-(S_eta_y(i,2))^2+(S_eta_x(i,2))^2))/Q_new(i,2,1);
        Q_v(i,1,4) = (1./(r_gas.*Q_new(i,2,1))).*((gamma-1).*Q_new(i,2,4)-0.5.*(gamma-1).*Q_new(i,2,1).*((Q_new(i,2,2)./Q_new(i,2,1)).^2+(Q_new(i,2,3)./Q_new(i,2,1)).^2));
    
        % Top wall BC
        Q_v(i,ny+1,1) = (gamma-1).*Q_new(i,ny,4)-0.5.*(gamma-1).*Q_new(i,ny,1).*((Q_new(i,ny,2)./Q_new(i,ny,1)).^2+(Q_new(i,ny,3)./Q_new(i,ny,1)).^2);
        Q_v(i,ny+1,2) = (((S_eta_y(i,ny))^2+(S_eta_x(i,ny))^2)^(-1))*(Q_new(i,ny,2)*((S_eta_y(i,ny))^2-(S_eta_x(i,ny))^2)+Q_new(i,ny,3)*(-2*S_eta_y(i,ny)*S_eta_x(i,ny)))/Q_new(i,ny,1);
        Q_v(i,ny+1,3) = (((S_eta_y(i,ny))^2+(S_eta_x(i,ny))^2)^(-1))*(Q_new(i,ny,2)*(-2*S_eta_y(i,ny)*S_eta_x(i,ny))+Q_new(i,ny,3)*(-(S_eta_y(i,ny))^2+(S_eta_x(i,ny))^2))/Q_new(i,ny,1);
        Q_v(i,ny+1,4) = (1./(r_gas.*Q_new(i,ny,1))).*((gamma-1).*Q_new(i,ny,4)-0.5.*(gamma-1).*Q_new(i,ny,1).*((Q_new(i,ny,2)./Q_new(i,ny,1)).^2+(Q_new(i,ny,3)./Q_new(i,ny,1)).^2));
    end
end

% Apply BC velocity slip condition to halo cells for Q at j = 1; j = ny+1
for i = 2:nx
    % Top wall halo cells
    Q(i,ny+1,1) = Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4));
    Q(i,ny+1,2) = (Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4)))*Q_v(i,ny+1,2);
    Q(i,ny+1,3) = (Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4)))*Q_v(i,ny+1,3);
    Q(i,ny+1,4) = (Q_v(i,ny+1,1)/(gamma-1))+0.5*(Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4)))*((Q_v(i,ny+1,2))^2+((Q_v(i,ny+1,3)))^2);

    % Bottom wall halo cells
    Q(i,1,1) = Q_v(i,1,1)/(r_gas*Q_v(i,1,4));
    Q(i,1,2) = (Q_v(i,1,1)/(r_gas*Q_v(i,1,4)))*Q_v(i,1,2);
    Q(i,1,3) = (Q_v(i,1,1)/(r_gas*Q_v(i,1,4)))*Q_v(i,1,3);
    Q(i,1,4) = (Q_v(i,1,1)/(gamma-1))+0.5*(Q_v(i,1,1)/(r_gas*Q_v(i,1,4)))*((Q_v(i,1,2))^2+((Q_v(i,1,3)))^2);
end

% Initialize flux vector matrices
E_xi = zeros(nx,ny-1,4);
F_eta = zeros(nx-1,ny,4);

% Function caluclating E fluxes using left Q values (Q(i)) at the i+1/2 cell face
function [E_left] = E_L(nx,ny,Q,S_xi_x,S_xi_y,S_xi,gamma)
    E_left = zeros(nx,ny-1,4);
    for i = 1:nx
        for j = 1:ny-1
            E_left(i,j,1) = Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1));
            E_left(i,j,2) = Q(i,j+1,2)^2/Q(i,j+1,1)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+(Q(i,j+1,2)*Q(i,j+1,3)/Q(i,j+1,1))*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)+(gamma-1)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))*(Q(i,j+1,4)-0.5*(Q(i,j+1,2)^2/Q(i,j+1,1)+Q(i,j+1,3)^2/Q(i,j+1,1)));
            E_left(i,j,3) = (Q(i,j+1,3)*Q(i,j+1,2)/Q(i,j+1,1))*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+(Q(i,j+1,3)^2/Q(i,j+1,1))*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))+(gamma-1)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*(Q(i,j+1,4)-0.5*(Q(i,j+1,2)^2/Q(i,j+1,1)+Q(i,j+1,3)^2/Q(i,j+1,1)));
            E_left(i,j,4) = (gamma*Q(i,j+1,4)-0.5*(gamma-1)*(Q(i,j+1,2)^2/Q(i,j+1,1)+Q(i,j+1,3)^2/Q(i,j+1,1)))*(Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))/Q(i,j+1,1)+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))/Q(i,j+1,1));
        end
    end
end

% Function caluclating E convective vectors using right Q values (Q(i+1)) at the i+1/2 cell face
function [E_right] = E_R(nx,ny,Q,S_xi_x,S_xi_y,S_xi,gamma)
    E_right = zeros(nx,ny-1,4);
    for i = 1:nx
        for j = 1:ny-1
            E_right(i,j,1) = Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1));
            E_right(i,j,2) = Q(i+1,j+1,2)^2/Q(i+1,j+1,1)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+(Q(i+1,j+1,2)*Q(i+1,j+1,3)/Q(i+1,j+1,1))*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)+(gamma-1)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))*(Q(i+1,j+1,4)-0.5*(Q(i+1,j+1,2)^2/Q(i+1,j+1,1)+Q(i+1,j+1,3)^2/Q(i+1,j+1,1)));
            E_right(i,j,3) = (Q(i+1,j+1,3)*Q(i+1,j+1,2)/Q(i+1,j+1,1))*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+(Q(i+1,j+1,3)^2/Q(i+1,j+1,1))*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))+(gamma-1)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*(Q(i+1,j+1,4)-0.5*(Q(i+1,j+1,2)^2/Q(i+1,j+1,1)+Q(i+1,j+1,3)^2/Q(i+1,j+1,1)));
            E_right(i,j,4) = (gamma*Q(i+1,j+1,4)-0.5*(gamma-1)*(Q(i+1,j+1,2)^2/Q(i+1,j+1,1)+Q(i+1,j+1,3)^2/Q(i+1,j+1,1)))*(Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))/Q(i+1,j+1,1)+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))/Q(i+1,j+1,1));
        end
    end
end

% Function caluclating F convective vectors using left Q values (Q(j)) at the j+1/2 cell face
function [F_left] = F_L(nx,ny,Q,S_eta_x,S_eta_y,S_eta,gamma)
    F_left = zeros(nx-1,ny,4);
    for i = 1:nx-1
        for j = 1:ny
            F_left(i,j,1) = Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1));
            F_left(i,j,2) = Q(i+1,j,2)^2/Q(i+1,j,1)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+(Q(i+1,j,2)*Q(i+1,j,3)/Q(i+1,j,1))*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)+(gamma-1)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*(Q(i+1,j,4)-0.5*(Q(i+1,j,2)^2/Q(i+1,j,1)+Q(i+1,j,3)^2/Q(i+1,j,1)));
            F_left(i,j,3) = (Q(i+1,j,3)*Q(i+1,j,2)/Q(i+1,j,1))*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+(Q(i+1,j,3)^2/Q(i+1,j,1))*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))+(gamma-1)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*(Q(i+1,j,4)-0.5*(Q(i+1,j,2)^2/Q(i+1,j,1)+Q(i+1,j,3)^2/Q(i+1,j,1)));
            F_left(i,j,4) = (gamma*Q(i+1,j,4)-0.5*(gamma-1)*(Q(i+1,j,2)^2/Q(i+1,j,1)+Q(i+1,j,3)^2/Q(i+1,j,1)))*(Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))/Q(i+1,j,1)+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))/Q(i+1,j,1));
        end
    end
end

% Function caluclating F convective vectors using right Q values (Q(j+1)) at the j+1/2 cell face
function [F_right] = F_R(nx,ny,Q,S_eta_x,S_eta_y,S_eta,gamma)
    F_right = zeros(nx-1,ny,4);
    for i = 1:nx-1
        for j = 1:ny
            F_right(i,j,1) = Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1));
            F_right(i,j,2) = Q(i+1,j+1,2)^2/Q(i+1,j+1,1)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+(Q(i+1,j+1,2)*Q(i+1,j+1,3)/Q(i+1,j+1,1))*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)+(gamma-1)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*(Q(i+1,j+1,4)-0.5*(Q(i+1,j+1,2)^2/Q(i+1,j+1,1)+Q(i+1,j+1,3)^2/Q(i+1,j+1,1)));
            F_right(i,j,3) = (Q(i+1,j+1,3)*Q(i+1,j+1,2)/Q(i+1,j+1,1))*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+(Q(i+1,j+1,3)^2/Q(i+1,j+1,1))*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))+(gamma-1)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*(Q(i+1,j+1,4)-0.5*(Q(i+1,j+1,2)^2/Q(i+1,j+1,1)+Q(i+1,j+1,3)^2/Q(i+1,j+1,1)));
            F_right(i,j,4) = (gamma*Q(i+1,j+1,4)-0.5*(gamma-1)*(Q(i+1,j+1,2)^2/Q(i+1,j+1,1)+Q(i+1,j+1,3)^2/Q(i+1,j+1,1)))*(Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))/Q(i+1,j+1,1)+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))/Q(i+1,j+1,1));
        end
    end
end              

% Function for local time matrix initialized to store cell-centered time-steps in the
% interior of the domain
function [tau] = local_timestep(nx,ny,CFL,Q,Q_v,S_xi_x,S_xi_y,S_xi,S_eta_x,S_eta_y,S_eta,gamma,V_cell)
    tau = zeros(nx-1,ny-1);
    for i = 1:nx-1
        for j = 1:ny-1
            tau(i,j) = CFL*min(((abs((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1))+(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((S_xi_x(i+1,j+1)/V_cell(i+1,j+1))^2+(S_xi_y(i+1,j+1)/V_cell(i+1,j+1))^2)^(1/2))^(-1), ...
            ((abs((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1))+(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((S_eta_x(i+1,j+1)/V_cell(i+1,j+1))^2+(S_eta_y(i+1,j+1)/V_cell(i+1,j+1))^2)^(1/2))^(-1));
        end
    end
end

% Function calculating (+) eigenvalues of A jacobian using left state values
function [eigenA_plus_left] = eigenA_P_L(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma,epsilon)
    eigen_hold = zeros(nx,ny-1,4);
    eigenA_plus_left = zeros(nx,ny-1,4,4);
    for i = 1:nx
        for j = 1:ny-1
            eigen_hold(i,j,1) = 0.5*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)-sqrt(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))+sqrt(((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)-(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,2) = 0.5*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)+sqrt(((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1))^2+epsilon^2));
            eigen_hold(i,j,3) = 0.5*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)+sqrt(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))+sqrt(((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)+(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,4) = 0.5*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)+sqrt(((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1))^2+epsilon^2));
        end
    end
    % Construct diagonalized eigenvalue matrix
    for i = 1:nx
        for j = 1:ny-1
            eigenA_plus_left(i,j,:,:) = [eigen_hold(i,j,1);eigen_hold(i,j,2);eigen_hold(i,j,3);eigen_hold(i,j,4)].*eye(4);
        end
    end
end

% Function calculating (-) eigenvalues of A jacobian using right state
% values
function [eigenA_minus_right] = eigenA_M_R(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma,epsilon)
    eigen_hold = zeros(nx,ny-1,4);
    eigenA_minus_right = zeros(nx,ny-1,4,4);
    for i = 1:nx
        for j = 1:ny-1
            eigen_hold(i,j,1) = 0.5*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)-sqrt(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))-sqrt(((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)-(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,2) = 0.5*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)-sqrt(((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1))^2+epsilon^2));
            eigen_hold(i,j,3) = 0.5*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)+sqrt(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))-sqrt(((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)+(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,4) = 0.5*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)-sqrt(((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1))^2+epsilon^2));
        end
    end
    % Construct diagonalized eigenvalue matrix
    for i = 1:nx
        for j = 1:ny-1
            eigenA_minus_right(i,j,:,:) = [eigen_hold(i,j,1);eigen_hold(i,j,2);eigen_hold(i,j,3);eigen_hold(i,j,4)].*eye(4);
        end
    end
end

% Function calculating right eigenvectors matrices of A jacobian using left state values
function [vectorA_right_left] = vectorA_R_L(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma)
    vectorA_right_left = zeros(nx,ny-1,4,4);
    for i = 1:nx
        for j = 1:ny-1
            vectorA_right_left(i,j,:,:) = [1,1,1,0;
                                           Q_v(i,j+1,2)-(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2), Q_v(i,j+1,2), Q_v(i,j+1,2)+(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2), S_xi_y(i+1,j+1)/S_xi(i+1,j+1);
                                           Q_v(i,j+1,3)-(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2), Q_v(i,j+1,3), Q_v(i,j+1,3)+(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2), -S_xi_x(i+1,j+1)/S_xi(i+1,j+1);
                                           gamma*Q_v(i,j+1,1)/((gamma-1)*Q(i,j+1,1))+0.5*(Q_v(i,j+1,2)^2+Q_v(i,j+1,3)^2)-((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)), 0.5*(Q_v(i,j+1,2)^2+Q_v(i,j+1,3)^2), gamma*Q_v(i,j+1,1)/((gamma-1)*Q(i,j+1,1))+0.5*(Q_v(i,j+1,2)^2+Q_v(i,j+1,3)^2)+((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1)), Q_v(i,j+1,2)*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)-Q_v(i,j+1,3)*S_xi_x(i+1,j+1)/S_xi(i+1,j+1)];
        end
    end
end

% Function calculating right eigenvectors matrices of A jacobian using right state values
function [vectorA_right_right] = vectorA_R_R(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma)
    vectorA_right_right = zeros(nx,ny-1,4,4);
    for i = 1:nx
        for j = 1:ny-1
            vectorA_right_right(i,j,:,:) = [1,1,1,0;
                                           Q_v(i+1,j+1,2)-(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), Q_v(i+1,j+1,2), Q_v(i+1,j+1,2)+(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), S_xi_y(i+1,j+1)/S_xi(i+1,j+1);
                                           Q_v(i+1,j+1,3)-(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), Q_v(i+1,j+1,3), Q_v(i+1,j+1,3)+(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), -S_xi_x(i+1,j+1)/S_xi(i+1,j+1);
                                           gamma*Q_v(i+1,j+1,1)/((gamma-1)*Q(i+1,j+1,1))+0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2)-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)), 0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2), gamma*Q_v(i+1,j+1,1)/((gamma-1)*Q(i+1,j+1,1))+0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2)+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1)), Q_v(i+1,j+1,2)*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)-Q_v(i+1,j+1,3)*S_xi_x(i+1,j+1)/S_xi(i+1,j+1)];
        end
    end
end

% Function calculating left eigenvector matrices of A jacobian using left
% state values
function [vectorA_left_left] = vectorA_L_L(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma)
    vectorA_left_left = zeros(nx,ny-1,4,4);
    for i = 1:nx
        for j = 1:ny-1
            vectorA_left_left(i,j,:,:) = [(1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((gamma-1)*(0.5*(Q_v(i,j+1,2)^2+Q_v(i,j+1,3)^2))+((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1))), (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((1-gamma)*Q_v(i,j+1,2)-((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*S_xi_x(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((1-gamma)*Q_v(i,j+1,3)-((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*(gamma-1);
                                        (1/((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))-(gamma-1)*(0.5*(Q_v(i,j+1,2)^2+Q_v(i,j+1,3)^2))), (1/((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((gamma-1)*Q_v(i,j+1,2)), (1/((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((gamma-1)*Q_v(i,j+1,3)), (1/((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*(1-gamma);
                                        (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((gamma-1)*(0.5*(Q_v(i,j+1,2)^2+Q_v(i,j+1,3)^2))-((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1))), (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((1-gamma)*Q_v(i,j+1,2)+((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*S_xi_x(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*((1-gamma)*Q_v(i,j+1,3)+((gamma*Q_v(i,j+1,1)/Q(i,j+1,1))^(1/2))*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i,j+1,1)/Q(i,j+1,1))))*(gamma-1);
                                        (S_xi(i+1,j+1)/S_xi_x(i+1,j+1))*(Q_v(i,j+1,3)-(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*((Q(i,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i,j+1,1))), S_xi_y(i+1,j+1)/S_xi(i+1,j+1), ((S_xi_y(i+1,j+1)/S_xi(i+1,j+1))^2-1)/(S_xi_x(i+1,j+1)/S_xi(i+1,j+1)), 0];
        end
    end
end

% Function calculating left eigenvector matrices of A jacobian using right
% state values
function [vectorA_left_right] = vectorA_L_R(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma)
    vectorA_left_right = zeros(nx,ny-1,4,4);
    for i = 1:nx
        for j = 1:ny-1
            vectorA_left_right(i,j,:,:) = [(1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*(0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2))+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1))), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,2)-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_xi_x(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,3)-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*(gamma-1);
                                        (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))-(gamma-1)*(0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2))), (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*Q_v(i+1,j+1,2)), (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*Q_v(i+1,j+1,3)), (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*(1-gamma);
                                        (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*(0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2))-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1))), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,2)+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_xi_x(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,3)+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_xi_y(i+1,j+1)/S_xi(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*(gamma-1);
                                        (S_xi(i+1,j+1)/S_xi_x(i+1,j+1))*(Q_v(i+1,j+1,3)-(S_xi_y(i+1,j+1)/S_xi(i+1,j+1))*((Q(i+1,j+1,2)*(S_xi_x(i+1,j+1)/S_xi(i+1,j+1))+Q(i+1,j+1,3)*(S_xi_y(i+1,j+1)/S_xi(i+1,j+1)))/Q(i+1,j+1,1))), S_xi_y(i+1,j+1)/S_xi(i+1,j+1), ((S_xi_y(i+1,j+1)/S_xi(i+1,j+1))^2-1)/(S_xi_x(i+1,j+1)/S_xi(i+1,j+1)), 0];
        end
    end
end

% Function calculating (+) eigenvalues of B jacobian using left state values
function [eigenB_plus_left] = eigenB_P_L(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma,epsilon)
    eigen_hold = zeros(nx-1,ny,4);
    eigenB_plus_left = zeros(nx-1,ny,4,4);
    for i = 1:nx-1
        for j = 1:ny
            eigen_hold(i,j,1) = 0.5*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)-sqrt(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))+sqrt(((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)-(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,2) = 0.5*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)+sqrt(((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1))^2+epsilon^2));
            eigen_hold(i,j,3) = 0.5*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)+sqrt(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))+sqrt(((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)+(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,4) = 0.5*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)+sqrt(((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1))^2+epsilon^2));
        end
    end
    % Construct diagonalized eigenvalue matrix
    for i = 1:nx-1
        for j = 1:ny
            eigenB_plus_left(i,j,:,:) = [eigen_hold(i,j,1);eigen_hold(i,j,2);eigen_hold(i,j,3);eigen_hold(i,j,4)].*eye(4);
        end
    end
end

% Function calculating (-) eigenvalues of B jacobian using right state
% values
function [eigenB_minus_right] = eigenB_M_R(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma,epsilon)
    eigen_hold = zeros(nx-1,ny,4);
    eigenB_minus_right = zeros(nx-1,ny,4,4);
    for i = 1:nx-1
        for j = 1:ny
            eigen_hold(i,j,1) = 0.5*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)-sqrt(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))-sqrt(((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)-(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,2) = 0.5*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)-sqrt(((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1))^2+epsilon^2));
            eigen_hold(i,j,3) = 0.5*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)+sqrt(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))-sqrt(((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)+(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))^2+epsilon^2));
            eigen_hold(i,j,4) = 0.5*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)-sqrt(((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1))^2+epsilon^2));
        end
    end
    % Construct diagonalized eigenvalue matrix
    for i = 1:nx-1
        for j = 1:ny
            eigenB_minus_right(i,j,:,:) = [eigen_hold(i,j,1);eigen_hold(i,j,2);eigen_hold(i,j,3);eigen_hold(i,j,4)].*eye(4);
        end
    end
end

% Function calculating right eigenvectors matrices of B jacobian using left state values
function [vectorB_right_left] = vectorB_R_L(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma)
    vectorB_right_left = zeros(nx-1,ny,4,4);
    for i = 1:nx-1
        for j = 1:ny
            vectorB_right_left(i,j,:,:) = [1,1,1,0;
                                           Q_v(i+1,j,2)-(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2), Q_v(i+1,j,2), Q_v(i+1,j,2)+(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2), S_eta_y(i+1,j+1)/S_eta(i+1,j+1);
                                           Q_v(i+1,j,3)-(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2), Q_v(i+1,j,3), Q_v(i+1,j,3)+(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2), -S_eta_x(i+1,j+1)/S_eta(i+1,j+1);
                                           gamma*Q_v(i+1,j,1)/((gamma-1)*Q(i+1,j,1))+0.5*(Q_v(i+1,j,2)^2+Q_v(i+1,j,3)^2)-((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)), 0.5*(Q_v(i+1,j,2)^2+Q_v(i+1,j,3)^2), gamma*Q_v(i+1,j,1)/((gamma-1)*Q(i+1,j,1))+0.5*(Q_v(i+1,j,2)^2+Q_v(i+1,j,3)^2)+((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1)), Q_v(i+1,j,2)*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)-Q_v(i+1,j,3)*S_eta_x(i+1,j+1)/S_eta(i+1,j+1)];
        end
    end
end

% Function calculating right eigenvectors matrices of B jacobian using right state values
function [vectorB_right_right] = vectorB_R_R(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma)
    vectorB_right_right = zeros(nx-1,ny,4,4);
    for i = 1:nx-1
        for j = 1:ny
            vectorB_right_right(i,j,:,:) = [1,1,1,0;
                                           Q_v(i+1,j+1,2)-(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), Q_v(i+1,j+1,2), Q_v(i+1,j+1,2)+(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), S_eta_y(i+1,j+1)/S_eta(i+1,j+1);
                                           Q_v(i+1,j+1,3)-(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), Q_v(i+1,j+1,3), Q_v(i+1,j+1,3)+(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2), -S_eta_x(i+1,j+1)/S_eta(i+1,j+1);
                                           gamma*Q_v(i+1,j+1,1)/((gamma-1)*Q(i+1,j+1,1))+0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2)-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)), 0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2), gamma*Q_v(i+1,j+1,1)/((gamma-1)*Q(i+1,j+1,1))+0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2)+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1)), Q_v(i+1,j+1,2)*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)-Q_v(i+1,j+1,3)*S_eta_x(i+1,j+1)/S_eta(i+1,j+1)];
        end
    end
end

% Function calculating left eigenvector matrices of B jacobian using left
% state values
function [vectorB_left_left] = vectorB_L_L(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma)
    vectorB_left_left = zeros(nx-1,ny,4,4);
    for i = 1:nx-1
        for j = 1:ny
            vectorB_left_left(i,j,:,:) = [(1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((gamma-1)*(0.5*(Q_v(i+1,j,2)^2+Q_v(i+1,j,3)^2))+((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1))), (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((1-gamma)*Q_v(i+1,j,2)-((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*S_eta_x(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((1-gamma)*Q_v(i+1,j,3)-((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*(gamma-1);
                                        (1/((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))-(gamma-1)*(0.5*(Q_v(i+1,j,2)^2+Q_v(i+1,j,3)^2))), (1/((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((gamma-1)*Q_v(i+1,j,2)), (1/((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((gamma-1)*Q_v(i+1,j,3)), (1/((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*(1-gamma);
                                        (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((gamma-1)*(0.5*(Q_v(i+1,j,2)^2+Q_v(i+1,j,3)^2))-((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*((Q(i+1,j,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j,1))), (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((1-gamma)*Q_v(i+1,j,2)+((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*S_eta_x(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*((1-gamma)*Q_v(i+1,j,3)+((gamma*Q_v(i+1,j,1)/Q(i+1,j,1))^(1/2))*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j,1)/Q(i+1,j,1))))*(gamma-1);
                                        (((S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*Q_v(i+1,j,2)+(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*Q_v(i+1,j,3))*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))-Q_v(i+1,j,2))/(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)), S_eta_y(i+1,j+1)/S_eta(i+1,j+1), -S_eta_x(i+1,j+1)/S_eta(i+1,j+1), 0];
        end
    end
end

% Function calculating left eigenvector matrices of B jacobian using right
% state values
function [vectorB_left_right] = vectorB_L_R(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma)
    vectorB_left_right = zeros(nx-1,ny,4,4);
    for i = 1:nx-1
        for j = 1:ny
            vectorB_left_right(i,j,:,:) = [(1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*(0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2))+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1))), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,2)-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_eta_x(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,3)-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*(gamma-1);
                                        (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))-(gamma-1)*(0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2))), (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*Q_v(i+1,j+1,2)), (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*Q_v(i+1,j+1,3)), (1/((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*(1-gamma);
                                        (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((gamma-1)*(0.5*(Q_v(i+1,j+1,2)^2+Q_v(i+1,j+1,3)^2))-((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*((Q(i+1,j+1,2)*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))+Q(i+1,j+1,3)*(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)))/Q(i+1,j+1,1))), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,2)+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_eta_x(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*((1-gamma)*Q_v(i+1,j+1,3)+((gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))^(1/2))*S_eta_y(i+1,j+1)/S_eta(i+1,j+1)), (1/(2*(gamma*Q_v(i+1,j+1,1)/Q(i+1,j+1,1))))*(gamma-1);
                                        (((S_eta_x(i+1,j+1)/S_eta(i+1,j+1))*Q_v(i+1,j+1,2)+(S_eta_y(i+1,j+1)/S_eta(i+1,j+1))*Q_v(i+1,j+1,3))*(S_eta_x(i+1,j+1)/S_eta(i+1,j+1))-Q_v(i+1,j+1,2))/(S_eta_y(i+1,j+1)/S_eta(i+1,j+1)), S_eta_y(i+1,j+1)/S_eta(i+1,j+1), -S_eta_x(i+1,j+1)/S_eta(i+1,j+1), 0];
        end
    end
end

% Iterate on discretized problem for convergence
% Initialize iteration counter
step = 1

% Initialize new conserved state vector
Q_new = zeros(nx+1,ny+1,4);

% Initialize error matrices
error_rho = zeros(nx+1,ny+1,1);
error_rho_u = zeros(nx+1,ny+1,1);
error_rho_v = zeros(nx+1,ny+1,1);
error_rho_et = zeros(nx+1,ny+1,1);


while step <= step_max
    % Construct local time-steps using Q_old and Q_v old
    % Local time-step performed for interior cells only
    tau = local_timestep(nx,ny,CFL,Q,Q_v,S_xi_x,S_xi_y,S_xi,S_eta_x,S_eta_y,S_eta,gamma,V_cell);

    % Eigenvalue and eigenvector matrices for A jacobian
    eigenA_plus_left = eigenA_P_L(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma,epsilon);
    eigenA_minus_right = eigenA_M_R(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma,epsilon);
    vectorA_right_left = vectorA_R_L(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma);
    vectorA_right_right = vectorA_R_R(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma);
    vectorA_left_left = vectorA_L_L(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma);
    vectorA_left_right = vectorA_L_R(nx,ny,Q,Q_v,S_xi_x,S_xi_y,S_xi,gamma);

    % E convective fluxes
    E_left = E_L(nx,ny,Q,S_xi_x,S_xi_y,S_xi,gamma);
    E_right = E_R(nx,ny,Q,S_xi_x,S_xi_y,S_xi,gamma);

    % Eigenvalue and eigenvector matrices for B jacobian
    eigenB_plus_left = eigenB_P_L(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma,epsilon);
    eigenB_minus_right = eigenB_M_R(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma,epsilon);
    vectorB_right_left = vectorB_R_L(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma);
    vectorB_right_right = vectorB_R_R(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma);
    vectorB_left_left = vectorB_L_L(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma);
    vectorB_left_right = vectorB_L_R(nx,ny,Q,Q_v,S_eta_x,S_eta_y,S_eta,gamma);

    % F covective fluxes
    F_left = F_L(nx,ny,Q,S_eta_x,S_eta_y,S_eta,gamma);
    F_right = F_R(nx,ny,Q,S_eta_x,S_eta_y,S_eta,gamma);

    % Construct finite volume E_x fluxes using split A Jacobians
    for i = 1:nx
        for j = 1:ny-1
            E_xi(i,j,:) = squeeze(vectorA_right_left(i,j,:,:))*squeeze(eigenA_plus_left(i,j,:,:))*squeeze(vectorA_left_left(i,j,:,:))*[Q(i,j+1,1);Q(i,j+1,2);Q(i,j+1,3);Q(i,j+1,4)] + squeeze(vectorA_right_right(i,j,:,:))*squeeze(eigenA_minus_right(i,j,:,:))*squeeze(vectorA_left_right(i,j,:,:))*[Q(i+1,j+1,1);Q(i+1,j+1,2);Q(i+1,j+1,3);Q(i+1,j+1,4)];
        end 
    end

    % Construct finite volume F_eta fluxes using split B Jacobians
    for i = 1:nx-1
        for j = 1:ny
                F_eta(i,j,:) = squeeze(vectorB_right_left(i,j,:,:))*squeeze(eigenB_plus_left(i,j,:,:))*squeeze(vectorB_left_left(i,j,:,:))*[Q(i+1,j,1);Q(i+1,j,2);Q(i+1,j,3);Q(i+1,j,4)] + squeeze(vectorB_right_right(i,j,:,:))*squeeze(eigenB_minus_right(i,j,:,:))*squeeze(vectorB_left_right(i,j,:,:))*[Q(i+1,j+1,1);Q(i+1,j+1,2);Q(i+1,j+1,3);Q(i+1,j+1,4)];
        end
    end 

    % Perform flux differencing in the Xi and Eta directions and solve for
    % new interior Q vector (Q_new)
    for i = 1:nx-1
        for j = 1:ny-1
            Q_new(i+1,j+1,:) = Q(i+1,j+1,:) - (tau(i,j)/V_cell(i+1,j+1))*(E_xi(i+1,j,:)*S_xi(i+2,j+1)-E_xi(i,j,:)*S_xi(i+1,j+1) + F_eta(i,j+1,:)*S_eta(i+1,j+2) - F_eta(i,j,:)*S_eta(i+1,j+1));
        end
    end
    
    % Enforce Inlet BC
    Q_new(1,1:ny+1,:) = Q(1,1:ny+1,:);

    % Enforce outlet
    Q_new(nx+1,1:ny+1,:) = Q(nx+1,1:ny+1,:);

    % Recalculate Q_v with Q_new
    Q_v = Q_v_new(nx,ny,gamma,Q_new,r_gas,S_eta_x,S_eta_y);

    % Apply BC velocity slip condition to halo cells at j = 1; j = ny+1
    for i = 2:nx
        % Top wall halo cells
        Q_new(i,ny+1,1) = Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4));
        Q_new(i,ny+1,2) = (Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4)))*Q_v(i,ny+1,2);
        Q_new(i,ny+1,3) = (Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4)))*Q_v(i,ny+1,3);
        Q_new(i,ny+1,4) = (Q_v(i,ny+1,1)/(gamma-1))+0.5*(Q_v(i,ny+1,1)/(r_gas*Q_v(i,ny+1,4)))*((Q_v(i,ny+1,2))^2+((Q_v(i,ny+1,3)))^2);
    
        % Bottom wall halo cells
        Q_new(i,1,1) = Q_v(i,1,1)/(r_gas*Q_v(i,1,4));
        Q_new(i,1,2) = (Q_v(i,1,1)/(r_gas*Q_v(i,1,4)))*Q_v(i,1,2);
        Q_new(i,1,3) = (Q_v(i,1,1)/(r_gas*Q_v(i,1,4)))*Q_v(i,1,3);
        Q_new(i,1,4) = (Q_v(i,1,1)/(gamma-1))+0.5*(Q_v(i,1,1)/(r_gas*Q_v(i,1,4)))*((Q_v(i,1,2))^2+((Q_v(i,1,3)))^2);
    end
    
    % Calculate maximum errors across grid for each conserved variable
    error_rho = max(max(abs(Q_new(:,:,1) - Q(:,:,1))/rho_ref));
    error_rho_u = max(max(abs(Q_new(:,:,2) - Q(:,:,2))/(rho_ref*u_ref)));
    error_rho_v = max(max(abs(Q_new(:,:,3) - Q(:,:,3))/(rho_ref*u_ref)));
    error_rho_et = max(max(abs(Q_new(:,:,4) - Q(:,:,4))/(rho_ref*et_ref)));
    
    % Calculate residual using most limiting variable error
    residual(step) = max([error_rho,error_rho_u,error_rho_v,error_rho_et]);
    current_residual = residual(step)
    density_resid(step) = error_rho;
    x_momentum_resid(step) = error_rho_u;
    y_momentum_resid(step) = error_rho_v;
    energ_resid(step) = error_rho_et;

    % Exit iterations on satisfying tolerance requirements
    if residual(step) <= tolerance
        break
    else
    end

    % Update Q
    Q = Q_new;
    step = step + 1
end


% Plotting
figure(1)
semilogy(residual,'-*','MarkerIndices',1:75:length(residual))
xlabel("Iterations")
ylabel("Residual")
title("L_\infty- Norm Convergence")
hold on
semilogy(density_resid)
semilogy(x_momentum_resid)
semilogy(y_momentum_resid)
semilogy(energ_resid)
yline(tolerance)
legend("Maximum Residual", "Density Residual", "X-Momentum Residual", "Y-Momentum Residual", "Energy Residual", "Tolerance","Location","northeast")
hold off

% Velocity
figure(2)
subplot(3,1,1)
contourf(x_cell_centered,y_cell_centered,Q_v(:,:,2)./u_ref,100,"LineStyle","none")
colormap("turbo")
cb = colorbar('northoutside');
title(cb,'u/U_r_e_f')
xlabel('x [m]')
ylabel('y [m]')

% Pressure
subplot(3,1,2)
contourf(x_cell_centered,y_cell_centered,(Q_v(:,:,1)-p_ref)./(rho_ref*u_ref^2),100,"LineStyle","none")
colormap("turbo")
cb = colorbar('northoutside');
title(cb,'(p-p_r_e_f)/(rho_r_e_f*U_r_e_f^2)')
xlabel('x [m]')
ylabel('y [m]')

% Temperature
subplot(3,1,3)
contourf(x_cell_centered,y_cell_centered,(Q_v(:,:,4)-T_ref)./(c_ref^2/c_p),100,"LineStyle","none")
colormap("turbo")
cb = colorbar('northoutside');
title(cb,'(T-T_r_e_f)/(c_r_e_f^2/c_p_,_r_e_f)')
xlabel('x [m]')
ylabel('y [m]')

% Density
figure(3)
subplot(2,1,1)
contourf(x_cell_centered,y_cell_centered,Q_new(:,:,1)./rho_ref,100,"LineStyle","none")
colormap("turbo")
cb = colorbar('northoutside');
title(cb,'\rho/\rho_r_e_f')
xlabel('x [m]')
ylabel('y [m]')

% Mach
subplot(2,1,2)
contourf(x_cell_centered,y_cell_centered,Q_v(:,:,2)./(sqrt(gamma*Q_v(:,:,1)./Q_new(:,:,1))),100,"LineStyle","none")
colormap("turbo")
cb = colorbar('northoutside');
title(cb,'Mach')
xlabel('x [m]')
ylabel('y [m]')
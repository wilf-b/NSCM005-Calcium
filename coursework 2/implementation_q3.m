% Set up parameters
params = struct();
params.k_flux = 8.1;
params.b = 0.111;
params.V1 = 0.889;
params.k1 = 0.7;
params.gamma = 2.0;
params.k_gamma = 0.1;
params.beta = 0.02;
params.tau_n = 2.0;
params.k2 = 0.7;
params.D_p = 300;
params.D_c = 20;
params.k_p = 0.1;
params.k_mu = 4.0;
params.mu0 = 0.567;
params.mu1 = 0.433;

% Set up spatial grid
Nx = 100;
Ny = 100;
dx = 5;
dy = 5;
x = linspace(0, Nx*dx, Nx);
y = linspace(0, Ny*dy, Ny);
[X, Y] = meshgrid(x, y);

% Set up initial conditions
c0 = 0.1 * ones(Nx, Ny);
n0 = 0.1 * ones(Nx, Ny);
P0 = zeros(Nx, Ny);

% Create initial condition vector
y0 = [c0(:); n0(:); P0(:)];

% Set up time span
tspan = [0 100];

% Solve the PDE system
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
[t, y] = ode15s(@(t, y) calcium_wave_ode(t, y, params, Nx, Ny, dx, dy), tspan, y0, options);

% Reshape the solution
c = reshape(y(:, 1:Nx*Ny), [length(t), Nx, Ny]);
n = reshape(y(:, Nx*Ny+1:2*Nx*Ny), [length(t), Nx, Ny]);
P = reshape(y(:, 2*Nx*Ny+1:end), [length(t), Nx, Ny]);

% Visualize the results
figure;
for i = 1:10:length(t)
    imagesc(x, y, squeeze(c(i, :, :)));
    colorbar;
    title(sprintf('Calcium concentration at t = %.2f', t(i)));
    xlabel('x');
    ylabel('y');
    pause(0.1);
end


function dydt = calcium_wave_ode(t, y, params, Nx, Ny, dx, dy)
    % Reshape y into c, n, and P
    c = reshape(y(1:Nx*Ny), [Nx, Ny]);
    n = reshape(y(Nx*Ny+1:2*Nx*Ny), [Nx, Ny]);
    P = reshape(y(2*Nx*Ny+1:end), [Nx, Ny]);
    
    % Calculate Laplacians
    Lc = laplacian2D(c, dx, dy);
    LP = laplacian2D(P, dx, dy);
    
    % Calculate mu
    mu = params.mu0 + (params.mu1 * P) ./ (params.k_mu + P);
    
    % Calculate fluxes
    J_channel = params.k_flux * mu .* n .* (params.b + (params.V1 * c) ./ (params.k1 + c));
    J_pump = (params.gamma * c) ./ (params.k_gamma + c);
    J_leak = params.beta * ones(size(c));
    
    % Calculate derivatives
    dcdt = params.D_c * Lc + J_channel - J_pump + J_leak;
    dndt = (1 - c.^2 ./ (params.k2^2 + c.^2) - n) / params.tau_n;
    dPdt = params.D_p * LP - params.k_p * P;
    
    % Reshape derivatives into a column vector
    dydt = [dcdt(:); dndt(:); dPdt(:)];
end



function L = laplacian2D(u, dx, dy)
    [Nx, Ny] = size(u);
    L = zeros(Nx, Ny);
    
    % Interior points
    L(2:end-1, 2:end-1) = (u(1:end-2, 2:end-1) + u(3:end, 2:end-1) - 2*u(2:end-1, 2:end-1)) / dx^2 + ...
                          (u(2:end-1, 1:end-2) + u(2:end-1, 3:end) - 2*u(2:end-1, 2:end-1)) / dy^2;
    
    % Boundary conditions (no-flux)
    L(1, :) = L(2, :);
    L(end, :) = L(end-1, :);
    L(:, 1) = L(:, 2);
    L(:, end) = L(:, end-1);
end



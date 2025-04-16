% Model parameters
k_flux = 8.1;
b = 0.111;
V1 = 0.889;
k1 = 0.7;
gamma = 2.0;
k_gamma = 0.1;
beta = 0.01;
tau_n = 2.0;
k2 = 0.7;

% Time span for integration
tspan = [0 50];

% Different initial conditions
initial_conditions = [
    [0.1, 0.9];
    [0.5, 0.5];
    [0.9, 0.1]
];

% Different mu values
mu_values = [0.58, 0.7, 0.95];

% Create figure
figure('Position', [100, 100, 1200, 800]);

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    for j = 1:size(initial_conditions, 1)
        y0 = initial_conditions(j, :);
        
        % Solve ODE
        [t, y] = ode45(@(t, y) calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2), tspan, y0);
        
        % Plot results
        subplot(3, 2, 2*i-1);
        plot(t, y(:,1));
        hold on;
        title(['Calcium concentration (μ = ', num2str(mu), ')']);
        xlabel('Time');
        ylabel('[Ca^{2+}]');
        
        subplot(3, 2, 2*i);
        plot(t, y(:,2));
        hold on;
        title(['n (μ = ', num2str(mu), ')']);
        xlabel('Time');
        ylabel('n');
    end
end

% Add legend to the top right plot
subplot(3, 2, 2);
legend('Ca = 0.1, n = 0.9', 'Ca = 0.5, n = 0.5', 'Ca = 0.9, n = 0.1', 'Location', 'southeast');

function dydt = calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2)
    c = y(1);
    n = y(2);
    
    dcdt = k_flux * mu * n * (b + (V1 * c) / (k1 + c)) - (gamma * c) / (k_gamma + c) + beta;
    dndt = (n_inf(c, k2) - n) / tau_n;
    
    dydt = [dcdt; dndt];
end

function n_inf = n_inf(c, k2)
    n_inf = 1 - (c^2) / (k2^2 + c^2);
end



% Model parameters
k_flux = 8.1;
b = 0.111;
V1 = 0.889;
k1 = 0.7;
gamma = 2.0;
k_gamma = 0.1;
beta = 0.01;
tau_n = 2.0;
k2 = 0.7;

% Time span for integration
tspan = [0 50];

% Different initial conditions
initial_conditions = [
    0.1, 0.9;
    0.5, 0.5;
    0.9, 0.1
];

% Different mu values
mu_values = [0.58, 0.7, 0.95];

% Create figure
figure('Position', [100, 100, 1200, 800]);

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    subplot(1, length(mu_values), i); % Create a subplot for each mu value
    hold on; % Allow multiple trajectories on the same plot
    
    for j = 1:size(initial_conditions, 1)
        y0 = initial_conditions(j, :); % Initial condition
        
        % Solve ODE using ode45
        [t, y] = ode45(@(t, y) calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2), tspan, y0);
        
        % Plot phase plane trajectory
        plot(y(:,1), y(:,2), 'LineWidth', 1.5); % Plot [Ca^{2+}] vs n
    end
    
    % Add labels and title
    xlabel('[Ca^{2+}]');
    ylabel('n');
    title(['Phase Plane (\mu = ', num2str(mu), ')']);
    legend({'c = 0.1 n= 0.9', 'c = 0.5, n= 0.5', 'c = 0.9 n = 0.1'}, 'Location', 'best');
    grid on; % Add grid for better visualization
    
    hold off; % Stop adding to the current plot
end

function dydt = calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2)
    % Unpack variables
    c = y(1); % Calcium concentration
    n = y(2); % Gating variable
    
    % Define ODEs
    dcdt = k_flux * mu * n * (b + (V1 * c) / (k1 + c)) - (gamma * c) / (k_gamma + c) + beta;
    dndt = (n_inf(c, k2) - n) / tau_n;
    
    % Return derivatives as a column vector
    dydt = [dcdt; dndt];
end

function n_inf_val = n_inf(c, k2)
    % Steady-state function for gating variable n
    n_inf_val = 1 - (c^2) / (k2^2 + c^2);
end

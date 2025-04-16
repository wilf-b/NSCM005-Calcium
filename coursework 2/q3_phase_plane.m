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
initial_conditions = [0.1, 0.9; 0.5, 0.5; 0.9, 0.1];

% Different mu values
mu_values = [0.58, 0.7, 0.95];

%% Time Series Plots (Calcium & n vs Time)
figure('Name', 'Time Series Plots', 'Position', [100, 100, 1400, 800]);

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    % Plot Calcium Concentration vs Time
    subplot(2, length(mu_values), i);
    hold on;
    
    for j = 1:size(initial_conditions, 1)
        y0 = initial_conditions(j, :);
        [t, y] = ode45(@(t, y) calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2), tspan, y0);
        plot(t, y(:, 1), 'LineWidth', 1.5); % Calcium concentration
    end
    
    title(['[Ca^{2+}] vs Time (\mu = ', num2str(mu), ')']);
    xlabel('Time (s)');
    ylabel('[Ca^{2+}]');
    legend({'c = 0.1, n = 0.9', 'c = 0.5, n = 0.5', 'c = 0.9, n = 0.1'}, 'Location', 'best');
    grid on;
    
    % Plot n vs Time
    subplot(2, length(mu_values), length(mu_values) + i);
    hold on;
    
    for j = 1:size(initial_conditions, 1)
        y0 = initial_conditions(j, :);
        [t, y] = ode45(@(t, y) calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2), tspan, y0);
        plot(t, y(:, 2), 'LineWidth', 1.5); % n variable
    end
    
    title(['n vs Time (\mu = ', num2str(mu), ')']);
    xlabel('Time (s)');
    ylabel('n');
    grid on;
    
    hold off;
end

%% Phase Plane Plots (Calcium vs n)
figure('Name', 'Phase Plane Plots', 'Position', [150, 150, 1200, 600]);

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    subplot(1, length(mu_values), i);
    hold on;
    
    for j = 1:size(initial_conditions, 1)
        y0 = initial_conditions(j, :);
        [t, y] = ode45(@(t, y) calcium_model(t, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2), tspan, y0);
        plot(y(:, 1), y(:, 2), 'LineWidth', 1.5);
    end
    
    title(['Phase Plane (\mu = ', num2str(mu), ')']);
    xlabel('[Ca^{2+}]');
    ylabel('n');
    legend({'c = 0.1, n = 0.9', 'c = 0.5, n = 0.5', 'c = 0.9, n = 0.1'}, 'Location', 'best');
    grid on;
    
    hold off;
end

%% ODE System Function
function dydt = calcium_model(~, y, mu, k_flux, b, V1, k1, gamma, k_gamma, beta, tau_n, k2)
    c = y(1);
    n = y(2);
    
    dcdt = k_flux * mu * n * (b + (V1 * c) / (k1 + c)) - ...
           (gamma * c) / (k_gamma + c) + beta;
    dndt = (n_inf(c, k2) - n) / tau_n;
    
    dydt = [dcdt; dndt];
end

%% Steady-State Function for n
function n_inf_val = n_inf(c, k2)
    n_inf_val = 1 - (c^2) / (k2^2 + c^2);
end


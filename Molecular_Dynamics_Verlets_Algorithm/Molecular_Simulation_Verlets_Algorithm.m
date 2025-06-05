% -- % Molecular Simulation using Verlet's Algorithm % -- %

% This script simulates the van der Waals interaction between two Argon atoms,
% where one atom is fixed and the other has a mass of 40 Dalton.
% The Lennard-Jones potential is used to model the interaction.
% The Velocity Verlet algorithm is implemented to solve for position and velocity.

% The Velocity Verlet Algorithm:
%   r(t + h) = 2 * r(t) + h * v(t) + 0.5 * a(t) * h^2
%   v(t + h) = v(t) + 0.5 * h * (a(t) + a(t + h))
% where a = F/m and F is the negative derivative of the Lennard-Jones potential.

% --- Main Simulation --- %

% Simulate with different time steps as instructed

% Time step = 1 ps (N = 3)
Verlet(3);

% Time step = 0.3 ps (N = 10)
Verlet(10);

% Time step = 0.15 ps (N = 20)
Verlet(20);

% --- Velocity Verlet Simulation Function --- %
function Verlet(N)
    % Verlet(N)
    % Simulates the motion of an Argon atom using the Velocity Verlet algorithm.
    % N: Number of time steps (the total simulated time is 3 ps)

    % --- Parameters and Initial Conditions --- %
    r0 = 4;                  % Initial position [Å]
    v0 = 4.9;                % Initial velocity [Å/ps]
    eps = 0.24 * 418.4;      % Depth of potential well [J/mol] (converted to appropriate units)
    sigma = 3.4;             % Distance at which potential is zero [Å]
    m = 40;                  % Mass [Dalton]
    h = 3 / N;               % Time step [ps] (total time = 3 ps)

    % Preallocate arrays for position, velocity, energy, and time
    R = zeros(1, N + 1);     % Position array
    V = zeros(1, N + 1);     % Velocity array
    E = zeros(1, N + 1);     % Energy array
    Z = zeros(1, N + 1);     % Time array

    % Set initial values
    R(1) = r0;
    V(1) = v0;

    % --- Velocity Verlet Algorithm --- %
    for i = 1:N
        % Define the negative derivative of the Lennard-Jones potential (force)
        neg_derivative = @(r) ((-48 * eps * sigma^12) / (r^13) + (24 * eps * sigma^6) / (r^7));

        % Calculate position at next time step
        R(i + 1) = 2 * R(i) + h * V(i) + (h^2 / (2 * m)) * neg_derivative(R(i));

        % Calculate velocity at next time step
        V(i + 1) = V(i) + (h / (2 * m)) * (neg_derivative(R(i)) + neg_derivative(R(i + 1)));

        % Calculate energy at current time step
        E(i) = 0.5 * m * V(i)^2 - neg_derivative(R(i));

        % Update time
        Z(i + 1) = Z(i) + h;
    end

    % Final energy calculation
    E(N + 1) = 0.5 * m * V(N + 1)^2 - neg_derivative(R(N + 1));

    % --- Plotting Results --- %
    % Plot position
    subplot(1, 3, 1)
    plot(Z, R, '-o')
    title('Position')
    xlabel('Time [ps]')
    ylabel('Position [Å]')
    legend('h = 1', 'h = 0.3', 'h = 0.15', 'Location', 'northwest')
    hold on;

    % Plot velocity
    subplot(1, 3, 2)
    plot(Z, V, '-o')
    title('Velocity')
    xlabel('Time [ps]')
    ylabel('Velocity [Å/ps]')
    legend('h = 1', 'h = 0.3', 'h = 0.15', 'Location', 'northwest')
    hold on;

    % Plot energy
    subplot(1, 3, 3)
    plot(Z, E, '-o')
    title('Energy')
    xlabel('Time [ps]')
    ylabel('Energy')
    legend('h = 1', 'h = 0.3', 'h = 0.15', 'Location', 'northwest')
    hold on;
end






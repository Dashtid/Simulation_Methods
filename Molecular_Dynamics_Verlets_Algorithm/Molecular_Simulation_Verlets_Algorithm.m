% -- % Answer to Quiz: Molecular Sim.  % --- %
% ------------------------------------------ %
% -------- Written by David Dashti --------- %
% ------------------------------------------ %

% The task is to simulate, using Verlet's Method, the van der Waal
% interactions of two Argon atoms. One of the atoms is fixed and the
% other one has a mass of 40 Dalton. We model the interactions with the
% Lennart-Jones equation of potential and we are solving for both position
% & velocity of the moving object, being the atoms. To solve this I am
% going to implement the Velocity Verlet Algorithm

% The Velocity Verlet's Algorithm is the following:
% r(t + h) = 2 * r(t) + h * v(t) + 1/2 * a(t) * h^2 
% v(t + h) = v(t) + h/2 * (a(t) + a(t + h)

% Also:
% a = F/m
% F = ((-4) * eps) * (((-12) * ((sigma^12)/(r(t)^13))) - ((-6) * ...
% *((sigma^6)/(r(t)^7))))
% F is the negative derivation of the Lennart-Jones potential 

% --- Main --- %

% Solving the system with different time steps as instructed %

% Time step = 1 ps
Verlet(3);

% Time step = 0.3 ps
Verlet(10);

% Time step = 0.15 ps
Verlet(20);


% --- Function --- % 
function Verlet(N)

% --- SETUP --- %
% Values given in the task %
r0 = 4;
v0 = 4.9;
eps = 0.24 * 418.4;
sigma = 3.4;
m = 40;
h = 3/N;

% Creating lists needed for simulation
R = zeros(1,N+1);
V = zeros(1,N+1);
E = zeros(1,N+1);
Z = zeros(1,N+1);

% Setting initial values
R(1) = r0;
V(1) = v0;

% --- ALGORITHM --- %
% Running the Algorithm
for i = 1:N

% Calculating the position
R(i + 1) = 2 * R(i) + h * V(i) + ...
+ (0.5 * h^2 * (-4) * eps) * (((-12) * ((sigma^12)/(R(i)^13))) - ((-6) * ((sigma^6)/(R(i)^7))));

% Calculating the velocity
V(i + 1) = V(i) + ...
+ (h/2 * 1/m * (-4) * eps) * (((-12) * ((sigma^12)/(R(i)^13))) - ((-6) * ((sigma^6)/(R(i)^7)))) + ...
+ (h/2 * 1/m * (-4) * eps) * (((-12) * ((sigma^12)/(R(i + 1)^13))) - ((-6) * ((sigma^6)/(R(i + 1)^7))));


% Calculating the energy
E(i) = 1/2 * m * V(i)^2 + ...
+ (4 * eps) * (((-12) * ((sigma^12)/(R(i)^13))) - ((-6) * ((sigma^6)/(R(i)^7))));    


%Updating time step
Z(i + 1) = Z(i) + h;
end

% Final energy calculation
E(N+1) = 1/2 * m * V(N+1)^2 + 4 * eps * ((sigma/R(N + 1))^12 - (sigma/R(N + 1))^6);

% --- PLOTTING --- %
% Plotting the positions
subplot(1,3,1)
plot(Z',R',"-o")
title('Position')
legend('h = 1','h = 0.3','h = 0.15','Location','northwest')
hold on;

% Plotting the velocities
subplot(1,3,2)
plot(Z,V,"-o")
title('Velocity')
legend('h = 1','h = 0.3','h = 0.15','Location','northwest')
hold on;

% PLotting the energy
subplot(1,3,3)
plot(Z,E,"-o")
title('Energy')
legend('h = 1','h = 0.3','h = 0.15','Location','northwest')
hold on;
end




    

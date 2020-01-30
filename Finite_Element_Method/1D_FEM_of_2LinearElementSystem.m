% ----- % Answer to Quiz Question 2B  % ---- %
% ------------------------------------------ %
% -------- Written by David Dashti --------- %
% ------------------------------------------ %

% This is a smaller problem which was a part of a quiz
% Normally, A = 1 but in this case it was A = x + 1 which made for a ...
% different stiffness matrix

% % Declaring the variables given in the task
% I also going to use a lot of formulas implicitly
E = 1; 
L = 1; 
h = L/2;  % Changed to 2 instead of 4 to accound for change in linear elements
F_L = 1;

% % Stiffness Matrixes
% Calculated by hand to verify
k = [E*(1/2 + 1/h), -E*(1/2 + 1/h), 0;
    -E*(1/2 + 1/h), E*(1 + 1/h), -E*(1/2 + 1/h);
    0, -E*(1/2 + 1/h), E*(1/2 + 1/h)];

% % Solving for the partial displacements
k_partial = [E*(1 + 1/h), -E*(1/2 + 1/h);
            -E*(1/2 + 1/h), E*(1/2 + 1/h)];
 
F_partial = [h;
             h/2 + F_L];

u_partial = k_partial\F_partial;

% % Final calculation
% Solving for force that is to be subtracted
F_0 = h/2 + E*(1/2 + 1/h)*u_partial(1);

% Creating a forces vector
F = [h/2 - F_0;
    F_partial];

u_final = k\F;

h_v = [0*h, 1*h, 2*h];

% % Plotting of the solution
plot(h_v, u_final);
xlabel('h');
ylabel('u');
title('Solution to the question 2 using 2 linear elements');
% ------ % Solution to Laboration 3 % ------ %
% ------- % Finite Element Problem % ------- %
% ------------------------------------------ %
% Written by David Dashti & Filip Söderquist %
% ------------------------------------------ %

% Given variables
I = (pi * 0.05^4)/4;            % Inertia
E = 100 * 10^6;                 % Young's modulus
L = 2/4;                        % Lenght of element in meters
A = pi*0.05^2;                  % Area of cross-section

% Angles given
beta = 70 * pi/180;             % Angle of 1st element
lambda = 55 * pi/180;           % Angle of 2nd element
delta = 30 * pi/180;            % Angle of 3rd element

%% Solving the task using Finite Element Method

% % Transformation matrixes

% Transformation matrix for 1st element
T_beta = [cos(beta), sin(beta), 0, 0, 0, 0;
    -sin(beta), cos(beta), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(beta), sin(beta), 0;
    0, 0, 0, -sin(beta), cos(beta), 0;
   0, 0, 0, 0, 0, 1];

% Transformation matrix for 2nd element
T_lamda = [cos(lambda), sin(lambda), 0, 0, 0, 0;
    -sin(lambda), cos(lambda), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(lambda), sin(lambda), 0;
    0, 0, 0, -sin(lambda), cos(lambda), 0;
   0, 0, 0, 0, 0, 1];    

% Transformation matrix for 3rd element
T_delta = [cos(delta), sin(delta), 0, 0, 0, 0;
    -sin(delta), cos(delta), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(delta), sin(delta), 0;
    0, 0, 0, -sin(delta), cos(delta), 0;
   0, 0, 0, 0, 0, 1];

% % Transformation matrix for 4th element
T_sista = [cos(0), sin(0), 0, 0, 0, 0;
    -sin(0), cos(0), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(0), sin(0), 0;
    0, 0, 0, -sin(0), cos(0), 0;
   0, 0, 0, 0, 0, 1];

% Transposed matrices 
T_beta_transp = T_beta';
T_lamda_transp = T_lamda';
T_delta_transp = T_delta';
T_sista_transp = T_sista';

% % Stiffness matrix
k_e = [E*A/L, 0, 0, -E*A/L, 0, 0;
    0, 12*E*I/L^3, 6*E*I/L^2, 0, -12*E*I/L^3, 6*E*I/L^2;
    0, 6*E*I/L^2, 4*E*I/L, 0, -6*E*I/L^2, 2*E*I/L;
    -E*A/L, 0 ,0 ,E*A/L 0, 0;
    0, -12*E*I/L^3, -6*E*I/L^2, 0, 12*E*I/L^3, -6*E*I/L^2;
    0, 6*E*I/L^2, 2*E*I/L, 0, -6*E*I/L^2, 4*E*I/L];

% % Creating the global K matrix

% Calculating the Ke matrices for each element
k1 = T_beta_transp * k_e * T_beta;
k2 = T_lamda_transp * k_e * T_lamda;
k3 = T_delta_transp * k_e * T_delta;
k4 = T_sista_transp * k_e * T_sista;

% Initiating the global K matrix
A = zeros(15,15);
B = zeros(15,15);
C = zeros(15,15);
D = zeros(15,15);

% Filling the values into the global K matrix

% Values from k1
A(1:3,1:3) = k1(1:3,1:3);
A(1:3,4:6) = k1(1:3,4:6);
A(4:6,1:3) = k1(4:6,1:3);
A(4:6,4:6) = k1(4:6,4:6);

% Values from k2
B(4:6,4:6) = k2(1:3,1:3);
B(4:6,7:9) = k2(1:3,4:6);
B(7:9,4:6) = k2(4:6,1:3);
B(7:9,7:9) = k2(4:6,4:6);

% Values from k3
C(7:9,7:9) = k3(1:3,1:3);
C(7:9,10:12) = k3(1:3,4:6);
C(10:12,7:9) = k3(4:6,1:3);
C(10:12,10:12) = k3(4:6,4:6);

% Values from k4
D(10:12,10:12) = k4(1:3,1:3);
D(10:12,13:15) = k4(1:3,4:6);
D(13:15,10:12) = k4(4:6,1:3);
D(13:15,13:15) = k4(4:6,4:6);

% Summing these matrices together to create a final matrix
K_global = A + B + C + D;

% Erase first three rows and columns 
K_global_1212 = K_global((4:15),(4:15));

% Creating a global force vector F
F_global_121 = zeros(12, 1);
F_global_121(11, 1) = -20;

% Calculating the displacements
D_global_121 = K_global_1212\F_global_121;

D_global = [zeros(3,1);
            D_global_121];

D_zero = zeros(15,1);

%Add displacements and plot
hold on;
rodPlot(beta, lambda, delta, L, D_global);
rodPlot(beta, lambda, delta, L, D_zero);
title('The fishing rod in the original position and then bent at the tip');
xlabel('x [m]');
ylabel('y [m]');


% % This is a function that plots out the rod
% Given that you know the angles and the displacements
function rodPlot(ang1, ang2, ang3, L, disp)
    pos = [0 + disp(1), 0 + disp(2)];
    pos = [pos; 
          L*cos(ang1) + disp(4), L*sin(ang1) + disp(5)];
    pos = [pos; 
          pos(2,1) + L*cos(ang2) + disp(7), pos(2,2) + L*sin(ang2) + disp(8)];
    pos = [pos;
          pos(3,1) + L*cos(ang3) + disp(10), pos(3,2) + L*sin(ang3) + disp(11)];
    pos = [pos;
          pos(4,1) + L + disp(13), pos(4,2) + disp(14)];
    plot(pos(:, 1),pos(:, 2), '-o')
end
% ------ % Solution to Laboration 2 % ------ %
% ------- % Boundary Value Problem % ------- %
% ------------------------------------------ %
% Written by David Dashti & Filip Söderquist %
% ------------------------------------------ %


Nvec = [9, 19, 39, 79];      % Number of steps                      
Vvec = [0, 0,1, 0,5, 1] ;    % Used when looing through v

%% Main

% Looping through all possible step sizes (N)
figure(1);
plotNvec(Nvec);

% Looping through all possible v-values
figure(2);
plotVvec(Vvec);


% Functions for solving and plotting
% Plots with N = 79 & different V-values
function plotVvec(Vvec)
for I = 1:length(Vvec)

% Defining constants given
L = 10;
a = 1;
b = 3;
Q0 = 50;
k = 0.5;
p = 1;
C = 1;
Tout = 300;  % End temperature
T0 = 400;    % Starting temperature
z0 = 0;

v = Vvec(I);                                 % Setting v value   
N = 79;                                      % Hardcoding 79 as N-value
boundaries = [0, L];                         % Setting the boundaries
h = (boundaries(end)-boundaries(1))/(N+1);   % Calculating step size

% All steps Zi with both the inner points and the boundary points x=x0, x=xN+1
z = boundaries(1):h:boundaries(end); 

% % Defining the tridiagonal matrix A % % 
% Based on the Central Difference Formula % 

% Creating the matrices
A = zeros(N, N);
Q = zeros(1,N);

% Values for the diagonal matrix A
m11 = 2*k;                           % Diagonal Elements
m12 = (h*v*p*C)/2 - k;               % Elements in the col. to the right of the diagnoal
m21 = -(h*v*p*C)/2 - k;              % Elements in the row. below the diagonal


% Insertion of values into the matrix
for i = 1:N-1
    A(i,i) = m11;
    A(i,i+1) = m12;
    A(i+1,i) = m21;
end

% This is the last value which need to be inserted at spot A(end,end)
A(N,N) = m11; 

% % Defining matrix Q: Vector of length N with 3 possible states % %

% Running the formula for Q in the open intervall a to b
for i = 1:N

% From 0 to a
if (z(i + 1) >= 0 && z(i + 1) < a)
Q(i) = 0;   

% Values between a and b
elseif (z(i + 1) >= a && z(i + 1) <= b) 

Q(i) = Q0 * sin(((z(i + 1) - a)*pi)/(b - a)); 

% From: b to L
elseif (z(i + 1) > b && z(i + 1) <= L)
    
Q(i) = 0;

end

end

% % Defining the boundary vector % %

% This is a vector that has a length of N and only ... 
% consists of the boundary condition values that are given ... 
% which is the first and last value. The rest are 0.

% Creating vector
BCvec = zeros(1,N);

% Boundary Conditions
BCvec(1) = T0 * (-(h*v*p*C/2) - k) ;
BCvec(end) = Tout * ((h*v*p*C/2) - k);
    
% % Obtaining the Au = b form % %

B = (h^2)*Q - BCvec; 

% The next step is to calculate the vector U.
% U is the vector with all the calculated temperatures in-between our ...
% boundary condition values.

T_approx = A\B';

% To obtain the vector T, we have to add the initial values to the vector U

T_approx = [T0;
    T_approx;
    Tout];

% %  Plotting the solution: Temperatures at points z

plot(z', T_approx);
ylabel('Approximated temperature (T)');
xlabel('Length (z)')
hold on

end
    
hold off
title('Solution to BVP Heat Conduction w. different v')

end

% Plots with v = 0 & different N-values
function plotNvec(Nvec)
for I = 1:length(Nvec)
    
% Defining constants given

L = 10;
a = 1;
b = 3;
Q0 = 50;
k = 0.5;
p = 1;
C = 1;
Tout = 300;  % End temperature
T0 = 400;    % Starting temperature
z0 = 0;



v = 0;                                       % Hardcoding V-value
N = Nvec(I);                                 % Choosing one N-value
boundaries = [0, L];                         % Setting the boundaries
h = (boundaries(end)-boundaries(1))/(N+1);   % Calculating step size

% All steps Zi with both the inner points and the boundary points x=x0, x=xN+1
z = boundaries(1):h:boundaries(end); 

% % Defining the tridiagonal matrix A % % 
% Based on the Central Difference Formula % 

% Creating the matrices
A = zeros(N, N);
Q = zeros(1,N);

% Values for the diagonal matrix A
m11 = 2*k;                           % Diagonal Elements
m12 = (h*v*p*C)/2 - k;               % Elements in the col. to the right of the diagnoal
m21 = -(h*v*p*C)/2 - k;              % Elements in the row. below the diagonal


% Insertion of values into the matrix
for i = 1:N-1
    A(i,i) = m11;
    A(i,i+1) = m12;
    A(i+1,i) = m21;
end

% This is the last value which need to be inserted at spot A(end,end)
A(N,N) = m11; 

% % Defining matrix Q: Vector of length N with 3 possible states % %

% Running the formula for Q in the open intervall a to b
for i = 1:N

% From 0 to a
if (z(i + 1) >= 0 && z(i + 1) < a)
Q(i) = 0;   

% Values between a and b
elseif (z(i + 1) >= a && z(i + 1) <= b) 

Q(i) = Q0 * sin(((z(i + 1) - a)*pi)/(b - a)); 

% From: b to L
elseif (z(i + 1) > b && z(i + 1) <= L)
    
Q(i) = 0;

end

end

% % Defining the boundary vector % %

% This is a vector that has a length of N and only ... 
% consists of the boundary condition values that are given ... 
% which is the first and last value. The rest are 0.

% Creating vector
BCvec = zeros(1,N);

% Boundary Conditions
BCvec(1) = T0 * (-(h*v*p*C/2) - k) ;
BCvec(end) = Tout * ((h*v*p*C/2) - k);
    
% % Obtaining the Au = b form % %

B = (h^2)*Q - BCvec; 

% The next step is to calculate the vector U.
% U is the vector with all the calculated temperatures in-between our ...
% boundary condition values.

T_approx = A\B';

% To obtain the vector T, we have to add the initial values to the vector U

T_approx = [T0;
    T_approx;
    Tout];

% %  Plotting the solution: Temperatures at points z

plot(z', T_approx);
ylabel('Approximated temperature (T)');
xlabel('Length (z)')
hold on

end
    
hold off
title('Solution to BVP Heat Conduction w. different N')

end
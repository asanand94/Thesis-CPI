function [ x, y ] = forward_euler ( f_ode, xRange, yInitial, numSteps )
% [ x, y ] = forward_euler ( f_ode, xRange, yInitial, numSteps ) uses
% Euler?s explicit method to solve a system of first-order ODEs
% dy/dx=f_ode(x,y).
% f = function handle for a function with signature
% fValue = f_ode(x,y)
% where fValue is a column vector
% xRange = [x1,x2] where the solution is sought on x1<=x<=x2
% yInitial = column vector of initial values for y at x1
% numSteps = number of equally-sized steps to take from x1 to x2
% x = row vector of values of x
% y = matrix whose k-th column is the approximate solution at x(k).
x(1) = xRange(1);
h = ( xRange(2) - xRange(1) ) / numSteps;
y(:,1) = yInitial;
for k = 1 : numSteps
x(1,k+1) = x(1,k) + h;
y(:,k+1) = y(:,k) + h * f_ode( x(k), y(:,k) );
end
%1D diffusion equation calculated discretely
%assumption. Left boundary = dirichlet. right = neumann
%initial conditions
close all
T = 12000;      %total time
m = 600;       %time steps
dt = T/m;      %delta t
a = 1;          %total space
n = 39;         %spatial steps; (n+2) grid points
dx = a/(n+1);   %delta x
D = .00001;         %diffusion coefficient
g1 = 1;         %left dirichlet boundary
g2 = 2;         %right neumann boundary

s = D*dt/(dx^2);    %calculate s

%create vector for heat for each point in space for each time step
uVector = zeros(n,m);

%create A vector for multiplication
A = zeros(n);
    for i = 1:n   
        for j = 1:n

            if i == j && i~=n
                A(i,j) = 1-2*s;
            elseif abs(i-j) == 1
                A(i,j) = s;
            elseif i == j && i == n
                A(i,j) = 1-s;
            end
            
        end
    end
    
%create b vector for addition
b = zeros(n,1);
b(1,1) = s*g1;
b(n,1) = s*dx*g2;
%{
%create and fill in initial column in vector
startVector = zeros(n,1);

for k = 1:n
    %use initial condition to fill in first column of uVector
    startVector(k,1) =.2; %2*dx*k + sin(2*pi*dx*k) + 1
end
%}
startVector = rand(n,1);
%calculate uValues for for each time step at each place
for l = 1:m
        uVector(:,l) = A*startVector + b;
        startVector = uVector(:,l); %refresh start vector
       
end
figure();
color = @(t) [t / m , 0, 0];
videoWriter = VideoWriter('diffusion.avi');
open(videoWriter);
for t=1:m
    plot(linspace(0, 1, n), uVector(:, t), 'Color', color(t));
    ylim([.3,max(uVector(:))]); xlabel('x'); ylabel('u(t; x)');
    writeVideo(videoWriter, getframe(gcf()));
end
%{
for t=round(linspace(1,m,32)) 
    plot(linspace(0, 1, n), uVector(:, t), 'Color', color(t));
    ylim([.3,2]);
    writeVideo(videoWriter, getframe(gcf()));
end
%}
close(videoWriter);

xlabel('x'); ylabel('u(t; x)');
%{
figure()
plot(linspace(0, 1, n), uVector(:,end));
xlabel('x'); ylabel('u(t; x)');
%}
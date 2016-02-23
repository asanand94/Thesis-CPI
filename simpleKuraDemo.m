% file simpleKuraDemo.m
function simpleKuraDemo

n = 64;

% Make a random "network".
A = randi(2, n) - 1;
p = .5;
A(A>p) = 1;
A(A<=p) = 0;

% Make the network symmetric (undirected).
for i=1:n
    for j=1:n
        A(i,j) = A(j,i);
    end
end

% Make a uniform natural frequency distribution and set the coupling
% strength.
om = rand(n, 1) * 1.5; om=om-mean(om);
K = 1.95;
function out=dThetaDt(~, thetaVec)
    out = zeros(n, 1);
    for i=1:n
        sumsini = 0;
        for j=1:n
            sumsini = sumsini + A(i,j) * sin(thetaVec(j) - thetaVec(i));
        end
        out(i) = om(i) + K/n * sumsini;
    end
end

theta0 = rand(n,1) * 3.1415;
tmax = 40;
[times, history] = ode45(@dThetaDt, 0:.20:tmax, theta0);
%{
% Plot the trajectory.
 figure();
 plot(times, history, 'k');
 xlabel('time');
 ylabel('\theta_i(t)');
%}
% Plot a couple low-dimensional states vs omegas.
%{
tis = [...
        round(.1*numel(times)), ...
        round(.4*numel(times)), ...
        round(.7*numel(times)), ...
        round(.9*numel(times))  ...
       ];
%}
figure('Position', [0, 0, 1600, 900]);
tis = 1:numel(times);
videoWriter = VideoWriter('kuramoto.avi');
open(videoWriter);
%theta vs omega plots for varying times
for i=1:numel(tis)
    ti = tis(i);
    %subplot(1, numel(tis), i);
    scatter(om, history(ti, :));
    xlabel('\omega_i');
    ylabel(sprintf('\\theta_i(t=%f)', times(ti)));
    thetaMin = min(min(history));
    thetaMax = max(max(history));
    ylim([thetaMin, thetaMax]);
    
    writeVideo(videoWriter, getframe(gcf()));
end
close(videoWriter);
end
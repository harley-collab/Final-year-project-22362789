function ring_sim_4_Regimes()

tic; % start timer
basePlotDir = 'plots_4_1_1';
if ~exist(basePlotDir, 'dir')
    mkdir(basePlotDir);
end

m = 1.0;      % mass of the ring
R = 1.025;    % radius of the ring
h = 0.88;     % ring height
w = 0.05;     % difference between inner and outer diameter
g = 9.81;     % gravity
I1 = 0.5772;  % inertia along e1
I3 = 1.0253;  % inertia along e3

r_Gc = [0; -h/2; R];
theta0     = 0.55;
phi0       = 0;
psi0       = 0;
theta_dot0 = 0;
phi_dot0   = 10;
psi_dot0   = 0;

tspan = [0 600];  

options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-2, ...
                 'Events', @(t,y) combinedEvents(t,y));


C1_values = 0.010;%drag
C2_values = 0.01;
C3_values =  0.05;
for C1 = C1_values
for C2 = C2_values
for C3 = C3_values

w10 = 0;
w20 = psi_dot0 + phi_dot0*sin(theta0);
w30 = phi_dot0*cos(theta0);


y0 = [ 0; 0; R*cos(theta0); w10; w20; w30; theta0; phi0; psi0 ];

[t, y] = ode45(@(t,y) eom(t,y, m, R, h, g, I1, I3, C1, C2, C3, r_Gc), ...
               tspan, y0, options);

Ff  = zeros(length(t),1);
rcdd = zeros(length(t),3);

for k = 1:length(t)
    [Ff(k), rcdd(k,:)] = turntrack(y(k,:),m,R,h,g,I1,I3,C1,C2,C3,r_Gc);
end

idx = find(Ff(1:end-1).*Ff(2:end) < 0, 1);

if ~isempty(idx)
    if idx > 1 && idx+1 <= length(t)
        sharpness = (Ff(idx+1) - Ff(idx-1)) / (t(idx+1) - t(idx-1));%eqn.5.4
        fprintf('C1=%.3f C2=%.3f C3=%.8f | Sharpness dFf/dt|t* = %.6f\n', C1, C2, C3, sharpness);
    else
        sharpness = NaN;
        fprintf('C1=%.3f C2=%.3f C3=%.8f | Sharpness: bracket at boundary, cannot compute\n', C1, C2, C3);
    end
else
    sharpness = NaN;
    fprintf('C1=%.3f C2=%.3f C3=%.8f | No transition found\n', C1, C2, C3);
end



tag = sprintf('C1_%g_C2_%g_C3_%g', C1, C2, C3);

fs = 16;  
ft = 18;   
lw = 2.5;  
 
fig = figure('visible','off');
plot(y(:,1), y(:,2), 'LineWidth', lw);
axis equal; grid on;
xlabel('x (m)', 'FontSize', fs);
ylabel('y (m)', 'FontSize', fs);
title('Trajectory in XY Plane', 'FontSize', ft);
set(gca, 'FontSize', fs);
print(fig, fullfile(basePlotDir, ['XY_' tag '.png']), '-dpng', '-r150');
close(fig);
 
fig = figure('visible','off');
plot(t, Ff, 'LineWidth', lw);
xlabel('t (s)',    'FontSize', fs);
ylabel('F_f / m', 'FontSize', fs);
title('Friction Force vs. Time', 'FontSize', ft);
grid on; box on;
set(gca, 'FontSize', fs);
print(fig, fullfile(basePlotDir, ['Ff_' tag '.png']), '-dpng', '-r150');
close(fig);
 
fig = figure('visible','off');
plot3(y(:,1), y(:,2), y(:,3), 'b', 'LineWidth', lw);
xlabel('x (m)', 'FontSize', fs);
ylabel('y (m)', 'FontSize', fs);
zlabel('z (m)', 'FontSize', fs);
title('r_A Trajectory', 'FontSize', ft);
grid on; box on;
set(gca, 'FontSize', fs);
print(fig, fullfile(basePlotDir, ['3d_' tag '.png']), '-dpng', '-r150');
close(fig);
 
fig = figure('visible','off');
plot3(y(:,4), y(:,5), y(:,6), 'LineWidth', lw);
xlabel('\omega_1 (rad/s)', 'FontSize', fs);
ylabel('\omega_2 (rad/s)', 'FontSize', fs);
zlabel('\omega_3 (rad/s)', 'FontSize', fs);
title('Angular Velocity Phase Portrait', 'FontSize', ft);
grid on; box on;
set(gca, 'FontSize', fs);
print(fig, fullfile(basePlotDir, ['w3_' tag '.png']), '-dpng', '-r150');
close(fig);
 
theta = y(:,7);
phi   = y(:,8);
psi   = y(:,9);
 
fig = figure('visible','off');
subplot(3,1,1);
plot(t, theta, 'LineWidth', lw);
ylabel('\theta (rad)', 'FontSize', fs);
set(gca, 'FontSize', fs); grid on;
subplot(3,1,2);
plot(t, phi, 'LineWidth', lw);
ylabel('\phi (rad)', 'FontSize', fs);
set(gca, 'FontSize', fs); grid on;
subplot(3,1,3);
plot(t, psi, 'LineWidth', lw);
ylabel('\psi (rad)', 'FontSize', fs);
set(gca, 'FontSize', fs); grid on;
xlabel('Time (s)', 'FontSize', fs);
sgtitle('Euler Angles vs. Time', 'FontSize', ft);
print(fig, fullfile(basePlotDir, ['ang_' tag '.png']), '-dpng', '-r150');
close(fig);
end
end
end
fprintf('\nCompleted in %.2f seconds\n', toc);
end
% Ornithopter flapping simulation
clc;
clear all;
close all;

%% Simulation parameters
dt = 0.01;                          % Time step (s)
t_end = 10;                         % End time (s)
t = 0:dt:t_end;                     % Time array

%% Ornithopter parameters
m = 0.003;                          % Mass (kg)
wing_span = 0.1;                    % Wingspan (m)
wing_area = 0.01;                   % Wing area (m^2)
cg = [0, 0, 0]';                    % Center of gravity (m)
I = diag([1, 1, 1]);                % Moment of inertia (kg*m^2)

%% Initial conditions
x0 = [0, 0, 0]';                    % Initial position (m)
v0 = [0, 0, 0]';                    % Initial velocity (m/s)
q0 = [1, 0, 0, 0]';                 % Initial orientation (quaternion)
w0 = [0, 0, 0]';                    % Initial angular velocity (rad/s)

%% Flapping wing parameters
freq = 25;                          % Flapping frequency (Hz)
amp = 30;                           % Flapping amplitude (deg)
phase = 0;                          % Flapping phase (deg)
Cl_max = 1.5;                       % Maximum lift coefficient
Cd = 0.1;                           % Drag coefficient
rho = 1.225;                        % Air density (kg/m^3)
wing_vel = [0, 0, 0]';              % Wing velocity (m/s)
wing_pos = [0, 0, 0]';              % Wing position (m)

%% Run simulation
x = x0;
v = v0;
q = q0;
w = w0;

figure;
hold on;
axis equal;
grid on;
view(3);

for i = 1:length(t)
    % Flapping motion
    wing_pitch = amp*sin(2*pi*freq*t(i) + phase);
    wing_vel = [0, 0, -2*pi*freq*amp*cos(2*pi*freq*t(i) + phase)]';
    wing_pos = [0, wing_span/2*sin(deg2rad(wing_pitch)), -wing_span/2*cos(deg2rad(wing_pitch))]';

    % Calculate aerodynamic forces and torques on the body
    [f_aero, tau_aero] = calc_aero_forces_torques(wing_vel, wing_pitch, Cl_max, Cd, rho, wing_area, cg, I);

    % Calculate acceleration and angular acceleration of the body
    a = f_aero/m;
    alpha = I\(tau_aero - cross(w, I*w));

    % Update position, velocity, orientation, and angular velocity
    [x, v, q, w] = update_states(x, v, q, w, a, alpha, dt);

    % Plot results
    plot3(x(1), x(2), x(3), 'bo');
    plot_quat(q, x, 0.05);
    pause(0.01);
end

%% Functions
% Calculate aerodynamic forces and torques on the body
function [f_aero, tau_aero] = calc_aero_forces_torques(wing_vel, wing_pitch, Cl_max, Cd, rho, wing_area, cg, I)
    % Rotation matrix from wing to body
R_wb = [1, 0, 0; 0, cos(deg2rad(wing_pitch)), -sin(deg2rad(wing_pitch)); 0, sin(deg2rad(wing_pitch)), cos(deg2rad(wing_pitch))];
% Velocity of the air relative to the wing
v_rel = wing_vel - cross([0, 0, 2*pi*freq*wing_pos(2)]', R_wb*[0, wing_span/2, 0]');

% Lift coefficient
alpha = atan2(v_rel(3), v_rel(1));
Cl = Cl_max*sin(2*alpha);

% Aerodynamic force
f_aero = 0.5*rho*norm(v_rel)^2*wing_area*Cl*R_wb'*[0, 1, 0]';

% Drag force
f_drag = -0.5*rho*norm(v_rel)^2*wing_area*Cd*R_wb'*[1, 0, 0]';

% Aerodynamic torque
tau_aero = cross(wing_pos - cg, f_aero + f_drag);
end

% Update position, velocity, orientation, and angular velocity
function [x, v, q, w] = update_states(x, v, q, w, a, alpha, dt)
% Update position and velocity
x = x + vdt + 0.5*adt^2;
v = v + adt;
% Update orientation and angular velocity
[q, w] = update_quat(q, w, alpha, dt);
end

% Update quaternion and angular velocity
function [q, w] = update_quat(q, w, alpha, dt)
q_dot = 0.5*[0, -w'; w, skew_sym(w)]*q;
q = q + q_dotdt;
w_dot = alpha - skew_sym(w)*I*w;
w = w + w_dot*dt;
end

% Plot quaternion as a vector with an arrow
function plot_quat(q, pos, scale)
r = quat2rotm(q);
v = r*[1, 0, 0]';
quiver3(pos(1), pos(2), pos(3), v(1), v(2), v(3), scale);
end

% Skew-symmetric matrix of a vector
function S = skew_sym(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

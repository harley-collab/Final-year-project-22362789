function [Ff, rcdd] = turntrack(y, m, R, h, g, I1, I3, C1, C2, C3, r_Gc)
r_cx = y(1); % position
r_cy = y(2); 
r_cz = y(3);
w1 = y(4); %omega
w2 = y(5); 
w3 = y(6);
theta = y(7); %angle
phi = y(8); 
psi = y(9);

phi_dot   = w3/cos(theta);%angular acceleration
theta_dot = w1;
psi_dot   = w2 - phi_dot*sin(theta);

omega = [w1; w2; w3];%eqn. 2.2

r_A = [0; -h/2; 0];% position 
Omega = omega - [0; psi_dot; 0];%orbital acceleration 


v_Cc = cross(omega, r_Gc);
v_Ce = Rc2e * v_Cc;
v_A  = v_Ce - Rc2e * cross(Omega, r_A); %eqn. 2.4

I = diag([I1 I1 I3]); %inertia tensor 
j = I + m*[ ((h^2)/4)+R^2, 0, 0;%augmented inertia tensor eqn. 2.12
            0, R^2, -h*R/2;
            0, -h*R/2, h^2/4 ];

Lg = I * omega;%momentum


tau_drag = -[C1*abs(w1)*w1; C2*abs(w2)*w2; C3*abs(w3)*w3]; %eqn. 2.18

f = -cross(Omega,Lg) - m*cross(r_Gc, cross(Omega,v_Cc)) + m*g*[R*sin(theta)-(h/2)*cos(theta);0;0];
omega_dot = j \ f + tau_drag;

es = [0; cos(theta); -sin(theta)];
rhs = -I*omega_dot - cross(Omega, Lg);
F = cross(r_Gc, rhs) / dot(r_Gc,r_Gc);%eqn. 2.14

rcdd = F/m - [0; g*sin(theta); g*cos(theta)];
Ff = -(dot(F, es)/m);%eqn. 2.13
end
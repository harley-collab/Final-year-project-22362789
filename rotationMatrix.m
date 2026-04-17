function R = rotationMatrix(theta,phi)%eqn. 2.3
Rz = [ cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1 ];
Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R = Rz*Rx;
end
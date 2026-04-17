function [value, isterminal, direction] = combinedEvents(~,y)
theta = y(7);

value(1) = theta - pi*(99/200); %ring effectively parrallel to ground
isterminal(1) = 1;   
direction(1) = +1;   

value(2) = theta; % ring fully upright
isterminal(2) = 1;   
direction(2) = -1; 
end
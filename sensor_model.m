%Sensor Model
%Returns sensor output for a state vector
function [y] = sensor_model(x)

 y = [wrapToPi(atan2((x(5)-x(2)),(x(4)-x(1))) - x(3));
      sqrt((x(1)-x(4))^2 + (x(2)-x(5))^2);
      wrapToPi(atan2((x(2) - x(5)),(x(1)-x(4))) - x(6));
      x(4);
      x(5)];
%   y = [atan2((x(5)-x(2)),(x(4)-x(1))) - x(3);
%       sqrt((x(1)-x(4))^2 + (x(2)-x(5))^2);
%       atan2((x(2) - x(5)),(x(1)-x(4))) - x(6);
%       x(4);
%       x(5)];

end
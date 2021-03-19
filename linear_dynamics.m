function [delta_x,delta_y] = linear_dynamics(F,H,delta_x_up,delta_t,steps)
    
    delta_x = zeros(6,steps+1);
    delta_y = zeros(5,steps);
    for k = 1:steps
            %Time update
            if k == 1
            delta_x(:,k) = delta_x_up;
               
            else
            delta_x_up = F(:,:,k-1)*delta_x_up;
            delta_y_up = H(:,:,k)*delta_x_up;
            delta_x(:,k) = delta_x_up;
            delta_y(:,k) = delta_y_up;
            end
    end
end

function [delta_x_estimate,delta_x_no_meas,delta_x_covariance,p_plus_val,s_val] = Linearized_kalman_filter(F,G,H,Omega,R,Q,delta_x_up,p0,delta_y,delta_t,steps)
    
    delta_x_estimate = zeros(6,steps+1);
    delta_x_no_meas = zeros(6,steps);    
    delta_x_covariance = zeros(6,6,steps+1);
    p_plus_val = zeros(6,6,steps+1);
    s_val = zeros(5,5,steps);
    I = eye(6);
    for k = 1:steps
            %Time update
            if k == 1
               delta_x_estimate(:,k) = delta_x_up;
               p_plus_val(:,:,k) = p0;

            else
            delta_x_minus = F(:,:,k-1)*delta_x_up;
            delta_x_minus(3) = wrapToPi(delta_x_minus(3)); 
            delta_x_minus(6) = wrapToPi(delta_x_minus(6));
            delta_x_no_meas(:,k-1) = delta_x_minus;    
            pi_minus = F(:,:,k-1)*p0*F(:,:,k-1)' + Omega*Q*Omega';
            S = (H(:,:,k)*pi_minus*H(:,:,k)') + R;
            k_update = pi_minus*H(:,:,k)'/S;
            %mean_minus_col = [mean_minus];
            %S_col(:,:,i) = S; 
            %Measurement update
            delta_y_minus = H(:,:,k)*delta_x_minus;
            delta_y_minus(1) =wrapToPi(delta_y_minus(1));
            delta_y_minus(3) =wrapToPi(delta_y_minus(3));
            innovation = delta_y(:,k-1) - delta_y_minus;
            innovation(1) = wrapToPi(innovation(1));
            innovation(3) = wrapToPi(innovation(3));
            delta_x_plus = delta_x_minus + k_update*(innovation);
            delta_x_plus(3) = wrapToPi(delta_x_plus(3)); 
            delta_x_plus(6) = wrapToPi(delta_x_plus(6));
            pk_plus = (I - k_update*H(:,:,k))*pi_minus;
            delta_x_estimate(:,k) = delta_x_plus;
            delta_x_covariance(:,:,k) = pk_plus; 
            delta_x_up = delta_x_plus;
            p0 = pk_plus;
            p_plus_val(:,:,k) = pk_plus;
            s_val(:,:,k-1) = S;
            end
    end
end

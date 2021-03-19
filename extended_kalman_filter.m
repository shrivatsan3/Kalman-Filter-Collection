function [x_estimate,p,y_no_meas,s_val] = extended_kalman_filter(L,x0,p0,u_nominal,Omega,Q,R,y_synthetic,delta_t,steps)
x_estimate = zeros(6,steps+1);
p = zeros(6,6,steps+1);
y_no_meas = zeros(5,steps);    
s_val = zeros(5,5,steps);
I = eye(6);
for k = 0:steps
    if k == 0
        x_estimate(:,k+1) = x0;
        p(:,:,k+1) = p0;
        
    else
        %for k = 1
        intermediate_t = delta_t*(k-1):0.01:k*delta_t; % Perform integration at each time step with addition of noise
        [discard, x] = ode45(@(t,x) non_linear_dynamics(t,x,u_nominal,L,zeros(6,1)),intermediate_t,x0);
        x_minus_estimate = x(11,:)'; % only need value at the end of integration
        x_minus_estimate(3) = wrapToPi(x_minus_estimate(3));
        x_minus_estimate(6) = wrapToPi(x_minus_estimate(6));
        [A_tilde, B_tilde, ~, D_tilde] = linearize(x0, u_nominal,L); 
        F = I + delta_t*A_tilde;
        G = delta_t*B_tilde;
        
        p_minus = F*p0*F' + Omega*Q*Omega';

        y_minus_estimate = sensor_model(x_minus_estimate);
        y_minus_estimate(1) = wrapToPi(y_minus_estimate(1));
        y_minus_estimate(3) = wrapToPi(y_minus_estimate(3));
        y_no_meas(:,k) = y_minus_estimate;    
        
        [~, ~, C_tilde, ~] = linearize(x_minus_estimate, u_nominal,L); 
        H = C_tilde;
        e_update = y_synthetic(:,k)-y_minus_estimate;
        e_update(1) = wrapToPi(e_update(1));
        e_update(3) = wrapToPi(e_update(3));
        S = (H*p_minus*H')+R;
        s_val(:,:,k) = S;
        k_update = (p_minus*H')/S;
        x_plus_estimate = x_minus_estimate + (k_update*e_update);
        x_plus_estimate(3) = wrapToPi(x_plus_estimate(3));
        x_plus_estimate(6) = wrapToPi(x_plus_estimate(6));
        
        p_plus = (I-k_update*H)*p_minus;
        x_estimate(:,k+1) = x_plus_estimate;
        p(:,:,k+1) = p_plus;
        x0 = x_plus_estimate;
        p0= p_plus;
    end
end
end

    
        
        
% Function to emulate the truth model
% Returns time stamped state vectors
% Calculates truth values by adding noise to states formed by propagation
% of nominal values

function [x_truth, y_synthetic] = truth_model(delta_t, steps, L, x0, u_nominal, Qtrue, Rtrue)
    S_V = chol(Qtrue); % Standard deviation for process noise
    R_V = chol(Rtrue);
    x_truth = zeros(6,steps+1);
    y_synthetic = zeros(5,steps);
    for k = 1:steps+1 % k = 1 corresponds to t=0
        Process_noise_sample = S_V*randn(6,1);
        meas_noise_sample = R_V*randn(5,1);
        
        if k == 1
            x_truth(:,k) = x0; 
           
        else
            intermediate_t = delta_t*(k-1):0.01:k*delta_t; % Perform integration at each time step with addition of noise
            [discard, x] = ode45(@(t,x) non_linear_dynamics(t,x,u_nominal,L,Process_noise_sample),intermediate_t,x0);
            x_truth(:,k) = x(11,:)'; % only need value at the end of integration
            x_truth(3,k) = wrapToPi(x_truth(3,k));
            x_truth(6,k) = wrapToPi(x_truth(6,k));
            y_synthetic(:,k-1) = sensor_model(x_truth(:,k)) + meas_noise_sample; % Find synthetic measurements
            x0=x(11,:)'; % update initial condition for next iteration
        end
    end
end

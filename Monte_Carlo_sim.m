%ASEN 5044
%Final Project
%NEES & NIS TESTING
function[error_x_UKF, error_y_UKF, p,s_val,steps] = Monte_Carlo_sim()

load('cooplocalization_finalproj_KFdata.mat');

% Discretization parameters
delta_t = 0.1; % Time for discretization
steps = 1000; % Discrete time steps for simulation
L = 0.5;
t_span = 0:delta_t:(delta_t*steps);
%nominal system parameters to be linearized about
x0_nominal = [10, 0, pi/2, -60, 0, -pi/2]';
u_nominal = [2,-pi/18, 12, pi/25]';
x_per = [0,1,0,0,0,0.1]';

% Get truth values for truth model testing
[x_truth, y_synthetic] = truth_model(delta_t, steps, L, x0_nominal + x_per, u_nominal, Qtrue, Rtrue);

% Get nominal trajectory
[~,x_nominal] = ode45(@(t,x) non_linear_dynamics(t,x,u_nominal,L,zeros(6,1)),t_span,x0_nominal);
x_nominal = x_nominal';
x_nominal(3,:) = wrapToPi(x_nominal(3,:));
x_nominal(6,:) = wrapToPi(x_nominal(6,:));

% Get sensor outputs for nominal 
y_nominal = zeros(5,steps);
for k = 2:steps+1 
    y_nominal(:,k-1) = sensor_model(x_nominal(:,k));
end
ydata = ydata(:,2:end);
delta_y = ydata - y_nominal; % Find deviation of synthetic data from nominal
delta_y(1) = wrapToPi(delta_y(3));
delta_y(3) = wrapToPi(delta_y(3));
Gamma = eye(6,6);
Omega = delta_t*Gamma;
I = eye(6,6);
F = zeros(6,6,1000);
G = zeros(6,4,1000); 
H = zeros(5,6,1000);
% % 
% %Find the jacobians at each time step
% for k = 1:steps+1
% [A_tilde, B_tilde, C_tilde, D_tilde] = linearize(x_nominal(:,k), u_nominal,L); 
% F(:,:,k) = I + delta_t*A_tilde;
% G(:,:,k) = delta_t*B_tilde;
% H(:,:,k) = C_tilde;
% end
% 
% p0 = Qtrue;
% [delta_x_estimate,delta_x_no_meas,delta_x_covariance,p_plus_val_LKF,s_val_LKF] = Linearized_kalman_filter(F,G,H,Omega,10*Rtrue,500*Qtrue,x_per,p0,delta_y,delta_t,steps);
% x_estimate = x_nominal + delta_x_estimate;
% x_estimate(3,:) = wrapToPi(x_estimate(3,:));
% x_estimate(6,:) = wrapToPi(x_estimate(6,:));
% 
% x_minus_estimate = x_nominal(:,2:end) + delta_x_no_meas;
% 
% y_minus_estimate = zeros(5,steps);
% 
% for k = 1:steps
%     y_minus_estimate(:,k) = sensor_model(x_minus_estimate(:,k));
% end
% 
% error_x_LKF = x_truth - x_estimate;
% error_y_LKF = ydata - y_minus_estimate;
% 
% delta_x_auto_cor = zeros(6,steps+1);
% for k=1:steps+1
%     delta_x_auto_cor(1,k) = delta_x_covariance(1,1,k);
%     delta_x_auto_cor(2,k) = delta_x_covariance(2,2,k);
%     delta_x_auto_cor(3,k) = delta_x_covariance(3,3,k);
%     delta_x_auto_cor(4,k) = delta_x_covariance(4,4,k);
%     delta_x_auto_cor(5,k) = delta_x_covariance(5,5,k);
%     delta_x_auto_cor(6,k) = delta_x_covariance(6,6,k);
% end
% 
% [x_estimate_ekf,p_ekf,y_minus_ekf,s_val_ekf] = extended_kalman_filter(L,x0_nominal,p0,u_nominal,Omega,0.1*Qtrue,100*Rtrue,ydata,delta_t,steps);
% x_estimate_ekf(3,:) = wrapToPi(x_estimate_ekf(3,:));
% x_estimate_ekf(6,:) = wrapToPi(x_estimate_ekf(6,:));
% % error_x_EKF = x_truth - x_estimate_ekf;
% % error_y_EKF = y_synthetic - y_minus_ekf;
% 
% 
% p_auto_cor = zeros(6,steps+1);
% for k=1:steps+1
%     p_auto_cor(1,k) = p_ekf(1,1,k);
%     p_auto_cor(2,k) = p_ekf(2,2,k);
%     p_auto_cor(3,k) = p_ekf(3,3,k);
%     p_auto_cor(4,k) = p_ekf(4,4,k);
%     p_auto_cor(5,k) = p_ekf(5,5,k);
%     p_auto_cor(6,k) = p_ekf(6,6,k);
% end

p_auto_cor = zeros(6,steps+1);
for k=1:steps+1
    p_auto_cor(1,k) = p(1,1,k);
    p_auto_cor(2,k) = p(2,2,k);
    p_auto_cor(3,k) = p(3,3,k);
    p_auto_cor(4,k) = p(4,4,k);
    p_auto_cor(5,k) = p(5,5,k);
    p_auto_cor(6,k) = p(6,6,k);
end


Q = Qtrue;
p0 = 1000*Qtrue;
p0(1,1) = 1*p0(1,1);
p0(2,2) = 0.01*p0(2,2);


Q(1,1) = 1000*Q(1,1);
% Q(1,2) = 10*Q(1,2);
% Q(1,3) = 0.01*Q(1,3);
Q(2,2) = 1000*Q(2,2);
% Q(2,1) = 10*Q(2,1);
% Q(2,3) = 0.01*Q(2,3);
Q(3,3) = 0.01*Q(3,3);
% Q(3,1) = 0.01*Q(3,1);
% Q(3,2) = 0.01*Q(3,2);
[x_estimate_ukf,p,y_no_meas, s_val] = unscented_kalman_filter(L,delta_t,steps,x0_nominal,p0,ydata,u_nominal,0.001*Q,Rtrue);
x_estimate_ukf(3,:) = wrapToPi(x_estimate_ukf(3,:));
x_estimate_ukf(6,:) = wrapToPi(x_estimate_ukf(6,:));
error_x_UKF = x_truth - x_estimate_ukf;
error_y_UKF = y_synthetic - y_no_meas;

end



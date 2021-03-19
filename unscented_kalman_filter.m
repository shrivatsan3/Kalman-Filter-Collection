function [x_estimate_ukf,p,y_no_meas, s_val] = unscented_kalman_filter(L,delta_t,steps,x0,p0,y_synthetic,u_nominal,Q,R)
% Routine to imlpement an unscented kalman filter
% Unscented Kalman filter is used to estimate states for non linear systems
% Don't make changes to kappa and beta
% Experiment with alpha (which controls spread about mean)

    kappa = 0;
    beta = 2;
    n = 6;
    alpha = 0.001;
    lambda = ((alpha^2)*(n+kappa)) - n;
    wm0 = lambda/(lambda+n); % Weights for recombining sigma and gamma points
    wc0 = wm0 + 1 - alpha^2 + beta; % Weights for covariance
    wmi = 0.5*(1/(lambda+n));
    wci = wmi;
    y_no_meas = zeros(5,steps);    
    s_val = zeros(5,5,steps);
    gammai_cur = zeros(5,12);
    x_estimate_ukf = zeros(6,steps+1);
    p = zeros(6,6,steps+1);

    for k = 1:steps+1

        if k==1

        x_estimate_ukf(:,k) = x0;
        p(:,:,k+1) = p0;


        else

        intermediate_t = delta_t*(k-2):0.01:(k-1)*delta_t; % Perform integration at each time step without noise
        
        [sig0,sigi,gamma0,gammai] = sigma_points(L,x0,p0,n,lambda,intermediate_t,u_nominal); %get the sigma points after propagation
                                                           % Also returns gamma at k

        x_minus = wm0*sig0 + wmi*sum(sigi,2); %process uppdate
        x_minus(3) = wrapToPi(x_minus(3)); 
        x_minus(6) = wrapToPi(x_minus(6));
        p_minus = zeros(6,6);
        for i = 0:12
            if i==0
            p_minus = wc0*(sig0-x_minus)*((sig0-x_minus)');
            else
            p_minus = p_minus+ wci*((sigi(:,i)-x_minus)*((sigi(:,i)-x_minus)'));
            end
        end
        p_minus = p_minus + Q; % covariance update

        S = chol(p_minus);
        sig0 = x_minus;
        sig0(3) = wrapToPi(sig0(3));
        sig0(6) = wrapToPi(sig0(6));
        gamma0_cur = sensor_model(sig0); % get gamma at k+1 and update sigma with x_minus and p_minus
        for i = 1:6
            sigi(:,i) = x_minus + (sqrt(n+lambda)*S(i,:)');
            sigi(3,i) = wrapToPi(sigi(3,i));
            sigi(6,i) = wrapToPi(sigi(6,i));
            gammai_cur(:,i) = sensor_model(sigi(:,i)); 
            sigi(:,i+n)= x_minus - (sqrt(n+lambda)*S(i,:)');
            sigi(3,i+n) = wrapToPi(sigi(3,i+n));
            sigi(6,i+n) = wrapToPi(sigi(6,i+n));
            gammai_cur(:,i+n) = sensor_model(sigi(:,i+n));
        end

        y_minus = wm0*gamma0_cur + wmi*sum(gammai_cur,2); % get y_minus at k+1 by using gamma at k+1
        y_minus(1) = wrapToPi(y_minus(1));
        y_minus(3) = wrapToPi(y_minus(3));
        y_no_meas(:,k-1) = y_minus;    
        
        py = zeros(5,5);
        for i = 0:12
            if i==0
            py = wc0*(gamma0_cur-y_minus)*((gamma0_cur-y_minus)');
            else
            py = py+(wci*(gammai_cur(:,i)-y_minus)*((gammai_cur(:,i)-y_minus)'));
            end
        end
        s_val(:,:,k) = py;

        py = py + R; % Covariance of y uses gamma at k+1

        pxy = zeros(6,5); % Cross correlation uses gamma at k

        for i = 0:12
            if i==0
            pxy = wc0*(sig0-x_minus)*((gamma0-y_minus)');
            else
            pxy = pxy+(wci*(sigi(:,i)-x_minus)*((gammai(:,i)-y_minus)'));
            end
        end

        k_update = pxy/py; %gain update

        x_plus = x_minus + k_update*(y_synthetic(:,k-1)-y_minus); %add innovation
        x_plus(3) = wrapToPi(x_plus(3)); 
        x_plus(6) = wrapToPi(x_plus(6));
        
        p_plus = p_minus - (k_update*py*k_update');
        x_estimate_ukf(:,k) = x_plus;
        x0 = x_plus; % update for next iteration
        p0 = p_plus;
        p(:,:,k+1) = p0;


        end
end
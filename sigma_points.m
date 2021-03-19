function [sig0,sigi,gamma0,gammai] = sigma_points(L,xu,p,n,lambda,intermediate_t,u_nominal)
S = chol(p);
sig0 = xu;
sigi = zeros(6,12);
gammai = zeros(5,12);
   
[discard, x] = ode45(@(t,x) non_linear_dynamics(t,x,u_nominal,L,zeros(6,1)),intermediate_t,sig0);
sig0 = x(11,:)';
sig0(3) = wrapToPi(sig0(3));
sig0(6) = wrapToPi(sig0(6));
gamma0 = sensor_model(sig0);

for i = 1:6
    sigi(:,i) = xu + (sqrt(n+lambda)*S(i,:)');
    sigi(3,i) = wrapToPi(sigi(3,i));
    sigi(6,i) = wrapToPi(sigi(6,i));
    [discard, x] = ode45(@(t,x) non_linear_dynamics(t,x,u_nominal,L,zeros(6,1)),intermediate_t,sigi(:,i));
    sigi(:,i) = x(11,:)';
    sigi(3,i) = wrapToPi(sigi(3,i));
    sigi(6,i) = wrapToPi(sigi(6,i));
    gammai(:,i) = sensor_model(sigi(:,i));
    sigi(:,i+n)= xu - (sqrt(n+lambda)*S(i,:)');
    sigi(3,i+n) = wrapToPi(sigi(3,i+n));
    sigi(6,i+n) = wrapToPi(sigi(6,i+n));
    [discard, x] = ode45(@(t,x) non_linear_dynamics(t,x,u_nominal,L,zeros(6,1)),intermediate_t,sigi(:,i+n));
    sigi(:,i+n) = x(11,:)';
    sigi(3,i+n) = wrapToPi(sigi(3,i+n));
    sigi(6,i+n) = wrapToPi(sigi(6,i+n));
    gammai(:,i+n) = sensor_model(sigi(:,i+n));
    end

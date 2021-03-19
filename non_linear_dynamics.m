% Function to Simulate high fidelity non linear dynamics
% Input: Initial state conditions, control input, Noise inputs(if any) 

function x_dot = non_linear_dynamics(t,x,u,L,Process_noise_sample)

x_dot = [u(1)*cos(x(3)) + Process_noise_sample(1); 
         u(1)*sin(x(3))+ Process_noise_sample(2);
         (u(1)/L)*tan(u(2))+ Process_noise_sample(3);
         u(3)*cos(x(6))+ Process_noise_sample(4);
         u(3)*sin(x(6))+ Process_noise_sample(5);
         u(4)+ Process_noise_sample(6)];
end



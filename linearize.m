function [jacob_F_x, jacob_F_u, jacob_H_x, jacob_H_u] = linearize(x,u,L) 

jacob_F_x = [0 0 -u(1)*sin(x(3)) 0 0 0;
           0 0 u(1)*cos(x(3)) 0 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 -u(3)*sin(x(6));
           0 0 0 0 0 u(3)*cos(x(6));
           0 0 0 0 0 0];

jacob_F_u = [cos(x(3)) 0 0 0;
             sin(x(3)) 0 0 0;
             (1/L)*tan(u(2)) (u(1)/L)*(sec(u(2))^2) 0 0;
             0 0 cos(x(6)) 0;
             0 0 sin(x(6)) 0;
             0 0 0 1];
         
jacob_H_x = zeros(5,6);         
jacob_H_x(1,1) = (x(5)-x(2))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(1,2) = (x(1)-x(4))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(1,3) = -1;
jacob_H_x(1,4) = (x(2)-x(5))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(1,5) = (x(4)-x(1))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(2,1) = (x(1)-x(4))/sqrt((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(2,2) = (x(2)-x(5))/sqrt((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(2,4) = (x(4)-x(1))/sqrt((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(2,5) = (x(5)-x(2))/sqrt((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(3,1) = (x(5)-x(2))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(3,2) = (x(1)-x(4))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(3,4) =(x(2)-x(5))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(3,5) = (x(4)-x(1))/((x(4)-x(1))^2+(x(5)-x(2))^2);
jacob_H_x(3,6) = -1;
jacob_H_x(4,4) = 1;
jacob_H_x(5,5) = -1;
            
jacob_H_u = zeros(5,4);

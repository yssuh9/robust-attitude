function [euler] = quaternion2euler(quat) 
% function euler = quaternion2euler(q)

N = size(quat,2);

euler = zeros(3,N);
for i = 1:N
    q0 = quat(1,i);
    q1 = quat(2,i);
    q2 = quat(3,i);
    q3 = quat(4,i);
    euler(3,i) = atan2( 2 * q1 * q2 + 2 * q0 * q3 , 2 * q0 * q0 + 2 * q1 * q1 - 1);
    euler(2,i) = asin(-2*q1*q3 + 2*q0*q2);
    euler(1,i) = atan2( 2 * q2 * q3 + 2 * q0 * q1 , 2 * q0 * q0 + 2 * q3 * q3 - 1);
end
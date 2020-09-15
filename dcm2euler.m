function [euler] = dcm2euler(R) 
% function euler = quaternion2euler(q)

euler = zeros(3,1);
euler(2) = -asin(R(1,3));
euler(1) = atan2(R(2,3),R(3,3));
euler(3) = atan2(R(1,2),R(1,1));
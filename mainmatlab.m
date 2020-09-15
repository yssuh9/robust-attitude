% ------------------------------------------------------------------------
%   robust convex attitude optimization
%        this program uses cvx (http://cvxr.com/cvx/)
%        
%      
% 
%   Young Soo Suh (yssuh@ulsan.ac.kr)
%   University of Ulsan, Korea
% ------------------------------------------------------------------------

% simulation mode 
% 1: simulation without magnetic disturbance
% 2: simulation with magnetic disturbance (at 5.3 second, 11.2 second)
% 3: large initial acclerometer noise (first five data)
simulation_mode = 2;

% data number (there are 50 data sets)
% 1 <= data_num <= 50
data_num = 33;

% ------------------------------------------------------------------------
% Basic parameters
% ------------------------------------------------------------------------
% sampling period
T = 0.01;
% gravitation magnitude
g = 9.8;
D2R = pi/180;
R2D = 180 / pi;

% magnetic dip angle about -50 degree
% Ulsan, Korea (lattitude 35.32N, longitude 129.19E 
mu = -50 * pi / 180;
mtilde = [ cos(mu) ; 0 ; sin(mu) ];
gtilde = [ 0 ; 0 ; g];

% Sensor noises
rg = (0.01)^2;
ra = (0.05)^2;
rm = (0.01)^2;

% ------------------------------------------------------------------------
% Simulation results
% ------------------------------------------------------------------------
fprintf('simulation mode %d  data number %d\n',simulation_mode,data_num);

error0 = zeros(3,729);
error1 = zeros(3,729);
error2 = zeros(3,729);
error3 = zeros(3,729);
error4 = zeros(3,729);
J1 = zeros(1,729);
J2 = zeros(1,729);

H = zeros(6,6);
H(1:3,4:6) = g * eye(3);
H(4:6,1:3) = cos(mu) * eye(3);
H(4:6,4:6) = sin(mu) * eye(3);
R = diag([ra ra ra rm rm rm]);
P0 = inv(H' * inv(R) * H);
invP0 = H' * inv(R) * H;

% data reading
%  ya : accelerometer
%  yg : gyroscope
%  ym : magnetic sensor 
filename = sprintf('simdata1\\%d.mat',data_num);
load(filename);

N = size(ya,2);
eulertrue = quaternion2euler(q);
if ( simulation_mode == 2 )
    ym = ymdist;
elseif ( simulation_mode == 3)
    for j = 1:5
        ya(:,j) = ya(:,j) + 2 * randn(3,1);
    end
end

for i = 1:N
    ym(:,i) = ym(:,i) / norm(ym(:,i));
    ya(:,i) = g * ya(:,i) / norm(ya(:,i));
end


% --------------------------------------------------------------------
% 0: static estimation
% --------------------------------------------------------------------
[q4] = static_estimation(ya,ym,ra,rm);
euler0 = quaternion2euler(q4);

% --------------------------------------------------------------------
% 1,2,3: forward-backward smoother
%   1: forward
%   2: backward
%   3: smoother
% --------------------------------------------------------------------
[q4f,q4b,q4s] =  compute_forward_backward(ya,yg,ym,ra,rg,rm,T,g,mu);

euler1 = quaternion2euler(q4f);
euler2 = quaternion2euler(q4b);
euler3 = quaternion2euler(q4s);

% --------------------------------------------------------------------
% 5: proposed method
% --------------------------------------------------------------------
alpha =1;
delta_equality = 0.00000001;
max_iteration = 5;

[x, mag1, mag2, orth, status, A, B, b] = convex_optimization4(ya,ym,yg,ra,rm,rg,g,mu,T,alpha);
xstar = x;
fprintf('status %d mag1 %7.4f mag2 %7.4f orth %7.4f\n',status,sum(mag1),sum(mag2),sum(orth));
J1 = norm(A* x - b);
fprintf('error: %7.4f constraint %7.4f\n', J1, 10000*(sum(mag1)+sum(mag2)+sum(orth)) );

iter_count = 0;
while ( (sum(mag1)+sum(mag2)+sum(orth) > delta_equality) && (iter_count <= max_iteration) )
    [x, mag1, mag2, orth, status] = convex_optimization_local4(A,B,b,x,alpha);
    fprintf('--- status %d mag1 %7.4f mag2 %7.4f orth %7.4f total %7.4f\n',status,sum(mag1),sum(mag2),sum(orth),sum(mag1)+sum(mag2)+sum(orth));
    J2 = norm(A* x - b);
    fprintf('--- error: %7.4f constraint %7.4f\n', J2, 10000*(sum(mag1)+sum(mag2)+sum(orth)) );
    iter_count = iter_count + 1;
end
euler4 = zeros(3,N);
for i = 1:N
    c1 = x(6*(i-1)+1:6*(i-1)+3) / norm(x(6*(i-1)+1:6*(i-1)+3));
    c3 = x(6*(i-1)+4:6*(i-1)+6) / norm(x(6*(i-1)+4:6*(i-1)+6));
    c2 = cross(c3,c1);
    C = [ c1 , c2 , c3 ];
    euler4(:,i) = dcm2euler(C);
end

% visualization (simulation mode == 2)
tt = 0:T:(N-1)*T;

figure(1);
subplot(3,1,1);
plot(tt,eulertrue(1,:),'r:',tt,euler3(1,:),'b');
title('standard smoother');
legend('true','estimated');

z1 = axis;
ylabel('pitch (rad)','FontSize',9,'FontName','Times');
subplot(3,1,2);
plot(tt,eulertrue(2,:),'r:',tt,euler3(2,:),'b');
z2 = axis;
ylabel('roll (rad)','FontSize',9,'FontName','Times');
subplot(3,1,3);
plot(tt,eulertrue(3,:),'r:',tt,euler3(3,:),'b');
z3 = axis;
ylabel('yaw (rad)','FontSize',9,'FontName','Times');

figure(2);
subplot(3,1,1);
plot(tt,eulertrue(1,:),'r:',tt,euler4(1,:),'b');
title('proposed method');
legend('true','estimated');
axis(z1);
ylabel('pitch (rad)','FontSize',9,'FontName','Times');
subplot(3,1,2);
plot(tt,eulertrue(2,:),'r:',tt,euler4(2,:),'b');
axis(z2);
ylabel('roll (rad)','FontSize',9,'FontName','Times');
subplot(3,1,3);
plot(tt,eulertrue(3,:),'r:',tt,euler4(3,:),'b');
axis(z3);
ylabel('yaw (rad)','FontSize',9,'FontName','Times');


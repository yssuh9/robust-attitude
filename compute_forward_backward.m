function [q4f,q4b,q4s] =  compute_forward_backward(ya,yg,ym,ra,rg,rm,T,g,mu)

N = size(ya,2);
Q = 0.25 * rg * eye(3);
gtilde = [ 0 ; 0; g];
mtilde = [ cos(mu) ; 0 ; sin(mu) ];

% forward filter
q4f = zeros(4,N);
Pff = zeros(3,3,N);
[C,P] = dcmfromyaym2(ya(:,1),ym(:,1),ra,rm);
q4f(:,1) = dcm2quaternion(C);
Pff(:,:,1) = P;
oldomega4 = compute_44(yg(:,1));
H = zeros(6,3);
for i = 2:N
    Cq = quaternion2dcm(q4f(:,i-1));
    A = -vec2product(yg(:,i-1));
    dA = eye(3) + A * T + A * A * T^2 / 2;
    Qd = Q*T + 0.5*T*T*A*Q + 0.5*T*T*Q*A';
    Qd = 0.5*(Qd + Qd');
    P = dA * P * dA' + Qd;

    yyg = yg(:,i-1);
    omega4 = compute_44(yyg);
    q4f(:,i) = ( eye(4) + 0.75 * omega4 * T - 0.25 * oldomega4 * T - (1/6) * norm(yyg)^2 * T^2 * eye(4) - (1/24) * omega4 * oldomega4 * T^2  - (1/48)*norm(yyg)^2 *omega4 * T^3) * q4f(:,i-1);
%     q4f(:,i) =  expm( 0.5 * omega4 * T) * q4f(:,i-1);
    q4f(:,i) = q4f(:,i) / norm(q4f(:,i));
    oldomega4 = omega4;
    Cq = quaternion2dcm(q4f(:,i));

    H = [  2 * vec2product(Cq*gtilde) ;  2 * vec2product(Cq*mtilde) ];
    z = [ ya(:,i) - Cq * gtilde ; ym(:,i) - Cq * mtilde ];
    K = P * H' * inv(H * P * H' + diag([ra ra ra rm rm rm]));
    x = K * z;
    qe = [ 1 ;  x];
    q4f(:,i) = quaternionmul(q4f(:,i),qe);
    q4f(:,i)=  q4f(:,i) / norm(q4f(:,i));
    P = (eye(3) - K*H) * P;
    P = (P + P') / 2;
    Pff(:,:,i) = P;
end


% backward filter
q4b = zeros(4,N);
[C,P] = dcmfromyaym2(ya(:,end),ym(:,end),ra,rm);
q4b(:,end) = dcm2quaternion(C);


oldomega4 = compute_44(yg(:,end));
H = zeros(6,3);
q4bb = zeros(4,N);
Pbb = zeros(3,3,N);
Pbb(:,:,end) = P;
for i = N-1:-1:1
    Cq = quaternion2dcm(q4b(:,i+1));
    A = vec2product(yg(:,i+1));
    dA = eye(3) + A * T + A * A * T^2 / 2;
    Qd = Q*T + 0.5*T*T*A*Q + 0.5*T*T*Q*A';
    Qd = 0.5*(Qd + Qd');
    P = dA * P * dA' + Qd;

    yyg = yg(:,i+1);
    omega4 = compute_44(yyg);
    q4b(:,i) = ( eye(4) - 0.75 * omega4 * T + 0.25 * oldomega4 * T - (1/6) * norm(yyg)^2 * T^2 * eye(4) - (1/24) * omega4 * oldomega4 * T^2  + (1/48)*norm(yyg)^2 *omega4 * T^3) * q4b(:,i+1);
%     q4b(:,i) = expm(-0.5*omega4*T) * q4b(:,i+1);
    q4b(:,i) = q4b(:,i) / norm(q4b(:,i));
    q4bb(:,i) = q4b(:,i);
    Pbb(:,:,i) = P;

    oldomega4 = omega4;
    Cq = quaternion2dcm(q4b(:,i));

    H = [  2 * vec2product(Cq*gtilde) ;  2 * vec2product(Cq*mtilde) ];
    z = [ ya(:,i) - Cq * gtilde ; ym(:,i) - Cq * mtilde ];
    K = P * H' * inv(H * P * H' + diag([ra ra ra rm rm rm]));
    x = K * z;
    qe = [ 1 ;  x];
    q4b(:,i) = quaternionmul(q4b(:,i),qe);
    q4b(:,i)=  q4b(:,i) / norm(q4b(:,i));
    P = (eye(3) - K*H) * P;
    P = (P + P') / 2;
end


% smoother
R = zeros(3,3,2);
q4s = zeros(4,N);
for i = 1:N
    q = [ q4f(:,i) , q4bb(:,i) ];
    R(:,:,1) = Pff(:,:,i);
    R(:,:,2) = Pbb(:,:,i);
    q4s(:,i) = quaternionaveraging2(q,R);
end



end


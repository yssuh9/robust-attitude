function [qbar,Rbar] = quaternionaveraging2(q,R)

N = size(q,2);

M = zeros(4,4);
for i = 1:N
    qi = q(:,i);
    Ri = R(:,:,i);
    if ( min(eig(Ri)) <= 0 )
        disp('something wrong');
        pause
    end
    Xi = [ -qi(2:4) , qi(1)*eye(3) + vec2product(qi(2:4)) ];
    M = M + Xi' * inv(Ri) * Xi;
end

M = (M+M')/ 2;
[V,D] = eig(-M);

[maxvalue,maxindex] = max(diag(D));

qbar = V(:,maxindex);

Xi = [ -qbar(2:4) , qbar(1)*eye(3) + vec2product(qbar(2:4)) ];
Rbar = inv( Xi * M * Xi' );




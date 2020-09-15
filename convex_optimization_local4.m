function [y, mag1, mag2, orth, status] = convex_optimization_local4(A,B,b,x,alpha)

    N = length(x) / 6;
    bb = b - A * x;
    
    cvx_begin quiet
        variables z(2*3*N) y(2*3*N);
        minimize( norm( A * z + B * y - bb) + alpha*norm(y,1) );
% 
%         variable z(2*3*N);
%         minimize( norm( A * z - bb) );

        subject to
            for i = 1:N
                [ 2 * x(6*(i-1)+1:6*(i-1)+3)' , zeros(1,3) ; ...
                    x(6*(i-1)+4:6*(i-1)+6)' , x(6*(i-1)+1:6*(i-1)+3)' ; ...
                    zeros(1,3) , 2 * x(6*(i-1)+4:6*(i-1)+6)' ] * [ z(6*(i-1)+1:6*(i-1)+3) ;  z(6*(i-1)+4:6*(i-1)+6) ] ...
                    == [ 1 - x(6*(i-1)+1:6*(i-1)+3)' * x(6*(i-1)+1:6*(i-1)+3) ; - x(6*(i-1)+1:6*(i-1)+3)' * x(6*(i-1)+4:6*(i-1)+6) ; 1 - x(6*(i-1)+4:6*(i-1)+6)' * x(6*(i-1)+4:6*(i-1)+6) ];
            end
    cvx_end
    
    
    y = x + z;
    status = strcmp(cvx_status,'Solved');
    
    mag1 = zeros(1,N);
    mag2 = zeros(1,N);
    for i = 1:N
        mag1(i) = (1 - norm( y(6*(i-1)+1:6*(i-1)+3) ))^2;
        mag2(i) = (1 - norm( y(6*(i-1)+4:6*(i-1)+6) ))^2; 
    end

    orth = zeros(1,N);
    for i = 1:N
        orth(i) = (y(6*(i-1)+1:6*(i-1)+3)' * y(6*(i-1)+4:6*(i-1)+6))^2; 
    end
end


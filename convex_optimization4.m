function [x, mag1, mag2, orth, status, A, B, b] = convex_optimization4(ya,ym,yg,ra,rm,rg,g,mu,T,alpha)

    % x in R^(6N)

    N = size(ya,2);
    % ya
    A1 = zeros(N*3,2*3*N);
    B1 = zeros(N*3,2*3*N);
    b1 = zeros(N*3,1);
    for i = 1:N
        A1((i-1)*3+1:i*3,6*(i-1)+4:6*(i-1)+6) = (1/sqrt(ra)) * eye(3) * g; 
        B1((i-1)*3+1:i*3,6*(i-1)+1:6*(i-1)+3) = (1/sqrt(ra)) * eye(3);
        b1((i-1)*3+1:i*3) = (1/sqrt(ra)) * ya(:,i);
    end

    % ym
    A2 = zeros(N*3,2*3*N );
    B2 = zeros(N*3,2*3*N );
    b2 = zeros(N*3,1);
    for i = 1:N
        A2((i-1)*3+1:i*3,6*(i-1)+1:6*(i-1)+3) = (1/sqrt(rm)) * eye(3) * cos(mu); 
        A2((i-1)*3+1:i*3,6*(i-1)+4:6*(i-1)+6) = (1/sqrt(rm)) * eye(3) * sin(mu); 
        B2((i-1)*3+1:i*3,6*(i-1)+4:6*(i-1)+6) = (1/sqrt(rm)) * eye(3);
        b2((i-1)*3+1:i*3) = (1/sqrt(rm)) * ym(:,i);
    end    
    
    % omega C1
    A3 = zeros((N-1)*3,2*3*N );
    b3 = zeros((N-1)*3,1);
    for i = 1:N-1
%         MMM = expm(-vec2product(yg(:,i)*T));
        MMM = expm(-vec2product(yg(:,i)*T));
        A3((i-1)*3+1:i*3,6*(i-1)+1:6*(i-1)+3) = - (1/(sqrt(rg)*sqrt(T))) * MMM;
        A3((i-1)*3+1:i*3,6*i+1:6*i+3) = (1/(sqrt(rg)*sqrt(T))) * eye(3); 
    end    
    
    % omega C3
    A4 = zeros( (N-1)*3 , 2*3*N );
    b4 = zeros( (N-1)*3 , 1);
    for i = 1:N-1
         MMM = expm(-vec2product(yg(:,i)*T));

        A4((i-1)*3+1:i*3,6*(i-1)+4:6*(i-1)+6) = - (1/(sqrt(rg)*sqrt(T))) *MMM; 
        A4((i-1)*3+1:i*3,6*i+4:6*i+6) = (1/(sqrt(rg)*sqrt(T))) * eye(3); 
    end    
    
    A = [ A1 ; A2 ; A3 ; A4 ];
    b = [ b1 ; b2 ; b3 ; b4 ];
    B = [ B1 ; B2 ; zeros(2*(N-1)*3,6*N) ];
    
    
    cvx_begin quiet
        variables x(2*3*N) z(2*3*N);
        minimize( norm( A * x + B * z - b ) + alpha * norm( z, 1));
%         variable x(2*3*N)
%         minimize( norm( A * x - b) );
        subject to
            for i = 1:N
                [ 1 0 x(6*(i-1)+1:6*(i-1)+3)' ; ... 
                  0 1 x(6*(i-1)+4:6*(i-1)+6)' ; ...
                  x(6*(i-1)+1:6*(i-1)+3) x(6*(i-1)+4:6*(i-1)+6) eye(3) ] == semidefinite(5);
            end
    cvx_end
    
    status = strcmp(cvx_status,'Solved');
    
    mag1 = zeros(1,N);
    mag2 = zeros(1,N);
    for i = 1:N
        mag1(i) = (1 - norm( x(6*(i-1)+1:6*(i-1)+3) ))^2;
        mag2(i) = (1 - norm( x(6*(i-1)+4:6*(i-1)+6) ))^2; 
    end

    orth = zeros(1,N);
    for i = 1:N
        orth(i) = (x(6*(i-1)+1:6*(i-1)+3)' * x(6*(i-1)+4:6*(i-1)+6))^2; 
    end
end


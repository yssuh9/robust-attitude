function [qest,Pqe] = quaternionfromyaym(ya,ym,ra,rm,mode)

%
% mode = 1 : triad
%      = 2  : quest
%

if ( mode == 1 )
    ya = ya / norm(ya);
    ym = ym / norm(ym);

    foo = cross(ya,ym);
    foo = foo / norm(foo);
    C = [ -cross(ya,foo) , foo, ya ];

    ra = ra / (9.8^2);
    qest = dcm2quaternion(C);
    s1 = ya;
    s2 = cross(ya,ym) / norm( cross(ya,ym) );
    s3 = cross(s1,s2);
    Pqe = (1 / norm(s2)^2 ) * ( ra * (ym*ym'+s2*s2') + rm * ya * ya');
%     Pqe = ra*eye(3) + (1/norm(cross(ya,ym))^2) * ( (rm - ra)*ya*ya' + ...
%         ra*(ya'*ya)*(ya*ym'+ ym*ya') );
elseif ( mode == 2 ) 
    gtilde = [ 0 ; 0 ; 9.8 ];
    alpha = 50 * pi / 180;
    mtilde = [ cos(alpha) ; 0 ; -sin(alpha) ];
    if ( (ra == 0) && (rm == 0) )
	a1 = 1;
	a2 = 1;
    else
	a1 = rm / (ra + rm);
    	a2 = ra / (ra + rm);
    end
    M = a1 * ya * gtilde' + a2 * ym * mtilde';
    [u,s,v] = svd(M);
    C = u * diag([1 1 det(u)*det(v)]) * v';
    qest = dcm2quaternion(C);
    f1 = ra*rm / (ra + rm);
    f2 = 1 / ( norm(cross(ya,ym))^2 * (ra+rm) );
    Pqe = 0.25 * (f1 * eye(3) + f2 * ( rm*rm*ya*ya' + ra*ra*ym*ym' + ra*rm*(ya'*ym)*(ya*ym'+ym*ya') ) );
end
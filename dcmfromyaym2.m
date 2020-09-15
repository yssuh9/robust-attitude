function [C,Q] = dcmfromyaym2(ya,ym,ra,rm)
% function [C,Q] = dcmfromyaym2(ya,ym,ra,rm)
% input
%     ra: scalar
%     rm: scalar
% output
%     C : rotation matrix
%     Q = E[qe qe'] 

D2R = (pi / 180);

g = 9.8;
gtilde = [ 0 ; 0 ; g];
mtilde = [ cos(50*D2R)  ; 0 ; -sin(50*D2R) ];


yahat = ya / norm(ya);
yam = cross(ya,ym);
yamhat = yam / norm(yam);
yaamhat = cross(yahat,yamhat);

Q1 = [yahat yamhat yaamhat];

foo1 = cross(gtilde,mtilde);
foo2 = foo1 / norm(foo1);

foo3 = cross(gtilde,foo1);
foo4 = foo3 / norm(foo3);

Q2 = [ gtilde/norm(gtilde) , foo2 , foo4];

C = Q1 * Q2';

a = 0.5 * g  / ( norm(cross(gtilde,mtilde)));
b = -0.5*(gtilde' * mtilde) /  norm(  cross(gtilde,cross(gtilde,mtilde))  );
c = 1 / (-2*g);
d = 1 / (2*g);

Q = Q1 * [ a*a*rm + b*b*ra , 0 , b * d* ra ; ...
             0 , c*c*ra , 0 ; ...
            b*d*ra , 0 , d*d*ra] * Q1';

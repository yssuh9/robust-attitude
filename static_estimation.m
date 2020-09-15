function [q4] = static_estimation(ya,ym,ra,rm)

N = size(ya,2);
q4 = zeros(4,N);
for i = 1:N
    [q4(:,i)] = quaternionfromyaym(ya(:,i),ym(:,i),ra,rm,2);
end

end


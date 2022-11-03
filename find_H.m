function [A,H]=find_H(window_size,order)
% evaluating design matrix A and matrix H from A
k=window_size;
for v=1:2*k+1
    A(v,:)=(v-(k+1)).^(0:order);
end
H=inv(A'*A)*A';
%%
%%Variables and constants:
clear all;
tresh = 0.01;
dif = 99999999;
%A = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];
A = [10 1 1 2; 2 10 2 3 ; 2 1 10 1 ; 3 4 1 10];

B = [10; 10; 10; 10];
X = [0; 0; 0; 0];

%%
%%DLU decomposition:

[m,n] = size(A); %TODO: m should be equal to n anyway, right?
D = zeros(m,n);
L = zeros(m,n);
U = zeros(m,n);

for i = 1:m
    for j = 1:n
        if i == j
            D(i,j) = A(i,j);
        end
        if i > j
            L(i,j) = -A(i,j);
        end
        if i < j
            U(i,j) = -A(i,j);
        end
    end
end

%(D-L)^(-1)
DmLi = inv(D-L);
%%
%%Processing:
idx = 1;
while dif > tresh
   old = X;
   X = DmLi*U*X + DmLi*B;
    
   dif = norm(X - old);   
   conv(idx) = dif;
   idx = idx+1;
end

plot(conv)
%%
%%Debug, if the solution is correct DBG should be near zero.
DBG = X - A\B
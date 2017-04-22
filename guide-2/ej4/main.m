%%
%%Variables and constants:
clear all;
tresh = 0.01;
dif = 99999999;
%Original non-convergent system.
A = [2 -1 10 0; 0 3 -1 8; 10 -1 2 0; -1 11 -1 3];
B = [-11; -11; 6; 25];

%Rearranged, convergent system.

A = [10 -1 2 0;  -1 11 -1 3; 2 -1 10 0; 0 3 -1 8];
B = [6; 25; -11; -11];


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
Di = inv(D);
%%
%%Method matrices:
idx = 1;
E = Di*(L+U);
C = Di*B;

%%Spectral Radii of E (if < 1, the method had convergence).
eigE = eig(E);
pE = max([abs(max(eigE)) abs(min(eigE))]);
disp('Spectral radii pE. If pE<1, the method converges.')

pE

while dif > tresh
   old = X;
   X = E*X + C;
   
   dif = norm(X - old);   
   conv(idx) = dif;
   idx = idx+1;
end

plot(conv)
%%
%%Checking result:
disp('If the solution is correct DBG should be near zero.')
DBG = X - A\B
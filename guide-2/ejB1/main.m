%%
%%Variables and constants:
clear all;
tresh = 0.1;
dif = 99999999;

cantBars = 5;
barWidth = 3;
sysWidth = (2*cantBars+1)*barWidth;
sysHeight = 9;
baseHeight = (4/5)*sysHeight;
barsMinReach = (1/5)*sysHeight;
N = sysWidth*sysHeight; %Matrix size
Tamb = 30;
Tcpu = 60;


%Matrices.
A = zeros(N,N);
B = zeros(N);
X = zeros(N);

idx = 0;
for i = 1: sysHeight
    for j = 1:sysWidth
        
        idx = idx +1;

        
        if i < barsMinReach;
            A(idx,idx) = 1;
            B(idx) = Tamb;
            map(i,j) = -1;
            continue
        end

        if i == sysHeight
            A(idx,idx) = 1;
            B(idx) = Tcpu;
            map(i,j) = 1;
            continue
        end
        

        modBW = rem(j,2*barWidth);
        if (modBW > barWidth) || (modBW == 0) || ( i >= baseHeight && j > 1 && j < sysWidth)
            A(idx,idx) = 4;
            if idx-1 > 1
            A(idx,idx-1) = -1;
            end
            if idx+1 < N
            A(idx,idx+1) = -1;
            end
            if idx - sysWidth > 1
            A(idx,idx-sysWidth) = -1;
            end
            if idx + sysWidth < N
            A(idx,idx+sysWidth) = -1;
            end
            B(idx) = 0;
            map(i,j) = 0;
            continue
        end

        if modBW <= barWidth
            A(idx,idx) = 1;
            B(idx) = Tamb;
            map(i,j) = -1;
            continue
        end

    end
end


%A


%%
%%DLU decomposition:
D = zeros(N,N);
L = zeros(N,N);
U = zeros(N,N);

for i = 1:N
    for j = 1:N
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


%(D)^(-1)
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

plot(conv);
figure



%%
%%Making the physical map.
idx = 0;
for i = 1:sysHeight
    for j = 1:sysWidth
        idx = idx+1;
        S(i,j) = X(idx);
    end
end


mesh(S)
surf(S)
HeatMap(map)



%%
%%Checking result:
%disp('If the solution is correct DBG should be near zero.')
%DBG = X - A\B

%%
%Misc
clear all;
%%
%Variables and constants:
DIF_TRESH = 0.0001;
DIF_CUR = 99999999; %TODO: inf constant?
TIME_MAX = 15;
CONC_MID_START = 0.01;
CONC_BOUND_HIGH = 1;
CONC_BOUND_LOW = 0;
MID_START = 1;
MID_END = 15;
DIFFUSIVITY = 0.01;
DELTA_X = 1;
BOUNDARY_CELLS  = 2;

MID_LONG = MID_END - MID_START;
C_SIZE = (MID_LONG / DELTA_X) + BOUNDARY_CELLS;

C = transpose(zeros(C_SIZE));
B = zeros(C_SIZE, C_SIZE);

%%
%Fill B and C:
for i = 1:C_SIZE
    C(i) = CONC_MID_START;
end
C(1) = CONC_BOUND_LOW;
C(C_SIZE) = CONC_BOUND_HIGH;



% B is filled with the following discretization
% Cn+1(i) = (1/2)*(Cn(i+1) + Cn(i-1))
for i = 2:(C_SIZE-1)
            B(i,i-1) = 0.5;
            B(i,i+1) = 0.5;

end
B(1,1) = 1;
B(C_SIZE,C_SIZE) = 1;

%%
%Compute:
idx = 0;
while DIF_CUR > DIF_TRESH && idx <= TIME_MAX
    idx = idx+1;
    for i = 1:C_SIZE
        graph(idx, i) = C(i);
    end
    
    old = C;
    
    %Direct method --> product computes C_(t+1).
    C = B*C;
    DIF_CUR = norm(C - old);   
    conv(idx) = DIF_CUR;
end

figure;
plot(conv);

figure;
mesh(graph);

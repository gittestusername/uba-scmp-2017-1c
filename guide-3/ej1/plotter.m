
clear all;
M = dlmread('out');
n = size(M);
nc = n(1);

for i = 1:nc
    %figure;

    plot(M(i,:)), axis([0 120 -0.4 0.4])
    drawnow

    %zlabel('Amplitude');
    %ylabel('time');
    %xlabel('x');
    %filename = sprintf('~/%d', i);
end
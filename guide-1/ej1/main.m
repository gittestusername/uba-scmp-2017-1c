%
clear all;

x0=0; y0=0;
a=0.5;
b=0;
c=0.125;
lim_inf=-5;
lim_sup=5;
step=0.1;

N= 1 + (lim_sup-lim_inf)/step;
for k=1:lim_sup
    for i=1:N
        X=i*step +lim_inf;
        for j=1:N
            Y=j*step +lim_inf;
            Z(i,j) = k*exp( - (a*(X-x0)^2 + 2*b*(X-x0)*(Y-y0) + c*(Y-y0)^2)) ;
        end
    end
    pcolor(Z)
    %shading interp;
    file = strcat('plot', num2str(k));
    print(file,'-dbmp')
    %pause(0.1) %in seconds 
    csvwrite('csvlist.dat',Z);
    hold on
     H = csvread('csvlist.dat');
     contour(H,5,'linecolor','w', 'linewidth', 2);
     hold off
end


%
clear all;

x0=0; y0=0;
lim_inf=-5;
lim_sup=5;
step=0.1;
%%Define base u function.
N= 1 + (lim_sup-lim_inf)/step;
    for i=1:N
        x=i*step +lim_inf;
        for j=1:N
            y=j*step +lim_inf;
            u(i,j) = x*exp( -x^2 - y^2 );
        end
    end
    mesh(u);
    shading interp;
    
    for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            dudx(i,j) = (u(i+1,j)-u(i-1,j))/(2*step);
        end
    end
    
    for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            dudy(i,j) = (u(i,j+1)-u(i,j-1))/(2*step);
        end
    end
    
    quiver(dudx,dudy)
%
clear all;

x0=0; y0=0;
lim_inf=-5;
lim_sup=5;
step=0.1;
%%Define base u function matrix
N= 1 + (lim_sup-lim_inf)/step;
    for i=1:N
        x=i*step +lim_inf;
        for j=1:N
            y=j*step +lim_inf;
            u(i,j) = x*exp( -x^2 - y^2 );
        end
    end

    
    %%dudx and dudy matrix
    for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            ux(i,j) = (u(i+1,j)-u(i-1,j))/(2*step);
            uy(i,j) = (u(i,j+1)-u(i,j-1))/(2*step);
        end
    end

        %%d2ud2x and d2ud2y matrix
    for i=2:N-2
        x=i*step +lim_inf;
        for j=2:N-2
            y=j*step +lim_inf;
            ux2(i,j) = (ux(i+1,j)-ux(i-1,j))/(2*step);
            uy2(i,j) = (uy(i,j+1)-uy(i,j-1))/(2*step);
        end
    end
    %%Laplacian matrix (d2x + d2y)
    for i=2:N-2
        x=i*step +lim_inf;
        for j=2:N-2
            y=j*step +lim_inf;
            lap_u(i,j) = ux2(i,j) + uy2(i,j);
        end
    end
    
    %%define analytical laplacian matrix
    for i=2:N-2
        x=i*step +lim_inf;
        for j=2:N-2
            y=j*step +lim_inf;
            an_lap_u(i,j) = (x^2 + y^2 - 2) * 4*x*exp(-x^2 - y^2);
        end
    end   

    surf(lap_u)
    
    
%
clear all;

x0=0; y0=0;
lim_inf=-2;
lim_sup=2;
step=0.05;
%%Define base u function matrixes, and volume mask matrix.
N= 1 + (lim_sup-lim_inf)/step;
    for i=1:N
        x=i*step +lim_inf;
        for j=1:N
            y=j*step +lim_inf;
            for k=1:N
                z = k*step +lim_inf;
                ux(i,j,k) = 2*x;
                uy(i,j,k) = y^2;
                uz(i,j,k) = z^2;
                if x^2 + y^2 + z^2 < 1
                    V(i,j,k) = 1;
                else
                    V(i,j,k) = 0;
                end
            end
        end
    end
 
%%Derivatives

 for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            for k=2:N-1
                z = k*step +lim_inf;
                uxdx(i,j,k) = (ux(i+1,j,k)-ux(i-1,j,k))/(2*step);
                uydy(i,j,k) = (uy(i,j+1,k)-uy(i,j-1,k))/(2*step);
                uzdz(i,j,k) = (uz(i,j,k+1)-uz(i,j,k-1))/(2*step);

            end
        end
 end
        
%%divergencia
  for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            for k=2:N-1
                z = k*step +lim_inf;
                
                div_u(i,j,k) = uxdx(i,j,k) + uydy(i,j,k) + uzdz(i,j,k);
                %div_u(i,j,k) = 2 + 2*y + 2*z;
            end
        end
  end 
  
  
  S = 0;
    for i=2:N-1
        for j=2:N-1
            for k=2:N-1
                S = S + div_u(i,j,k)*step^3*V(i,j,k);
            end
        end
  end   

  S
  
  
%
clear all;

x0=0; y0=0;
lim_inf=-2;
lim_sup=2;
step=0.1;
%%Define base u function matrixes
N= 1 + (lim_sup-lim_inf)/step;
    for i=1:N
        x=i*step +lim_inf;
        for j=1:N
            y=j*step +lim_inf;
            for k=1:N
                z = k*step +lim_inf;
                ux(i,j,k) = -y / ( (x^2 + y^2)^(3/2) );
                uy(i,j,k) = -x / ( (x^2 + y^2)^(3/2) );
                uz(i,j,k) = 0 / ( (x^2 + y^2)^(3/2) );
            end
        end
    end
%%Function plot.
%%Transform into a quiver3 format.
%     idx = 1
%     for i=1:N
%         x=i*step +lim_inf;
%         for j=1:N
%             y=j*step +lim_inf;
%             for k=1:N
%                 z = k*step +lim_inf;
%                 xs(idx) = x;
%                 ys(idx) = y;
%                 zs(idx) = z;
%                 qux(idx) = -y / ( (x^2 + y^2)^(3/2) );
%                 quy(idx) = -x / ( (x^2 + y^2)^(3/2) );
%                 quz(idx) = 0 / ( (x^2 + y^2)^(3/2) );
%                 idx = idx+1;
%             end
%         end
%     end
%quiver3(xs, ys, zs, qux, quy, quz, 10);

%%Derivatives

 for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            for k=2:N-1
                z = k*step +lim_inf;
                uxdy(i,j,k) = (ux(i,j+1,k)-ux(i,j-1,k))/(2*step);
                uxdz(i,j,k) = (ux(i,j,k+1)-ux(i,j,k-1))/(2*step);
                
                uydx(i,j,k) = (uy(i+1,j,k)-uy(i-1,j,k))/(2*step);
                uydz(i,j,k) = (uy(i,j,k+1)-uy(i,j,k-1))/(2*step);

                uzdx(i,j,k) = (uy(i+1,j,k)-uy(i-1,j,k))/(2*step);
                uzdy(i,j,k) = (uy(i,j+1,k)-uy(i,j-1,k))/(2*step);
            end
        end
 end
        
%%rotor
idx = 1;
  for i=2:N-1
        x=i*step +lim_inf;
        for j=2:N-1
            y=j*step +lim_inf;
            for k=2:N-1
                z = k*step +lim_inf;
                
                rot_ux(i,j,k) = uzdy(i,j,k) - uydz(i,j,k);
                rot_uy(i,j,k) = uxdz(i,j,k) - uzdx(i,j,k);
                rot_uz(i,j,k) = uydx(i,j,k) - uxdy(i,j,k);
                
                qrotx(idx) = rot_ux(i,j,k);
                qroty(idx) = rot_uy(i,j,k);
                qrotz(idx) = rot_uz(i,j,k);
                
                xs(idx) = x;
                ys(idx) = y;
                zs(idx) = z;
                
                idx = idx+1;
            end
        end
  end     
        
  %quiver3(xs, ys, zs, qrotx, qroty, qrotz, 100);
  
  %%Norm matrix

  for i=2:N-1
        for j=2:N-1
            for k=2:N-1
                rot_norm(i,j,k) = (rot_ux(i,j,k)^2 + rot_uy(i,j,k)^2 + rot_uz(i,j,k)^2)^(1/2);
            end
        end
  end 
        
  
  %-isosurface(rot_norm, 5)

  
  
      surf1 =   isosurface(rot_norm, 0.5);
      p2 = patch(surf1);
      set(p2,'FaceColor','red','EdgeColor','none','FaceAlpha',0.8);
      
      surf2 =   isosurface(rot_norm, 1);
      p2 = patch(surf2);
      set(p2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.6);
      
      surf3 =   isosurface(rot_norm, 3);
      p2 = patch(surf3);
      set(p2,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.4);
      
      surf3 =   isosurface(rot_norm, 8);
      p2 = patch(surf3);
      set(p2,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.1);
      

camlight; lighting gouraud;
  
  
  
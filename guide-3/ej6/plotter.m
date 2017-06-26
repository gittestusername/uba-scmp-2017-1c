clear all;
data = dlmread('out');
n = size(data);
 data_rows = n(1);
 data_cols = n(2);
 mat_rows = 11;
 mat_cols = data_cols;
 skip = 1;
 for i=1:skip*mat_rows*2:data_rows
     U = data(i:i+mat_rows-1, 1:mat_cols);
     V = data(i+mat_rows: i+2*mat_rows-1, 1:mat_cols);
     
   idx = 16;
   for i=1:mat_rows
       x = i;
       for j=1:data_cols
           y=j;
               %create the arrays for the quiver3 plot.
               xs(idx) = x;
               ys(idx) = y;
               qux(idx) = U(i,j);
               quy(idx) = V(i,j);
               %create matrix for the heatmap plot.
               hm(j,i) = sqrt(qux(idx)*qux(idx) + quy(idx)*quy(idx));
 
               idx = idx+1;
       end
   end
quiver(xs, ys, qux, quy, 5);  
%surf(hm);
title('Function plot');
xlabel('x');
ylabel('y');
 
%figure
hold on 
drawnow
hold off
 end
% for iter = 1:32:(data_rows/mat_rows)
%   ui = data(1+mat_rows*iter-mat_rows:mat_rows*iter, 1:mat_rows);
%   vi = data(mat_rows*iter+1:mat_rows*iter+1 + mat_rows, 1:mat_rows);
%   idx = 1;
%   for i=1:mat_rows
%       x = i;
%       for j=1:data_cols
%           y=j;
%               %create the arrays for the quiver3 plot.
%               xs(idx) = x;
%               ys(idx) = y;
%               qux(idx) = ui(i,j);
%               quy(idx) = vi(i,j);
%               %create matrix for the heatmap plot.
%               hm(j,i) = sqrt(qux(idx)*qux(idx) + quy(idx)*quy(idx));
% 
%               idx = idx+1;
%       end
%   end
%    iter
%    %pause(0.5)
% 
% quiver(xs, ys, qux, quy, 5);
% %surf(hm);
% title('Function plot')
% xlabel('x')
% ylabel('y')
% 
% %figure
% hold on 
% drawnow
% hold off
% 
% 
% end











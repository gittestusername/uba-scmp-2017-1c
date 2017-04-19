%
clear all;

sem = [0; 0; 0];
tresh = 0.001;


dif = 10;

idx = 1;
while dif > tresh
   it = [(2 + sem(2))/4; (6 + sem(1)+ sem(3))/4; (2 + sem(2))/4];
   
   dif = norm(it - sem);
   
   sem = it;
   erJ(idx) = dif;
   idx = idx+1;
end
plot(erJ)



sem = [0; 0; 0];
idx = 1;






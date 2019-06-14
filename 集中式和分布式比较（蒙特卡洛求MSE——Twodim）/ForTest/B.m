function [B]= B(t)
f=0.01*sin(t);
B=[0.3+f 0;
       0 0.3+f];
end
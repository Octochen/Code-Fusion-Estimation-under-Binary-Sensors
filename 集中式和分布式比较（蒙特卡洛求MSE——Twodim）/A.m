function [A]= A(t)
f=0.1*sin(t);
A=[0.85+f 0;
        0 0.75+f];
end


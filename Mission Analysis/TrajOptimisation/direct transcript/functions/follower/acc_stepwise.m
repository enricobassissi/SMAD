%MY ACC STEPWISE MM
function out = acc_stepwise(u,t,t_vec)
k=0;
for i=1:(length(t_vec)-1)
if and(t_vec(i)<=t,t<t_vec(i+1))
    k=i;
end  
end
if k==0
    out=0;
else
    out=(u(k));
end
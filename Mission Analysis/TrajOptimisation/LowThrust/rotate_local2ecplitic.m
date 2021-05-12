
function [Rglobal] = rotate_local2ecplitic(RI,R_sdrp,n_sol,Href)

% R_sdrp = [output.r.*cos(output.l),output.r.*sin(output.l), output.z]

RIv = RI/norm(RI);
hv  = Href/norm(Href);
yv  = cross(hv,RIv);

Arot = [RIv(1) RIv(2) RIv(3)
        yv(1)   yv(2)  yv(3)
        hv(1)   hv(2)  hv(3)];
    
Rglobal = R_sdrp;

for i = 1:n_sol
    %prof: 
    %Rglobal(:,i) = Arot'*R_sdrp(:,i);
    Rglobal(i,:) = Arot'*R_sdrp(i,:)';
end

end

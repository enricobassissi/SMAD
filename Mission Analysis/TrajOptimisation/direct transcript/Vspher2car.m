function [vx,vy,vz] = Vspher2car(az,el,ro,az_,el_,ro_)
vx=cos(el)*cos(az)*ro_+ro*cos(el)*sin(az)*az_+ro*sin(el)*cos(az)*el_;
vy=cos(el)*sin(az)*ro_-ro*cos(el)*cos(az)*az_+ro*sin(el)*sin(az)*el_;
vz=sin(el)*ro_-ro*cos(el)*el_;
end


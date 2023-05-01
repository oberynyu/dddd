function [ffogx]=func_FCRA(thetax,delta,Av,vidx);
for i = 1:length(vidx)
    Ffog_(i) = Av(vidx(i))./(thetax-delta(vidx(i)));
end
Ffog = sum(Ffog_);
Vfog = length(vidx);
eps2     = 0.001;
theta_dw = max(delta(vidx)) ;
for i = 1:length(vidx)
    tmps(i) = Av(vidx(i)) * Vfog/Ffog +  delta(vidx(i));
end
theta_up = sum(tmps);

while abs(theta_up-theta_dw)>eps2
      theta = (theta_up+theta_dw)/2;
      if sum(Av(vidx)/theta -  delta(vidx)) > Ffog
         theta_dw = theta;
      else
         theta_up = theta;   
      end
end
 
thetax = abs(theta_up-theta_dw)/2;

%Output: ffogx
 
ffogx = abs(Av'./(thetax-delta))/4e4;


end






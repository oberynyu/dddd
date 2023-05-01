function [Tv,Tv_local,Tv_fog,Tv_cloud,Av]=func_Tv(fv_fogs,fv_clouds,fv_locals,xv,yv,zv)
global Wf;
global Pv;
global N0;
global Dv;
global lemdav;
global Fflog;
global fvcloud;
global ufuc;
global rfc;
global C;
global speed;
global I;
global M;
global m_;
global L;
global V;
%Calculate the estimated service delay of V-UE v using (13)
gvf       =-10*rand;
taovf     = Wf*log2(1+Pv*gvf/N0);
taofc     = rfc;
Tvf_trans =Dv/taovf;%4
Tvfc_trans=Dv/taofc;%5
Av        =Dv*lemdav;
Tvf_proc  =Av/fv_fogs;%11
lemdaf    =1+2*rand;
uf        =ufuc;
pf        =lemdaf/uf;
Tf_wait   =pf/(uf - lemdaf);%7
Tvc_trans =Tvf_trans + Tvfc_trans;%6
Tvc_proc  =Av/fv_clouds;%12
lemda_cloud= 0.1+rand/2;
C_cloud   = 0.5+rand/2;
u_cloud   = 0.5+rand/2;
p_cloud   =lemda_cloud/C_cloud;
a_cloud   =lemda_cloud/u_cloud; 
ttt=[];
for n = 1:C
    ttt(n)=a_cloud^(n-1)/factorial(n-1);
end
p0        = 1/(sum(ttt) + a_cloud^C/(factorial(C)*(1-p_cloud)));%9
Tc_wait   =a_cloud^C*p0/(factorial(C)*C*u_cloud*(1-p_cloud)^2);%8

 
%
Tv_local  = Av/fv_locals;%10
Tv_fog    = Tvf_trans + Tvf_proc + Tf_wait;%13
Tv_cloud  = Tvc_trans + Tvc_proc + Tc_wait;%13
%
Tv        = xv*Tv_local + yv*Tv_fog + zv*Tv_cloud;%13
pause(0.001)
clc;
clear;
close all;
warning off;
addpath(genpath(pwd));
rng('default');
rng(1);

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
 

% Transmit bandwidth, Wf 180 kHz
Wf = 180;
% Transmit power of V-UE v, Pv 200 mW
Pv = 200;
% The noise power density at the FN, N0 -174 dBm/Hz
N0 = -174;
% Data size of an application of V-UE v, Dv 0.42 MB
Dv = 0.42;
% Processing density of the application of V-UE v, λv 297.62 cycles/bit
lemdav = 297.62;
% Total computation capability of the FN, F% fog 2 G cycles/s
Fflog = 2;
% Cloud processing capability for V-UE v, f% cloudv 5 G cycles/s
fvcloud = 5;
% The service rate of the FN/cloud server, µf /µc 6
ufuc    = 6;
% Wired link rate between the FN and the cloud, rf,c 1 Mb/s
rfc     = 1;
% The number of cloud servers, C 2
C       = 2;
% The average vehicular velocity 70 km/h
speed   = 70;
% The number of fireworks, I 6
I       = 6;
% The number of total explosion sparks, M 4
M       = 4;
% The number of mutation sparks, mˆ 1
m_      = 1;
% The maximum number of iterations of FA, L 100
L       = 50;
V       = 3;

%1：Generate I random fireworks
Theta = zeros(V,3,I,L);
for i=1:I
    for j = 1:V
        tmps            = rand(1,3);
        Theta(j,1:3,i,1)  = [tmps/sum(tmps)];
    end
end
%2：For each firework, randomly allocate fog computation resources
tmps    = rand(1,V);
fv_fog  = Fflog*tmps/sum(tmps);
tmps    = rand(1,V);
fv_cloud= fvcloud*tmps/sum(tmps);
fv_local= 10*ones(1,V);
sv      = 5+5*rand(1,V);
dv      = 500+500*rand(1,V);
xv      = 1/3;
yv      = 1/3;
zv      = 1/3;
%3~16
for v = 1:V
    %Calculate the service delay threshold of V-UE v using (15)
    taov(v)   = dv(v)/sv(v); 
    for i = 1:I
        xv      = Theta(v,1,i,1);
        yv      = Theta(v,2,i,1);
        zv      = Theta(v,3,i,1);
        %Calculate the estimated service delay of V-UE v using (13)
        [Tv,Tv_local,Tv_fog,Tv_cloud,Av]=func_Tv(fv_fog(v),fv_cloud(v),fv_local(v),xv,yv,zv);
        if min(Tv_fog,Tv_cloud) > taov(v)
           Theta(v,:,i,1) = [1,0,0];
        elseif Tv_fog<Tv_cloud
           Theta(v,:,i,1) = [0,1,0];
        else
           Theta(v,:,i,1) = [0,0,1];
        end
        Avs(v,i) = Av;
    end
end



l    = 1;
F0   = 0;
F(1) = F0;
es   = 1e-6;
s    = zeros(L,I);
for l = 1:L
    l
    for v = 1:V
        for i = 1:I
           xv      = Theta(v,1,i,l);
           yv      = Theta(v,2,i,l);
           zv      = Theta(v,3,i,l);
           [Tv(v,i),Tv_local(v,i),Tv_fog(v,i),Tv_cloud(v,i),Avv(v,i)]=func_Tv(fv_fog(v),fv_cloud(v),fv_local(v),xv,yv,zv);
        end
    end
    for i = 1:I
        Ftheta(i)=max(Tv(:,i));
    end
    
    for i = 1:I
        %For firework O(l)i, run Algorithm 2 and calculate its fitness value 
        vidx=[];
        for vv=1:V
            tmps = Theta(v,:,i,l);
            if length(find(tmps==0))==2
               vidx=[vidx,vv];  
            end
        end
        if isempty(vidx)==1
           vidx=1; 
        end
        delta  = 1e-5*rand(1,V);
        tpm=[]; 
        for vv = 1:length(vidx)
            tpm(vv) = Avs(vidx(vv),i)/fv_fog(vidx(vv)) + delta(vidx(vv));
        end
        thetax = max(tpm); 
        [ffogx]=func_FCRA(thetax,delta,Avv(:,i),vidx); 
        %------------------------------------------------------------------------------------
        %Calculate sˆ(l)iaccording to (16).
        %Calculate the estimated service delay of V-UE v using (13)
        %Generate sˆ(l)iexplosion sparks from firework O(l)
        Ftmps = [];
        es1   = 0.000001;
        Fmax  = max(Ftheta);
        for ij = 1:I
            Ftmps(ij) = Fmax - Ftheta(ij) ;
        end
        s(l,i)=(M * (Fmax - Ftheta(i)  + es1) / (sum(Ftmps) + es1));%ceil
        
        %For each explosion spark, run Algorithm 2.%Calculate the fitness value of each explosion spark
        [ffogx2(i,:)]=func_FCRA(s(l,i),delta,Avv(:,i),vidx); 
    end    
    %Generate mˆ mutation sparks
    tmpx  = mean(ffogx2,2);
    tmpx2 = mean2(ffogx2);
    m_   = find(tmpx<2*tmpx2);
    %For each mutation spark, run Algorithm 2
    %Calculate the fitness value of each mutation spark
    for i = 1:length(m_);
        ffogx3(i,:)=func_FCRA(s(l,m_(i)),delta,Avv(:,m_(i)),vidx); 
    end
    %The firework, explosion spark or mutation spark with the
    %smallest fitness value is chosen as O(l+1), and the smallest
    %fitness value is denoted by F(l).
    F(l)=min(mean(ffogx3));
end



figure;
plot(smooth(F,4),'b-*');
xlabel('Numbero of iteration');
ylabel('Objective value(s)');
grid on
ylim([6,11]);

 

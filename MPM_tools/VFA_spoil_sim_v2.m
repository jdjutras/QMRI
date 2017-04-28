function [A,B]=VFA_spoil_sim_v2(phi0,N_spins,TR,angles,file_name)
%file_name='Parameters4MPM';
load(file_name)
%phi0=50;
%N_spins=1000;
%angles
%TR=mpm_par.TR1;
%ang1=mpm_par.ang1;
%ang2=mpm_par.ang2;

cRF=[0.7 0.8 0.9 1.0 1.1 1.2 1.3];
T1=[700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800];
T2=[80 80 80 80 80 90 100 110 110 110 110 110];
%T2=[110 110 110 110 110 110 110 110 110 110 110 110];
%T2=[250 250 250 250 250 250 250 250 250 250 250 250];
N_T1=length(T1);
N_cRF=length(cRF);
N_angles=length(angles);

T1prime=zeros(N_T1,N_cRF,'single');
entries=length(T1prime(:));
h=waitbar(0,'Calculating coefficients A and B');
t=0;

for j=1:N_cRF
alpha1=angles(1,1)*cRF(j);
alpha2=angles(1,2)*cRF(j);
if N_angles>=3
alpha3=angles(1,3)*cRF(j);
end
if N_angles==4
alpha4=angles(1,4)*cRF(j);
end


for k=1:N_T1
[Msteady1,~,~,~,~,~]=SPGR_sim_KW(alpha1,phi0,TR,T1(k),T2(k),N_spins);
y1=Msteady1/sind(alpha1);
x1=Msteady1/tand(alpha1);
[Msteady2,~,~,~,~,~]=SPGR_sim_KW(alpha2,phi0,TR,T1(k),T2(k),N_spins);
y2=Msteady2/sind(alpha2);
x2=Msteady2/tand(alpha2);

if N_angles==2
E1=(y2-y1)/(x2-x1);
T1prime(k,j)=-TR./log(E1);
elseif N_angles==3
[Msteady3,~,~,~,~,~]=SPGR_sim_KW(alpha3,phi0,TR,T1(k),T2(k),N_spins);
y3=Msteady3/sind(alpha3);
x3=Msteady3/tand(alpha3);
x=[x1 x2 x3];
y=[y1 y2 y3];
p=polyfit(x,y,1);
T1prime(k,j)=-TR./log(p(1));
elseif N_angles==4
[Msteady3,~,~,~,~,~]=SPGR_sim_KW(alpha3,phi0,TR,T1(k),T2(k),N_spins);
y3=Msteady3/sind(alpha3);
x3=Msteady3/tand(alpha3);
[Msteady4,~,~,~,~,~]=SPGR_sim_KW(alpha4,phi0,TR,T1(k),T2(k),N_spins);
y4=Msteady4/sind(alpha4);
x4=Msteady4/tand(alpha4);
x=[x1 x2 x3 x4];
y=[y1 y2 y3 y4];
p=polyfit(x,y,1);
T1prime(k,j)=-TR./log(p(1));
end        
t=t+1;
waitbar(t/entries);
end
end
close(h)

ABcoeff=zeros(N_cRF,2);
for k=1:N_cRF
ABcoeff(k,:)=polyfit(T1prime(:,k)',T1,1);
end

B=polyfit(cRF,ABcoeff(:,1)',2);
A=polyfit(cRF,ABcoeff(:,2)',2);
save(file_name,'A','B','-append');
clearvars -except A B
end
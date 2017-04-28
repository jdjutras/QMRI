
function [M_steady,Mxp,Myp,Mzp,Mxm,Mym,Mzm]=SPGR_sim_KW(flip_ang,phi0,TR,T1,T2,Nspin)
%#codegen
M0=100;
Nrf=max(round(6*T1/TR),150);
alpha=flip_ang*pi/180;
d_phi=phi0*pi/180;

Mp=zeros(3,Nrf,Nspin,'double');
Mm=zeros(3,Nrf,Nspin,'double');
Mm(3,1,:)=M0;
psi_vec=(0:2*pi/(Nspin-1):2*pi);
psi_vec3=permute(psi_vec,[1 3 2]);

n=(1:1:Nrf);
phi_vec=0.5*d_phi*(n.^2+n+2*ones(1,Nrf));
phi_vec3=permute(phi_vec,[1 3 2]);

coeff_mat1=[cos(alpha)+(1-cos(alpha)).*cos(phi_vec3).^2 (1-cos(alpha)).*sin(phi_vec3).*cos(phi_vec3) -sin(alpha).*sin(phi_vec3);
                (1-cos(alpha)).*sin(phi_vec3).*cos(phi_vec3) cos(alpha)+(1-cos(alpha)).*sin(phi_vec3).^2 sin(alpha).*cos(phi_vec3);
                sin(alpha).*sin(phi_vec3) -sin(alpha).*cos(phi_vec3) cos(alpha).*ones(1,1,Nrf)]; %3x3xNrf
            
coeff_mat2=[exp(-TR/T2).*cos(psi_vec3) exp(-TR/T2).*sin(psi_vec3) zeros(1,1,Nspin);
                   -exp(-TR/T2).*sin(psi_vec3) exp(-TR/T2).*cos(psi_vec3) zeros(1,1,Nspin);
                   zeros(1,1,Nspin) zeros(1,1,Nspin) exp(-TR/T1).*ones(1,1,Nspin)]; %3x3xNspin
    
for j=1:Nrf-1        
    Mp(:,j,:)=repmat(coeff_mat1(:,1,j),[1 1 Nspin]).*repmat(Mm(1,j,:),[3 1 1])...
        + repmat(coeff_mat1(:,2,j),[1 1 Nspin]).*repmat(Mm(2,j,:),[3 1 1])...
        + repmat(coeff_mat1(:,3,j),[1 1 Nspin]).*repmat(Mm(3,j,:),[3 1 1]);
    Mm(:,j+1,:)=coeff_mat2(:,1,:).*repmat(Mp(1,j,:),[3 1 1])...
        + coeff_mat2(:,2,:).*repmat(Mp(2,j,:),[3 1 1])...
        + coeff_mat2(:,3,:).*repmat(Mp(3,j,:),[3 1 1]);
    Mm(3,j+1,:)=Mm(3,j+1,:) + ones(1,1,Nspin).*M0*(1-exp(-TR/T1));
    %Mm(3,j+1,:)=Mm(3,j+1,:) + ones(1,1,Nspin).*M0*(TR/T1);
end

mxp=squeeze(Mp(1,:,:));
myp=squeeze(Mp(2,:,:));
mzp=squeeze(Mp(3,:,:));
mxm=squeeze(Mm(1,:,:));
mym=squeeze(Mm(2,:,:));
mzm=squeeze(Mm(3,:,:));
    
Mxp=sum(mxp,2)/Nspin; 
Myp=sum(myp,2)/Nspin; 
Mzp=sum(mzp,2)/Nspin; 
Mxm=sum(mxm,2)/Nspin; 
Mym=sum(mym,2)/Nspin; 
Mzm=sum(mzm,2)/Nspin; 
Mxy=sqrt(Mxp.^2+Myp.^2);
M_steady=mean(Mxy(Nrf-101:Nrf-1));
%S_ideal=M0*sind(alpha)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(alpha));
end
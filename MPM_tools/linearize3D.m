function [R2,sigma_m,m,b]=linearize3D(Y,X)

if size(Y)~=size(X)
    error('Unequal matrix size');
end
N=size(Y,4);
Sx=sum(X,4);
Sxx=sum(X.^2,4);
Sy=sum(Y,4);
Syy=sum(Y.^2,4);
Sxy=sum(X.*Y,4);
Delta=N*Sxx-Sx.^2;

r=(N*Sxy-Sx.*Sy)./sqrt((N*Sxx-Sx.^2).*(N*Syy-Sy.^2));
R2=r.^2;
m=(N*Sxy-Sx.*Sy)./Delta;
b=(Sxx.*Sy-Sx.*Sxy)./Delta;

summation=zeros(size(X),'single');
for i=1:N
summation(:,:,:,i)=((Y(:,:,:,i)-b-m.*X(:,:,:,i)).^2)/(N-2);
end
sigma_y=sqrt(sum(summation,4));
sigma_m=sigma_y.*sqrt(N/Delta);
end
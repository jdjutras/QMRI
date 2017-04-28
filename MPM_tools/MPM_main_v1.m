%MPM_Reg_vX (c) Jean-David Jutras
%This program is designed to reconstruct
%quantitative maps following image registration
%in 3D Slicer. The datasets to input must be
%mha or nii/nii.gz format. Up to 4 flip angles are allowed for the VFA T1
%curve-fit (DESPOT1), and up to 2 flip angles by 4 phase offsets, 
%or 3 flip angles by 3 phase offsets for VFA T2(DESPOT2) See the user guide.
clear
close all
disp('Enter bSSFP datasets: ang1-phi1, ang2-phi1, etc');
disp('Also enter phase offsets in this order: 2 phases: phi1,phi2,');
disp('3 phases: phi1,phi2,phi3, and 4 phases: phi1,phi3,phi2,phi4'); 

%calibration='M0'; 'M0' if you only want to normalize the final M0 map
calibration='M0&T1';
New_T1M0_mask='no'; %'yes' if you want a better mask to estimate the Psi0 bias field

%Define the input parameters and their options
prompt={'Parameters file','Image mask','Bias norm mask','FLASH dataset1','FLASH dataset2',...
    'FLASH dataset3 or MT-FLASH dataset1','FLASH dataset4 or MT-FLASH dataset2',...
    'T_{2}^{*} map1','T_{2}^{*} map2','T_{2}^{*} map3','T_{2}^{*} map4',...
    'bSSFP dataset1','bSSFP dataset2','bSSFP dataset3','bSSFP dataset4','bSSFP dataset5',...
    'bSSFP dataset6','bSSFP dataset7','bSSFP dataset8','bSSFP dataset9',...
    'RF correction type: [BiasField / B1map / no]','c_{RF}^{+} or \Psi_{1} file','c_{RF}^{-} or \Psi_{0} file',...
    'BiasField calibration \Omega: [brain / head / custom]','Bias Field c_{cal0}, T_{1}^{ref} (ms)',...
    'Output maps: [M0appB1 / M0&T1 / M0&T1&T2* / all]','T_{2}^{*} filter: [median / other / no]',...
    'FLASH T_{RF2}, bSSFP T_{RF2}(ms), \phi_{0}(deg)','T_{1}, T_{2}, T_{2}^{*} clip value(ms)','T_{2}^{*} combination: [wave1 / wave2]',...
    '\sigma^{2}_{T2*} map1','\sigma^{2}_{T2*} map2','Calibration mode [manual / auto / WMauto]:'};

dlg_title='Inputs for Multi-Parameter Mapping';
%Set all default parameter file names here
def={'Parameters4MPM','mask_im.mha','mask_im.mha','ME1.mha',...
    'ME2.mha','ME1MT.mha','ME3MT_reg.mha','T2s1.mha','T2s2.mha',...
    'T2s3_4.mha','T2s4_4.mha','bffe1_0.mha','bffe2_0.mha',...
    'bffe1_180.mha','bffe2_180.mha','bffe1_90.mha','bffe2_90.mha',...
    'bffe1_270.mha','bffe2_270.mha','bffe_last.mha','B1map',...
    'cRF_AFI.mha','cRF_AFI.mha','brain','0.97 900','all','median',...
    '0.00 0.55 50','5500 2500 500','wave1','sigmaT2s1.mha','sigmaT2s2.mha','WMauto'};

largest_question_length = size((strvcat(prompt')),2);
lineNo=[1,largest_question_length+1];
AddOpts.Interpreter='tex'; AddOpts.WindowStyle='normal'; AddOpts.Resize='on';

answer=inputdlgcol(prompt,dlg_title,lineNo,def,AddOpts,3);
%Pass on the answers to the respective parameter variables
parameter_filename=answer{1};
load(parameter_filename); 
mask1_im=answer{2}; mask2_im=answer{3};
ME_filename1=answer{4}; ME_filename2=answer{5}; ME_filename3=answer{6};
ME_filename4=answer{7}; T2s_filename1=answer{8}; T2s_filename2=answer{9};
T2s_filename3=answer{10}; T2s_filename4=answer{11};
bFFE_filename1=answer{12}; bFFE_filename2=answer{13};
bFFE_filename3=answer{14}; bFFE_filename4=answer{15};
bFFE_filename5=answer{16}; bFFE_filename6=answer{17};
bFFE_filename7=answer{18}; bFFE_filename8=answer{19};
bFFE_filename9=answer{20}; B1_corr=answer{21}; B1plus_file=answer{22};
B1minus_file=answer{23}; N4ITKnorm=answer{24};
cal_factors=str2num(answer{25});%~0.97 for brain or ~0.87 for head;
calculation=answer{26}; T2s_filt=answer{27};
RF_info=str2num(answer{28}); clip=str2num(answer{29});
T2s_comb=answer{30}; sigma2_filename1=answer{31}; sigma2_filename2=answer{32};
cal_mode=answer{33};
Bias_cal=cal_factors(1,1); T1_ref=cal_factors(1,2);
T_RFE=RF_info(1:2); phi0=RF_info(1,3);

if (strcmpi(calculation,'M0appB1') && strcmpi(B1_corr,'B1map'))
M0appB1='yes'; %If you want to estimate cRF- using a B1 map
calculation='M0&T1&T2*'; %Set the output parameters to this
else
M0appB1='no'; %If you want to estimate cRF- using a B1 map
end    
T1_clip=clip(1,1); %Pass-on the clipping values here
T2_clip=clip(1,2); 
T2s_clip=clip(1,3);
epsilon_cutoff=0; %Cut-off value for dealing with singularity in DESPOT2.
%Can also choose -0.45 based on mean WM T2=50ms.
tic
[~,~,ext] = fileparts(ME_filename1); %Check file extension
if strcmpi(ext,'.mha')
h=mha_read_volume(ME_filename1); clear ME_filename1 %Read-in mha ME1 file
merge1=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.nii.gz')
h=load_nii(ME_filename1); clear ME_filename1  %Read-in nifti ME1 file
merge1=mat2nifti(h.img,'nifti');
end

[~,~,ext] = fileparts(ME_filename2);%Check file extension
if strcmpi(ext,'.mha')
h=mha_read_volume(ME_filename2); clear ME_filename2
merge2=h.pixelData; %2 percent increase
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.nii.gz')
h=load_nii(ME_filename2);  clear ME_filename2      
merge2=mat2nifti(h.img,'nifti');
end

if mpm_par.N_T1FFE_ang>=3 %If 3 or more VFA T1 angles, open the ME3 dataset   
[~,~,ext] = fileparts(ME_filename3); %Check file extension
 if strcmpi(ext,'.mha')    
 h=mha_read_volume(ME_filename3); clear ME_filename3
 merge3=h.pixelData; 
 elseif strcmpi(ext,'.nii') || strcmpi(ext,'.nii.gz')
 h=load_nii(ME_filename3); clear ME_filename3       
 merge3=mat2nifti(h.img,'nifti');
 end
end

if mpm_par.N_T1FFE_ang==4 %If 4 or more VFA T1 angles, open the ME4 dataset   
[~,~,ext] = fileparts(ME_filename4);
if strcmpi(ext,'.mha')        
h=mha_read_volume(ME_filename4); clear ME_filename4
merge4=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.nii.gz')
h=load_nii(ME_filename4);  clear ME_filename4      
merge4=mat2nifti(h.img,'nifti');
end    
end

if strcmpi(mpm_par.MT_scan,'yes') %Open the ME1MT dataset 
[~,~,ext] = fileparts(ME_filename3); %(ME1 with same echoes as ME3MT)
if strcmpi(ext,'.mha')            
h=mha_read_volume(ME_filename3); clear ME_filename3
merge1MT=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
h=load_nii(ME_filename3);  clear ME_filename3
merge1MT=mat2nifti(h.img,'nifti');
end    
[~,~,ext] = fileparts(ME_filename4); %Open the ME3MT dataset 
if strcmpi(ext,'.mha')         
h=mha_read_volume(ME_filename4); clear ME_filename4
merge3MT=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
h=load_nii(ME_filename4);  clear ME_filename4
merge3MT=mat2nifti(h.img,'nifti');
end    
end
x=size(merge1,1); y=size(merge1,2); no_slices=size(merge1,3);
%% Masking using Otsu thresholding
[~,~,ext] = fileparts(mask1_im);
if strcmpi(ext,'.mha')         
h=mha_read_volume(mask1_im);
image_mask=h.pixelData; %clear mask1_im
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
h=load_nii(mask1_im);        
image_mask=mat2nifti(h.img,'nifti');
end    
gray_im=mat2gray(image_mask); clear image_mask
level=graythresh(gray_im);
object=zeros(size(gray_im),'int8');
    for i=1:mpm_par.no_slices
    object(:,:,i)=imbinarize(squeeze(gray_im(:,:,i)),0.3*level); %change from 0.5 to 0.2
    end
mask=logical(1-object); clear gray_im level
object=logical(object);
voxel_size=mpm_par.fov(1,1)*mpm_par.fov(1,2)*mpm_par.fov(1,3)/x/y/mpm_par.no_slices;
Volume=sum(sum(sum(object)))*voxel_size/100^3; %volume in Litres

if strcmpi(N4ITKnorm,'head') %if the calibration factor is unity, use the following
    %Bias_cal=0.941-0.0122*Volume; %based on SDAM
    Bias_cal=0.9142-0.0123*Volume;  %based on AFI
brain_mask=mask;
brain_object=object;
elseif strcmpi(N4ITKnorm,'brain')
%hist_corr1=0.97;    
[~,~,ext] = fileparts(mask2_im);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(mask2_im); clear mask2_im
    brain_mask=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(mask2_im);  clear mask2_im      
    brain_mask=mat2nifti(h.img,'nifti');
    end
gray_im=mat2gray(brain_mask); clear brain_mask
level=graythresh(gray_im);
brain_obj=zeros(size(gray_im),'int8');
    for i=1:mpm_par.no_slices
    brain_obj(:,:,i)=imbinarize(squeeze(gray_im(:,:,i)),0.5*level);
    end
brain_mask=logical(1-brain_obj); clear gray_im level 
brain_object=logical(brain_obj); clear brain_obj  
else
brain_mask=mask;
brain_object=object;
end

if strcmpi(B1_corr,'BiasField') 
    [~,~,ext] = fileparts(B1plus_file);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(B1plus_file); clear B1plus_file %First iteration
    c=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(B1plus_file); clear B1plus_file       
    c=mat2nifti(h.img,'nifti');
    end 
    cRF=zeros(size(c),'single');
    c_orig=sqrt(c);
    
    cRF(brain_object)=sqrt(c(brain_object)); %dump c values into the empty matrix (only where mask=1)
    vec=cRF(:);
    ind=vec<=0.001;
    vec(ind)=NaN;
    vec(isnan(vec))=[];
    mean_value=mean(vec);
    c=c_orig/mean_value*Bias_cal; clear cRF c_orig
    
    cRF2=zeros(size(c),'single');   
    cRF2(brain_object)=c(brain_object); 
    vec=cRF2(:);
    ind=vec<=0.001;
    vec(ind)=NaN;
    vec(isnan(vec))=[];
    mean_value=mean(vec);
    disp('<cRF+> after normalization:'); disp(mean_value); %display mean value after first iteration
    
    [~,~,ext] = fileparts(B1minus_file);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(B1minus_file); clear B1minus_file
    cm=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(B1minus_file);  clear B1minus_file       
    cm=mat2nifti(h.img,'nifti');
    end
    
    cRFm=zeros(size(cm),'single');
    cRFm(brain_object)=cm(brain_object);
    vec=cRFm(:);
    ind=vec==0;
    vec(ind)=NaN;
    vec(isnan(vec))=[];
    mean_value=mean(sqrt(vec));
    %cRF2=sqrt(cm/mean_value)*Bias_cal; clear cRFm
    %c_minus=(cRF2.^2)./c; clear cRF2
    c_minus=Bias_cal^2*cm./(c*mean_value^2);
    cal_old=Bias_cal;
    
    cRF2=zeros(size(c),'single');   
    cRF2(brain_object)=c_minus(brain_object); 
    vec=cRF2(:);
    ind=vec<=0.001;
    vec(ind)=NaN;
    vec(isnan(vec))=[];
    mean_value=mean(vec);
    disp('<cRF-> after normalization:'); disp(mean_value); %display mean value after first iteration

elseif strcmpi(B1_corr,'no')
    c=ones(size(merge1),'single');
    c_minus=c;
    
elseif strcmpi(B1_corr,'B1map')
    [~,~,ext] = fileparts(B1plus_file);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(B1plus_file); clear B1plus_file 
    c=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(B1plus_file); clear B1plus_file       
    c=mat2nifti(h.img,'nifti');
    end
    [~,~,ext] = fileparts(B1minus_file);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(B1minus_file); clear B1minus_file
    c_minus=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(B1minus_file);  clear B1minus_file       
    c_minus=mat2nifti(h.img,'nifti');
    end
    
    cRF2=zeros(size(c),'single');   
    cRF2(brain_object)=c_minus(brain_object); 
    vec=cRF2(:);
    ind=vec<=0.001;
    vec(ind)=NaN;
    vec(isnan(vec))=[];
    mean_value=mean(vec);
    c_minus=c_minus/mean_value; clear mean_value vec cRF2
end
%% DESPOT1 Calculation
if mpm_par.N_T1FFE_ang==2 %This section is just to look up the TR.
    TR=(mpm_par.TR1+mpm_par.TR2)/2;
    TE=(mpm_par.TE1+mpm_par.TE2)/2;
elseif mpm_par.N_T1FFE_ang==3
    TR=(mpm_par.TR1+mpm_par.TR2+mpm_par.TR3)/3;
    TE=(mpm_par.TE1+mpm_par.TE2+mpm_par.TE3)/3;
else
    TR=(mpm_par.TR1+mpm_par.TR2+mpm_par.TR3+mpm_par.TR4)/4;
    TE=(mpm_par.TE1+mpm_par.TE2+mpm_par.TE3+mpm_par.TE4)/4;
end
if size(c,1)~=size(merge1,1) %Check the dimensions of the c matrix   
    c_new=resize(c,[x y mpm_par.no_slices]);
    c=c_new; clear c_new
end

Y=zeros(x,y,mpm_par.no_slices,mpm_par.N_T1FFE_ang,'single'); X1=zeros(x,y,mpm_par.no_slices,mpm_par.N_T1FFE_ang,'single'); 
Y(:,:,:,1)=merge1./sind(mpm_par.ang1.*c);
X1(:,:,:,1)=merge1./tand(mpm_par.ang1.*c);
if T_RFE(1,1)<=0.1   
Y(:,:,:,2)=merge2./sind(mpm_par.ang2.*c);
X1(:,:,:,2)=merge2./tand(mpm_par.ang2.*c);
else
    T2star=50; %use a global T2*
    unity=ones(size(c));
    tau=T_RFE(1,1)*unity./(2*T2star);
phi=sqrt((mpm_par.ang2.*c*pi/180).^2-tau.^2);
b=exp(-tau).*mpm_par.ang2.*c.*pi./180./phi.*sin(phi);
a=exp(-tau).*(cos(phi)+tau./phi.*sin(phi));
Y(:,:,:,2)=merge2./b;
X1(:,:,:,2)=merge2.*a./b;
end

if mpm_par.N_T1FFE_ang>=3
Y(:,:,:,3)=merge3./sind(mpm_par.ang3.*c);
X1(:,:,:,3)=merge3./tand(mpm_par.ang3.*c);
end
if mpm_par.N_T1FFE_ang==4
Y(:,:,:,4)=merge4./sind(mpm_par.ang4.*c);
X1(:,:,:,4)=merge4./tand(mpm_par.ang4.*c);
end

[~,~,slope,int]=linearize3D(Y,X1); clear Y X1

T1=-TR./real(log(slope)); ind=T1>T1_clip; T1(ind)=T1_clip; ind=T1<0; T1(ind)=T1_clip;
T1(mask)=1; 
    
if (~strcmpi(calibration,'M0') && strcmpi(B1_corr,'BiasField'))
T1_uncal=T1;    

if strcmpi(mask1_im,'mask_im.mha') %Check if we need to exclude the neck
disp('Select a cutoff location below the brain.');
figure(1); imshow(T1(:,:,87),[0 3000]);
[~,yi]=ginput(1); clf;
last_slice=round(yi);
else
last_slice=size(T1,2); %If image is skull-stripped, don't worry about cutoff
end

T1_vec=(0:1:5000);
T1_crop=T1(1:last_slice,:,:);
mask_crop=mask(1:last_slice,:,:);

T1_hist=hist(T1_crop(:),T1_vec);
Npixels=size(T1_crop,1)*size(T1_crop,2)*size(T1_crop,3)-sum(sum(sum(mask_crop)));
T1_hist_norm=T1_hist/Npixels;
T1_hist_norm(1,1:10)=0;
T1_hist_norm(1,4997:end)=0;
T1_max=max(T1_hist_norm);
if strcmpi(cal_mode,'manual')
figure(1); plot(T1_vec,T1_hist_norm); ylim([0 T1_max]); xlim([0 2500]); 
    title('Normalized T1 histogram');
disp('Select two points on the T1 histogram corresponding to the left/right side of the WM peak');
[xi,~]=ginput(2);
T1_meas=(xi(2)+xi(1))/2; clf;
else
    if strcmpi(mask1_im,'mask_im_strip.nii.gz')   
    options1=fitoptions('gauss5');
    options1.StartPoint=[0.7e-3 890 60 0.7e-3 1550 275 0.5e-3 990 260 0.15e-3 1800 600 0.03e-3 3000 800];
    options1.Lower=[0.7e-3 800 40 0.6e-3 1300 220 0.4e-3 940 200 0.1e-3 1600 400 0.02e-3 2500 600];
    options1.Upper=[2.0e-3 1000 90 1.2e-3 1700 330 0.8e-3 1060 300 0.2e-3 2000 800 0.04e-3 3500 1000];
    option1.Normalize='on';
    option1.Robust='on';
    f1=fit(T1_vec',T1_hist_norm','gauss5',options1);
    disp(f1);
        if strcmpi(cal_mode,'WMauto')
        T1_meas=f1.b1;
        elseif strcmpi(cal_mode,'auto')
        T1_meas=(f1.b1+f1.b2)/2;
        end    
    else 
    options1=fitoptions('gauss7');
    options1.StartPoint=[0.7e-3 890 60 0.7e-3 1550 275 0.6e-4 475 85 0.5e-3 990 260 0.3e-3 590 150 0.15e-3 1800 600 0.03e-3 3000 800];
    options1.Lower=[0.4e-3 750 40 0.5e-3 1300 200 0 350 65 0.4e-3 900 230 0.2e-3 540 130 0.1e-3 1700 400 0.02e-3 2500 600];
    options1.Upper=[1.7e-3 1100 100 1.5e-3 1700 400 0.7e-3 550 105 0.6e-3 1100 290 0.4e-3 650 170 0.2e-3 2000 800 0.04e-3 3500 1000];
    option1.Normalize='on';
    option1.Robust='on';
    f1=fit(T1_vec',T1_hist_norm','gauss7',options1);
    disp(f1);
        if strcmpi(cal_mode,'WMauto')
        T1_meas=f1.b1;
        elseif strcmpi(cal_mode,'auto')
        T1_meas=(f1.b1+f1.b2)/2;
        end
    end
end    
cal_new=cal_old-(T1_ref-T1_meas)/(2*T1_ref);
    cRF=zeros(size(c),'single');
    cRF(brain_object)=c(brain_object); %dump c values into the empty matrix (only where mask=1)
    vec=cRF(:);
    ind=vec<=0.001;
    vec(ind)=NaN;
    vec(isnan(vec))=[];
    mean_value=mean(vec);
    c_new=c/mean_value*cal_new; clear cRF 
disp(['Assuming GM/WM reference T1=',num2str(T1_ref),'ms, cal_new=',num2str(cal_new),', cal_old=',num2str(cal_old)]);

Y=zeros(x,y,mpm_par.no_slices,mpm_par.N_T1FFE_ang,'single'); X1=zeros(x,y,mpm_par.no_slices,mpm_par.N_T1FFE_ang,'single'); 
Y(:,:,:,1)=merge1./sind(mpm_par.ang1.*c_new);
X1(:,:,:,1)=merge1./tand(mpm_par.ang1.*c_new);
if T_RFE(1,1)<=0.1   
Y(:,:,:,2)=merge2./sind(mpm_par.ang2.*c_new);
X1(:,:,:,2)=merge2./tand(mpm_par.ang2.*c_new);
else
    T2star=50; %use a global T2*
    unity=ones(size(c_new));
    tau=T_RFE(1,1)*unity./(2*T2star);
phi=sqrt((mpm_par.ang2.*c_new*pi/180).^2-tau.^2);
b=exp(-tau).*mpm_par.ang2.*c_new.*pi./180./phi.*sin(phi);
a=exp(-tau).*(cos(phi)+tau./phi.*sin(phi));
Y(:,:,:,2)=merge2./b;
X1(:,:,:,2)=merge2.*a./b;
end
if mpm_par.N_T1FFE_ang>=3
Y(:,:,:,3)=merge3./sind(mpm_par.ang3.*c_new);
X1(:,:,:,3)=merge3./tand(mpm_par.ang3.*c_new);
end
if mpm_par.N_T1FFE_ang==4
Y(:,:,:,4)=merge4./sind(mpm_par.ang4.*c_new);
X1(:,:,:,4)=merge4./tand(mpm_par.ang4.*c_new);
end

[~,~,slope,int]=linearize3D(Y,X1); clear Y X1
T1=-TR./real(log(slope)); ind=T1>T1_clip; T1(ind)=T1_clip; ind=T1<0; T1(ind)=T1_clip;
T1(mask)=1;
c=c_new;
end

    if (exist('A','var')==1 && exist('B','var') && strcmpi(B1_corr,'B1map'))
    T1corr=A(1)*c.^2+A(2)*c+A(3)+T1.*(B(1)*c.^2+B(2).*c+B(3)); %Correction for 
    T1_uncorr=T1; %non-ideal RF spoiling based on Preibisch and Deichmann's paper
    T1=T1corr; clear T1corr
    elseif phi0>0 && strcmpi(B1_corr,'B1map')
        %If the A,B coeff cannot be found but phi0>0, calculate them here
        if mpm_par.N_T1FFE_ang==2
        angles=[mpm_par.ang1 mpm_par.ang2];
        [A,B]=VFA_spoil_sim_v2(phi0,1000,TR,angles,parameter_filename);
        elseif mpm_par.N_T1FFE_ang==3
         angles=[mpm_par.ang1 mpm_par.ang2 mpm_par.ang3];
        [A,B]=VFA_spoil_sim_v2(phi0,1000,TR,angles,parameter_filename);         
        elseif mpm_par.N_T1FFE_ang==4
         angles=[mpm_par.ang1 mpm_par.ang2 mpm_par.ang3 mpm_par.ang4];
         [A,B]=VFA_spoil_sim_v2(phi0,1000,TR,angles,parameter_filename);
        end
    T1corr=A(1)*c.^2+A(2)*c+A(3)+T1.*(B(1)*c.^2+B(2).*c+B(3)); %Correction for 
    T1_uncorr=T1; %non-ideal RF spoiling based on Preibisch and Deichmann's paper
    T1=T1corr; clear T1corr
    end    
xp=mpm_par.fov(1,1)/x; yp=mpm_par.fov(1,2)/y; zp=mpm_par.fov(1,3)/mpm_par.no_slices;
if strcmpi(B1_corr,'no')
    T1app=T1; %Create an apparent T1 map in order to obtain the cRF+ field in N4ITK
    T1app(mask)=1;
    
    R1app=1./T1app;
    ind=R1app<=(1/5500);
    R1app(ind)=(1/5500);
    ind=R1app>=0.01;
    R1app(ind)=0.01;
    R1app(mask)=0;
    ind=isnan(R1app); R1app(ind)=0; R1app=R1app*1000;

v=[10 450 950 1500 3000];
options=struct('maxit',3,'epsilon',1e-5,'alpha',0.5,'sigma',0.5,'p',2);
[~,segments]=BCFCM3D(T1app,v,options);
CSF_mask=segments(:,:,:,5)>=0.1;
T1_mask=segments(:,:,:,1)<=0.1;
T1_mask(CSF_mask)=0; %clear segments v options

metaData=struct('CompressedData','False','ObjectType','Image','NumOfDimensions',3,...
    'BinaryData','True','BinaryDataByteOrderMSB','False','ImageAxesOrientation',[0 0 -1 0 1 0 -1 0 0],...
    'Origin',[mpm_par.fov(1,1)/2 -mpm_par.fov(1,2)/2 mpm_par.fov(1,3)/2],'CenterOfRotation',...
    [mpm_par.fov(1,1)/2 -mpm_par.fov(1,2)/2 mpm_par.fov(1,3)/2],'AnatomicalOrientation','LPS',...
    'Spacing',[xp yp zp],'Size',[x y mpm_par.no_slices],'ElementType','MET_FLOAT','ElementDataFile','Local');
pixelData=single(T1app);
pixelWithMetaStruct=struct('metaData',metaData,'pixelData',pixelData);
mha_write_volume('T1app.mha',pixelWithMetaStruct,true);

pixelData=single(R1app);
pixelWithMetaStruct=struct('metaData',metaData,'pixelData',pixelData);
mha_write_volume('R1app.mha',pixelWithMetaStruct,true);

nii1=mat2nifti(R1app,'mat');
R1_nii=make_nii(nii1,[zp xp yp],[mpm_par.fov(1,3)/2 mpm_par.fov(1,1)/2 mpm_par.fov(1,2)/2],16);
save_nii(R1_nii,'R1app.nii'); clear nii1 R1_nii R1app

pixelData=single(T1_mask);
pixelWithMetaStruct=struct('metaData',metaData,'pixelData',pixelData);
mha_write_volume('T1_mask.mha',pixelWithMetaStruct,true);
end
%% M0 and T2* Calculation
if strcmpi(calculation,'M0&T1&T2*')||strcmpi(calculation,'M0&T1')||strcmpi(calculation,'all')
    
    [~,~,ext] = fileparts(T2s_filename1);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(T2s_filename1); clear T2s_filename1
    T2s1=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(T2s_filename1); clear T2s_filename1       
    T2s1=mat2nifti(h.img,'nifti');
    end     
    ind=T2s1<1; T2s1(ind)=1;
    
    [~,~,ext] = fileparts(T2s_filename2);
    if strcmpi(ext,'.mha')         
    h=mha_read_volume(T2s_filename2); clear T2s_filename2
    T2s2=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(T2s_filename2); clear T2s_filename2       
    T2s2=mat2nifti(h.img,'nifti');
    end     
    ind=T2s2<1; T2s2(ind)=1;
    
    if mpm_par.N_T1FFE_ang>=3        
    [~,~,ext] = fileparts(T2s_filename3);
    if strcmpi(ext,'.mha') 
    h=mha_read_volume(T2s_filename3); clear T2s_filename3
    T2s3=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(T2s_filename3); clear T2s_filename3       
    T2s3=mat2nifti(h.img,'nifti');
    end     
    ind=T2s3<1; T2s3(ind)=1;
    
    end
    if mpm_par.N_T1FFE_ang==4
    [~,~,ext] = fileparts(T2s_filename4);
    if strcmpi(ext,'.mha')     
    h=mha_read_volume(T2s_filename4); clear T2s_filename4
    T2s4=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(T2s_filename4); clear T2s_filename4       
    T2s4=mat2nifti(h.img,'nifti');
    end
    ind=T2s4<1; T2s4(ind)=1;    
    end
%Combine T2* maps into a final map by weighted average
    unity=ones(size(T1),'single');
    E1=exp(-TR./T1);
if strcmpi(T2s_comb,'wave1')
    w1=sind(mpm_par.ang1.*c).^2./(unity-E1.*cosd(mpm_par.ang1.*c)).^2;
    w2=sind(mpm_par.ang2.*c).^2./(unity-E1.*cosd(mpm_par.ang2.*c)).^2;
    T2s=(w1.*T2s1+w2.*T2s2)./(w1+w2);
    if mpm_par.N_T1FFE_ang>=3
    w3=sind(mpm_par.ang3.*c).^2./(unity-E1.*cosd(mpm_par.ang3.*c)).^2;    
    T2s=(w1.*T2s1+w2.*T2s2+w3.*T2s3)./(w1+w2+w3);
    end
    if mpm_par.N_T1FFE_ang==4
    w4=sind(mpm_par.ang4.*c).^2./(unity-E1.*cosd(mpm_par.ang4.*c)).^2;    
    T2s=(w1.*T2s1+w2.*T2s2+w3.*T2s3+w4.*T2s4)./(w1+w2+w3+w4);
    end    
    clear w1 w2 w3 w4 E1
else    
[~,~,ext] = fileparts(sigma2_filename1);     
if strcmpi(ext,'.mha')           
h=mha_read_volume(sigma2_filename1); clear sigma2_filename1
sigmaT2s1=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
h=load_nii(sigma2_filename1);  clear sigma2_filename1      
sigmaT2s1=mat2nifti(h.img,'nifti');
end
ind=sigmaT2s1<=0; sigmaT2s1(ind)=4*T2s_clip;

[~,~,ext] = fileparts(sigma2_filename2);     
if strcmpi(ext,'.mha')           
h=mha_read_volume(sigma2_filename2); clear sigma2_filename2 
sigmaT2s2=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
h=load_nii(sigma2_filename2);  clear sigma2_filename2       
sigmaT2s2=mat2nifti(h.img,'nifti');
end

ind=sigmaT2s2<=0; sigmaT2s2(ind)=4*T2s_clip;

w1=sigmaT2s1.^-2; w2=sigmaT2s2.^-2;
T2s=(T2s1.*w1+T2s2.*w2)./(w1+w2);
clear indA indB indAB indC indD w1 w2
end

if strcmpi(T2s_filt,'median')
T2s_filt=zeros(size(T2s),'single'); %Filter the T2* map
for i=1:no_slices
    T2s_filt(:,:,i)=medfilt2(T2s(:,:,i),[3 3]); %Use 2D median filter
end
elseif strcmpi(T2s_filt,'other')
    T2s_filt=T2s;
end
%T2s=T2s_filt;
%Clip T2s values that are negative
ind1=T2s<=1; ind2=T2s>=T2s_clip;
T2s(ind1)=T2s_filt(ind1); T2s(ind2)=T2s_filt(ind2);
ind=isnan(T2s); T2s(ind)=1;

% Calculate the SNR gain
no_echoes=length(TE);
exp_sum=zeros(x,y,mpm_par.no_slices,mpm_par.last_echo(1,1),'single');
    
for i=1:mpm_par.last_echo(1,1)
    exp_sum(:,:,:,i)=exp(-2*TE(i)./T2s_filt);
end
if no_echoes>=6 %Only compensate for SNR gain if there are 6 echoes or more
    G=sqrt(sum(exp_sum,4));
ind=G<1; G(ind)=1; ind=isnan(G); G(ind)=1; %Constrain the SNR gain
elseif no_echoes==1
    G=sqrt(squeeze(exp_sum));
ind=G<0.88; G(ind)=0.88; ind=isnan(G); G(ind)=0.88; %Constrain the SNR gain
else
    G=ones(size(merge1),'single');
end

if strcmpi(M0appB1,'no')
PD=real(int./((unity-slope).*c_minus.*G)); %Compute PD map with B1- corr.
elseif (strcmpi(M0appB1,'yes') && strcmpi(B1_corr,'B1map'))
PD=real(int./((unity-slope).*G)); %Compute PD map without B1- corr.
else
PD=real(int./((unity-slope).*c_minus.*G)); %Compute PD map with B1- corr.    
end    
clear X1 Y int ind1 ind2 ind3 exp_sum G

T2s(mask)=1; PD(mask)=1; %Mask out the air cavities
%Scale PD map
vec=PD(:); ind=vec<=1; vec(ind)=NaN; vec(isnan(vec))=[]; vec(isinf(vec))=[]; 
PD_mean=mean(vec);

if (strcmpi(B1_corr,'no')||strcmpi(M0appB1,'yes'))
PD_max=3*PD_mean;
PD_new=PD/PD_max; 
ind=PD_new<0;
PD_new(ind)=0;
ind=PD_new>1;
PD_new(ind)=1;
ind=isnan(PD_new);
PD_new(ind)=0;
M0app=4096*PD_new;
R0app=1./M0app; ind=R0app>0.002; R0app(ind)=0; R0app=R0app*1000;

metaData=struct('CompressedData','False','ObjectType','Image','NumOfDimensions',3,...
    'BinaryData','True','BinaryDataByteOrderMSB','False','ImageAxesOrientation',[0 0 -1 0 1 0 -1 0 0],...
    'Origin',[mpm_par.fov(1,1)/2 -mpm_par.fov(1,2)/2 mpm_par.fov(1,3)/2],'CenterOfRotation',...
    [mpm_par.fov(1,1)/2 -mpm_par.fov(1,2)/2 mpm_par.fov(1,3)/2],'AnatomicalOrientation','LPS',...
    'Spacing',[xp yp zp],'Size',[x y mpm_par.no_slices],'ElementType','MET_FLOAT','ElementDataFile','Local');
pixelData=M0app;
pixelWithMetaStruct=struct('metaData',metaData,'pixelData',pixelData);

if strcmpi(M0appB1,'no')
mha_write_volume('M0app.mha',pixelWithMetaStruct,true);
else
mha_write_volume('M0appB1.mha',pixelWithMetaStruct,true);
end

pixelData=single(R0app);
pixelWithMetaStruct=struct('metaData',metaData,'pixelData',pixelData);

if strcmpi(M0appB1,'no')
mha_write_volume('R0app.mha',pixelWithMetaStruct,true);
nii0=mat2nifti(R0app,'mat');
R0_nii=make_nii(nii0,[zp xp yp],[mpm_par.fov(1,3)/2 mpm_par.fov(1,1)/2 mpm_par.fov(1,2)/2],16);
save_nii(R0_nii,'R0app.nii'); clear nii0 R0_nii %R0app
else
mha_write_volume('R0appB1.mha',pixelWithMetaStruct,true);
nii0=mat2nifti(R0app,'mat');
R0_nii=make_nii(nii0,[zp xp yp],[mpm_par.fov(1,3)/2 mpm_par.fov(1,1)/2 mpm_par.fov(1,2)/2],16);
save_nii(R0_nii,'R0appB1.nii'); clear nii0 R0_nii %R0app
end
    

else
PD_max=2*PD_mean;
PD_new=PD/PD_max;
ind=PD_new>=PD_max;
PD_new(ind)=PD_max;

ind=PD_new<0;
PD_new(ind)=0;
ind=PD_new>1;
PD_new(ind)=1;
ind=isnan(PD_new);
PD_new(ind)=0;
M0=PD_new*145.80; %A scaling factor that will bring PD close to its normalized value
end
clear vec ind PD_mean PD_max PD_new PD          
end
           
if (strcmpi(mpm_par.MT_scan,'yes') && strcmpi(calculation,'all'))
MTR=(merge1MT-merge3MT)./merge1MT;
MTR(mask)=0;
ang1=mpm_par.ang1*pi/180;

MTsat=(TR./T1+unity.*(ang1.^2)/2)./(unity./MTR-unity);
MTsat_corr=MTsat.*(unity-0.4*c)/0.6;

MTsat=MTsat_corr; clear MTsat_corr
MTsat(mask)=0;
end
%% Open bFFE datasets
if strcmpi(calculation,'T1&T2')||strcmpi(calculation,'all')
if mpm_par.N_offsets==2  
    if mpm_par.N_bFFE_ang==2        
    [~,~,ext] = fileparts(bFFE_filename1);     
    if strcmpi(ext,'.mha')               
    h=mha_read_volume(bFFE_filename1); clear bFFE_filename1
    bffe1_0=h.pixelData;
    h=mha_read_volume(bFFE_filename2); clear bFFE_filename2
    bffe2_0=h.pixelData;
    h=mha_read_volume(bFFE_filename3); clear bFFE_filename3
    bffe1_180=h.pixelData;
    h=mha_read_volume(bFFE_filename4); clear bFFE_filename4
    bffe2_180=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(bFFE_filename1);        
    bffe1_0=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename2);        
    bffe2_0=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename3);        
    bffe1_180=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename4);        
    bffe2_180=mat2nifti(h.img,'nifti');
    end

    elseif mpm_par.N_bFFE_ang==3        
    [~,~,ext] = fileparts(bFFE_filename1);     
    if strcmpi(ext,'.mha')                  
    h=mha_read_volume(bFFE_filename1); clear bFFE_filename1
    bffe1_0=h.pixelData;
    h=mha_read_volume(bFFE_filename2); clear bFFE_filename2
    bffe2_0=h.pixelData;
    h=mha_read_volume(bFFE_filename3); clear bFFE_filename3
    bffe3_0=h.pixelData;
    h=mha_read_volume(bFFE_filename4); clear bFFE_filename4
    bffe1_180=h.pixelData;
    h=mha_read_volume(bFFE_filename5); clear bFFE_filename5
    bffe2_180=h.pixelData;
    h=mha_read_volume(bFFE_filename6); clear bFFE_filename6
    bffe3_180=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(bFFE_filename1);        
    bffe1_0=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename2);        
    bffe2_0=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename3);        
    bffe3_0=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename4);        
    bffe1_180=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename5);        
    bffe2_180=mat2nifti(h.img,'nifti');
    h=load_nii(bFFE_filename6);        
    bffe3_180=mat2nifti(h.img,'nifti');    
    end

    elseif mpm_par.N_bFFE_ang==4
        
    [~,~,ext] = fileparts(bFFE_filename1);     
    if strcmpi(ext,'.mha')               
    h=mha_read_volume(bFFE_filename1); clear bFFE_filename1
    bffe1_0=h.pixelData;
    h=mha_read_volume(bFFE_filename2); clear bFFE_filename2
    bffe2_0=h.pixelData;
    h=mha_read_volume(bFFE_filename3); clear bFFE_filename3
    bffe3_0=h.pixelData;
    h=mha_read_volume(bFFE_filename4); clear bFFE_filename4
    bffe4_0=h.pixelData;
    h=mha_read_volume(bFFE_filename5); clear bFFE_filename5
    bffe1_180=h.pixelData;
    h=mha_read_volume(bFFE_filename6); clear bFFE_filename6
    bffe2_180=h.pixelData;
    h=mha_read_volume(bFFE_filename7); clear bFFE_filename7
    bffe3_180=h.pixelData;
    h=mha_read_volume(bFFE_filename8); clear bFFE_filename8
    bffe4_180=h.pixelData;    
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(bFFE_filename1);        
    bffe1_0=mat2nifti(h.img,'nifti'); clear bFFE_filename1   
    h=load_nii(bFFE_filename2);      
    bffe2_0=mat2nifti(h.img,'nifti'); clear bFFE_filename2
    h=load_nii(bFFE_filename3);         
    bffe3_0=mat2nifti(h.img,'nifti'); clear bFFE_filename3
    h=load_nii(bFFE_filename4);        
    bffe4_0=mat2nifti(h.img,'nifti'); clear bFFE_filename4 
    h=load_nii(bFFE_filename5);       
    bffe1_180=mat2nifti(h.img,'nifti'); clear bFFE_filename5
    h=load_nii(bFFE_filename6);        
    bffe2_180=mat2nifti(h.img,'nifti'); clear bFFE_filename6
    h=load_nii(bFFE_filename7);         
    bffe3_180=mat2nifti(h.img,'nifti'); clear bFFE_filename7     
    h=load_nii(bFFE_filename8);       
    bffe4_180=mat2nifti(h.img,'nifti'); clear bFFE_filename8        
    end
    end
    
im1=zeros(x,y,mpm_par.no_slices,mpm_par.N_offsets,'single');
im2=im1;
im1(:,:,:,1)=bffe1_0;   clear bffe1_0
im1(:,:,:,2)=bffe1_180; clear bffe1_180
im2(:,:,:,1)=bffe2_0;   clear bffe2_0
im2(:,:,:,2)=bffe2_180; clear bffe2_180
if mpm_par.N_bFFE_ang>=3
im3(:,:,:,1)=bffe3_0;   clear bffe3_0
im3(:,:,:,2)=bffe3_180; clear bffe3_180
end
if mpm_par.N_bFFE_ang==4
im4(:,:,:,1)=bffe4_0;   clear bffe4_0
im4(:,:,:,2)=bffe4_180; clear bffe4_180
end

elseif mpm_par.N_offsets==3 
    if mpm_par.N_bFFE_ang==2
    [~,~,ext] = fileparts(bFFE_filename1);     
    if strcmpi(ext,'.mha')                       
    h=mha_read_volume(bFFE_filename1); clear bFFE_filename1
    bffe1_0=h.pixelData;
    h=mha_read_volume(bFFE_filename2); clear bFFE_filename2
    bffe2_0=h.pixelData;
    h=mha_read_volume(bFFE_filename3); clear bFFE_filename3
    bffe1_120=h.pixelData;
    h=mha_read_volume(bFFE_filename4); clear bFFE_filename4
    bffe2_120=h.pixelData;
    h=mha_read_volume(bFFE_filename5); clear bFFE_filename5
    bffe1_240=h.pixelData;
    h=mha_read_volume(bFFE_filename6); clear bFFE_filename6
    bffe2_240=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(bFFE_filename1);        
    bffe1_0=mat2nifti(h.img,'nifti'); clear bFFE_filename1   
    h=load_nii(bFFE_filename2);      
    bffe2_0=mat2nifti(h.img,'nifti'); clear bFFE_filename2
    h=load_nii(bFFE_filename3);         
    bffe1_120=mat2nifti(h.img,'nifti'); clear bFFE_filename3
    h=load_nii(bFFE_filename4);         
    bffe2_120=mat2nifti(h.img,'nifti'); clear bFFE_filename4
    h=load_nii(bFFE_filename5);         
    bffe1_240=mat2nifti(h.img,'nifti'); clear bFFE_filename5
    h=load_nii(bFFE_filename6);         
    bffe2_240=mat2nifti(h.img,'nifti'); clear bFFE_filename6
    end
        
    elseif mpm_par.N_bFFE_ang==3
    [~,~,ext] = fileparts(bFFE_filename1);     
    if strcmpi(ext,'.mha')                               
    h=mha_read_volume(bFFE_filename1); clear bFFE_filename1
    bffe1_0=h.pixelData;
    h=mha_read_volume(bFFE_filename2); clear bFFE_filename2
    bffe2_0=h.pixelData;
    h=mha_read_volume(bFFE_filename3); clear bFFE_filename3
    bffe3_0=h.pixelData;
    h=mha_read_volume(bFFE_filename4); clear bFFE_filename4
    bffe1_120=h.pixelData;
    h=mha_read_volume(bFFE_filename5); clear bFFE_filename5
    bffe2_120=h.pixelData;
    h=mha_read_volume(bFFE_filename6); clear bFFE_filename6
    bffe3_120=h.pixelData;    
    h=mha_read_volume(bFFE_filename7); clear bFFE_filename7
    bffe1_240=h.pixelData;
    h=mha_read_volume(bFFE_filename8); clear bFFE_filename8
    bffe2_240=h.pixelData;
    h=mha_read_volume(bFFE_filename9); clear bFFE_filename9
    bffe3_240=h.pixelData;
    elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
    h=load_nii(bFFE_filename1);        
    bffe1_0=mat2nifti(h.img,'nifti'); clear bFFE_filename1   
    h=load_nii(bFFE_filename2);      
    bffe2_0=mat2nifti(h.img,'nifti'); clear bFFE_filename2
    h=load_nii(bFFE_filename3);
    bffe3_0=mat2nifti(h.img,'nifti'); clear bFFE_filename3    
    h=load_nii(bFFE_filename4);  
    bffe1_120=mat2nifti(h.img,'nifti'); clear bFFE_filename4
    h=load_nii(bFFE_filename5);         
    bffe2_120=mat2nifti(h.img,'nifti'); clear bFFE_filename5
    h=load_nii(bFFE_filename6);
    bffe3_120=mat2nifti(h.img,'nifti'); clear bFFE_filename6
    h=load_nii(bFFE_filename7);         
    bffe1_240=mat2nifti(h.img,'nifti'); clear bFFE_filename7
    h=load_nii(bFFE_filename8);         
    bffe2_240=mat2nifti(h.img,'nifti'); clear bFFE_filename8
    h=load_nii(bFFE_filename9);         
    bffe3_240=mat2nifti(h.img,'nifti'); clear bFFE_filename9    
    end
    end
        
im1=zeros(x,y,mpm_par.no_slices,mpm_par.N_offsets,'single');
im2=im1;
im1(:,:,:,1)=bffe1_0;   clear bffe1_0
im1(:,:,:,2)=bffe1_120; clear bffe1_120
im1(:,:,:,3)=bffe1_240; clear bffe1_240
im2(:,:,:,1)=bffe2_0;   clear bffe2_0
im2(:,:,:,2)=bffe2_120; clear bffe2_120
im2(:,:,:,3)=bffe2_240; clear bffe2_240
if mpm_par.N_bFFE_ang==3
im3(:,:,:,1)=bffe3_0;   clear bffe3_0
im3(:,:,:,2)=bffe3_120; clear bffe3_120
im3(:,:,:,3)=bffe3_240; clear bffe3_240
end
elseif mpm_par.N_offsets==4
    
[~,~,ext] = fileparts(bFFE_filename1);     
if strcmpi(ext,'.mha')                       
h=mha_read_volume(bFFE_filename1); clear bFFE_filename1
bffe1_0=h.pixelData;
h=mha_read_volume(bFFE_filename2); clear bFFE_filename2
bffe2_0=h.pixelData;
h=mha_read_volume(bFFE_filename3); clear bFFE_filename3
bffe1_180=h.pixelData;
h=mha_read_volume(bFFE_filename4); clear bFFE_filename4
bffe2_180=h.pixelData;
h=mha_read_volume(bFFE_filename5); clear bFFE_filename5
bffe1_90=h.pixelData;
h=mha_read_volume(bFFE_filename6); clear bFFE_filename6
bffe2_90=h.pixelData;
h=mha_read_volume(bFFE_filename7); clear bFFE_filename7
bffe1_270=h.pixelData;
h=mha_read_volume(bFFE_filename8); clear bFFE_filename8
bffe2_270=h.pixelData;
elseif strcmpi(ext,'.nii') || strcmpi(ext,'.gz')
h=load_nii(bFFE_filename1);        
bffe1_0=mat2nifti(h.img,'nifti'); clear bFFE_filename1   
h=load_nii(bFFE_filename2);        
bffe2_0=mat2nifti(h.img,'nifti'); clear bFFE_filename2   
h=load_nii(bFFE_filename3);        
bffe1_180=mat2nifti(h.img,'nifti'); clear bFFE_filename3   
h=load_nii(bFFE_filename4);        
bffe2_180=mat2nifti(h.img,'nifti'); clear bFFE_filename4   
h=load_nii(bFFE_filename5);        
bffe1_90=mat2nifti(h.img,'nifti'); clear bFFE_filename5   
h=load_nii(bFFE_filename6);        
bffe2_90=mat2nifti(h.img,'nifti'); clear bFFE_filename6   
h=load_nii(bFFE_filename7);        
bffe1_270=mat2nifti(h.img,'nifti'); clear bFFE_filename7   
h=load_nii(bFFE_filename8);        
bffe2_270=mat2nifti(h.img,'nifti'); clear bFFE_filename8   
end

im1=zeros(x,y,mpm_par.no_slices,mpm_par.N_offsets,'single');
im2=im1;
im1(:,:,:,1)=bffe1_0;   clear bffe1_0
im1(:,:,:,2)=bffe1_90; clear bffe1_90
im1(:,:,:,3)=bffe1_180;  clear bffe1_180
im1(:,:,:,4)=bffe1_270; clear bffe1_270

im2(:,:,:,1)=bffe2_0;    clear bffe2_0
im2(:,:,:,2)=bffe2_90;  clear bffe2_90
im2(:,:,:,3)=bffe2_180;   clear bffe2_180
im2(:,:,:,4)=bffe2_270;  clear bffe2_270
end
%% Calculate T2 maps via the DESPOT2 technique
if mpm_par.N_offsets==2   
Y1(:,:,:,1)=abs(im1(:,:,:,1))./sind(mpm_par.bffe_ang1.*c);
X1(:,:,:,1)=abs(im1(:,:,:,1))./tand(mpm_par.bffe_ang1.*c);
Y1(:,:,:,2)=abs(im2(:,:,:,1))./sind(mpm_par.bffe_ang2.*c);
X1(:,:,:,2)=abs(im2(:,:,:,1))./tand(mpm_par.bffe_ang2.*c);
if mpm_par.N_bFFE_ang>=3
Y1(:,:,:,3)=abs(im3(:,:,:,1))./sind(mpm_par.bffe_ang3.*c);
X1(:,:,:,3)=abs(im3(:,:,:,1))./tand(mpm_par.bffe_ang3.*c);
end
if mpm_par.N_bFFE_ang==4
Y1(:,:,:,4)=abs(im4(:,:,:,1))./sind(mpm_par.bffe_ang4.*c);
X1(:,:,:,4)=abs(im4(:,:,:,1))./tand(mpm_par.bffe_ang4.*c);
end
[~,~,slope1,~]=linearize3D(Y1,X1);
%ind=int1<0; int1(ind)=-int1(ind); %ind=slope1-E1>=0; slope1(ind)=1; int1(ind)=0;
clear Y1 X1

Y2(:,:,:,1)=abs(im1(:,:,:,2))./sind(mpm_par.bffe_ang1.*c);
X2(:,:,:,1)=abs(im1(:,:,:,2))./tand(mpm_par.bffe_ang1.*c);
Y2(:,:,:,2)=abs(im2(:,:,:,2))./sind(mpm_par.bffe_ang2.*c);
X2(:,:,:,2)=abs(im2(:,:,:,2))./tand(mpm_par.bffe_ang2.*c);
if mpm_par.N_bFFE_ang>=3
Y2(:,:,:,3)=abs(im3(:,:,:,2))./sind(mpm_par.bffe_ang3.*c);
X2(:,:,:,3)=abs(im3(:,:,:,2))./tand(mpm_par.bffe_ang3.*c);
end
if mpm_par.N_bFFE_ang==4
Y2(:,:,:,4)=abs(im4(:,:,:,2))./sind(mpm_par.bffe_ang4.*c);
X2(:,:,:,4)=abs(im4(:,:,:,2))./tand(mpm_par.bffe_ang4.*c);
end
[~,~,slope2,~]=linearize3D(Y2,X2);
%ind=int2<0; int2(ind)=-int2(ind); %ind=slope2-E1>=0; slope2(ind)=1; int2(ind)=0;
clear X2 Y2

E1=exp(-mpm_par.bffe_TR./T1);
%ind1=slope1-E1>=0; ind2=slope2-E1>=0;

unity=ones(x,y,no_slices,'single');
E2_1=(slope1-E1)./(slope1.*E1-unity);
ind=E2_1>1; E2_1(ind)=-E2_1(ind); ind=E2_1<0; E2_1(ind)=epsilon_cutoff; ind=E2_1>1;
E2_1(ind)=1;
T2_1=-mpm_par.bffe_TR./real(log(E2_1));

E2_2=(slope2-E1)./(slope2.*E1-unity);
ind=E2_2>1;E2_2(ind)=-E2_2(ind); ind=E2_2<0; E2_2(ind)=epsilon_cutoff; ind=E2_2>1;
E2_2(ind)=1;
T2_2=-mpm_par.bffe_TR./real(log(E2_2));

T2_1(mask)=0; T2_2(mask)=0;

ind=T2_1>T2_clip;
T2_1(ind)=T2_clip;
ind=isnan(T2_1);
T2_1(ind)=0;

ind=T2_2>T2_clip;
T2_2(ind)=T2_clip;
ind=isnan(T2_2);
T2_2(ind)=0;

E2_anal=real(sqrt((2.*E2_2.*E2_1-E2_1-E2_2)./(E2_1+E2_2-2*unity)));
arg=(E2_1-E2_anal.^2)./(E2_anal.*(1-E2_1));
ind=isinf(arg);
arg(ind)=0;

T2_anal=-(mpm_par.bffe_TR-0.498*T_RFE(1,2))./log(E2_anal);
T2_anal=T2_anal.*logical(1-mask);

ind=T2_anal>T2_clip;
T2_anal(ind)=T2_clip;
ind=isnan(T2_anal);
T2_anal(ind)=0;
T2=T2_anal;
clear slope1 slope2 E1 E2_1 E2_2 E2_anal
elseif mpm_par.N_offsets==3
E1=exp(-mpm_par.bffe_TR./T1);
unity=ones(x,y,mpm_par.no_slices,'single');

Y1(:,:,:,1)=abs(im1(:,:,:,1))./sind(mpm_par.bffe_ang1.*c);
X1(:,:,:,1)=abs(im1(:,:,:,1))./tand(mpm_par.bffe_ang1.*c);
Y1(:,:,:,2)=abs(im2(:,:,:,1))./sind(mpm_par.bffe_ang2.*c);
X1(:,:,:,2)=abs(im2(:,:,:,1))./tand(mpm_par.bffe_ang2.*c);
if mpm_par.N_bFFE_ang==3
Y1(:,:,:,3)=abs(im3(:,:,:,1))./sind(mpm_par.bffe_ang3.*c);
X1(:,:,:,3)=abs(im3(:,:,:,1))./tand(mpm_par.bffe_ang3.*c);
end
[~,~,slope1,~]=linearize3D(Y1,X1);
%ind=int1<0; int1(ind)=-int1(ind);  
clear Y1 X1

Y2(:,:,:,1)=abs(im1(:,:,:,2))./sind(mpm_par.bffe_ang1.*c);
X2(:,:,:,1)=abs(im1(:,:,:,2))./tand(mpm_par.bffe_ang1.*c);
Y2(:,:,:,2)=abs(im2(:,:,:,2))./sind(mpm_par.bffe_ang2.*c);
X2(:,:,:,2)=abs(im2(:,:,:,2))./tand(mpm_par.bffe_ang2.*c);
if mpm_par.N_bFFE_ang==3
Y2(:,:,:,3)=abs(im3(:,:,:,2))./sind(mpm_par.bffe_ang3.*c);
X2(:,:,:,3)=abs(im3(:,:,:,2))./tand(mpm_par.bffe_ang3.*c);
end
[~,~,slope2,~]=linearize3D(Y2,X2);
%ind=int2<0; int2(ind)=-int2(ind);

clear Y2 X2

Y3(:,:,:,1)=abs(im1(:,:,:,3))./sind(mpm_par.bffe_ang1.*c);
X3(:,:,:,1)=abs(im1(:,:,:,3))./tand(mpm_par.bffe_ang1.*c);
Y3(:,:,:,2)=abs(im2(:,:,:,3))./sind(mpm_par.bffe_ang2.*c);
X3(:,:,:,2)=abs(im2(:,:,:,3))./tand(mpm_par.bffe_ang2.*c);
if mpm_par.N_bFFE_ang==3
Y3(:,:,:,3)=abs(im3(:,:,:,3))./sind(mpm_par.bffe_ang3.*c);
X3(:,:,:,3)=abs(im3(:,:,:,3))./tand(mpm_par.bffe_ang3.*c);
end
[~,~,slope3,~]=linearize3D(Y3,X3);
%ind=int3<0; int3(ind)=-int3(ind); 
clear Y3 X3

T2_1=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope1-E1)./(slope1.*E1-unity)));
T2_2=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope2-E1)./(slope2.*E1-unity)));
T2_3=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope3-E1)./(slope3.*E1-unity)));

ind=T2_1>T2_clip;
T2_1(ind)=T2_clip;
ind=isnan(T2_1);
T2_1(ind)=0;

ind=T2_2>T2_clip;
T2_2(ind)=T2_clip;
ind=isnan(T2_2);
T2_2(ind)=0;

ind=T2_3>T2_clip;
T2_3(ind)=T2_clip;
ind=isnan(T2_3);
T2_3(ind)=0;

ind1=(slope1-E1)>=0; ind2=(slope2-E1)>=0; ind3=(slope3-E1)>=0;
T2_1(ind1)=0; T2_2(ind2)=0; T2_3(ind3)=0;

T2_1(mask)=0; T2_2(mask)=0; T2_3(mask)=0;

T2=sqrt(8/9)*sqrt(T2_1.^2+T2_2.^2+T2_3.^2);
ind=T2>T2_clip;
T2(ind)=T2_clip;
           
elseif mpm_par.N_offsets==4
E1=exp(-mpm_par.bffe_TR./T1);    
Y1(:,:,:,1)=abs(im1(:,:,:,1))./sind(mpm_par.bffe_ang1.*c);
X1(:,:,:,1)=abs(im1(:,:,:,1))./tand(mpm_par.bffe_ang1.*c);
Y1(:,:,:,2)=abs(im2(:,:,:,1))./sind(mpm_par.bffe_ang2.*c);
X1(:,:,:,2)=abs(im2(:,:,:,1))./tand(mpm_par.bffe_ang2.*c);
[~,~,slope1,~]=linearize3D(Y1,X1);
%ind=int1<0; int1(ind)=-int1(ind); indA=slope1-E1>=0; %slope1(ind)=1; int1(ind)=0;
clear Y1 X1 ind

Y2(:,:,:,1)=abs(im1(:,:,:,2))./sind(mpm_par.bffe_ang1.*c);
X2(:,:,:,1)=abs(im1(:,:,:,2))./tand(mpm_par.bffe_ang1.*c);
Y2(:,:,:,2)=abs(im2(:,:,:,2))./sind(mpm_par.bffe_ang2.*c);
X2(:,:,:,2)=abs(im2(:,:,:,2))./tand(mpm_par.bffe_ang2.*c);
[~,~,slope2,~]=linearize3D(Y2,X2);
%ind=int2<0; int2(ind)=-int2(ind); indB=slope2-E1>=0; %slope2(ind)=1; int2(ind)=0;
clear Y2 X2 ind

Y3(:,:,:,1)=abs(im1(:,:,:,3))./sind(mpm_par.bffe_ang1.*c);
X3(:,:,:,1)=abs(im1(:,:,:,3))./tand(mpm_par.bffe_ang1.*c);
Y3(:,:,:,2)=abs(im2(:,:,:,3))./sind(mpm_par.bffe_ang2.*c);
X3(:,:,:,2)=abs(im2(:,:,:,3))./tand(mpm_par.bffe_ang2.*c);
[~,~,slope3,~]=linearize3D(Y3,X3);
%ind=int3<0; int3(ind)=-int3(ind); indC=slope3-E1>=0; %slope3(ind)=1; int3(ind)=0;
clear Y3 X3 ind

Y4(:,:,:,1)=abs(im1(:,:,:,4))./sind(mpm_par.bffe_ang1.*c);
X4(:,:,:,1)=abs(im1(:,:,:,4))./tand(mpm_par.bffe_ang1.*c);
Y4(:,:,:,2)=abs(im2(:,:,:,4))./sind(mpm_par.bffe_ang2.*c);
X4(:,:,:,2)=abs(im2(:,:,:,4))./tand(mpm_par.bffe_ang2.*c);
[~,~,slope4,~]=linearize3D(Y4,X4);
%ind=int4<0; int4(ind)=-int4(ind); indD=slope4-E1>=0; %slope4(ind)=1; int4(ind)=0;
clear Y4 X4 ind

unity=ones(x,y,mpm_par.no_slices,'single');

E2_1=(slope1-E1)./(slope1.*E1-unity);
ind=E2_1>1; E2_1(ind)=epsilon_cutoff; ind=E2_1<epsilon_cutoff; E2_1(ind)=epsilon_cutoff; ind=isnan(E2_1); E2_1(ind)=epsilon_cutoff;
%T2_1=-(bffe_TR-0.498*T_RFE)./real(log(E2_1));

E2_2=(slope3-E1)./(slope3.*E1-unity);
ind=E2_2>1; E2_2(ind)=epsilon_cutoff; ind=E2_2<epsilon_cutoff; E2_2(ind)=epsilon_cutoff; ind=isnan(E2_2); E2_2(ind)=epsilon_cutoff;
%T2_2=-(bffe_TR-0.498*T_RFE)./real(log(E2_2));

E2_3=(slope2-E1)./(slope2.*E1-unity);
ind=E2_3>1; E2_3(ind)=epsilon_cutoff; ind=E2_3<epsilon_cutoff; E2_3(ind)=epsilon_cutoff; ind=isnan(E2_3); E2_3(ind)=epsilon_cutoff;
%T2_3=-(bffe_TR-0.498*T_RFE)./real(log(E2_3));

E2_4=(slope4-E1)./(slope4.*E1-unity);
ind=E2_4>1; E2_4(ind)=epsilon_cutoff; ind=E2_4<epsilon_cutoff; E2_4(ind)=epsilon_cutoff; ind=isnan(E2_4); E2_4(ind)=epsilon_cutoff;
%T2_4=-(bffe_TR-0.498*T_RFE)./real(log(E2_4));

E2_anal1=real(sqrt((2*E2_1.*E2_2-E2_1-E2_2)./(E2_1+E2_2-2*unity)));
E2_anal2=real(sqrt((2*E2_3.*E2_4-E2_3-E2_4)./(E2_3+E2_4-2*unity)));

ind1=isnan(E2_anal1); ind2=isnan(E2_anal2); E2_anal1(ind1)=0; E2_anal2(ind2)=0;
ind1=E2_anal1<=0.1; ind2=E2_anal2<=0.1;

E2_anal1(ind1)=E2_anal2(ind1); E2_anal2(ind2)=E2_anal1(ind2);

arg_1=(E2_1-E2_anal1.^2)./(E2_anal1.*(1-E2_1)); ind=isnan(arg_1); arg_1(ind)=0;
arg_2=(E2_3-E2_anal2.^2)./(E2_anal2.*(1-E2_3)); ind=isnan(arg_2); arg_2(ind)=0;

T2_anal1=-(mpm_par.bffe_TR-0.498*T_RFE(1,2))./real(log(E2_anal1));
ind1=isnan(T2_anal1); T2_anal1(ind1)=0; clear ind1

T2_anal2=-(mpm_par.bffe_TR-0.498*T_RFE(1,2))./real(log(E2_anal2));
ind2=isnan(T2_anal2); T2_anal2(ind2)=0; clear ind2

T2_anal=((arg_2.^2).*T2_anal1+(arg_1.^2).*T2_anal2)./(arg_1.^2+arg_2.^2);

T2_1=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope1-E1)./(slope1.*E1-unity)));
T2_2=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope2-E1)./(slope2.*E1-unity)));
T2_3=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope3-E1)./(slope3.*E1-unity)));
T2_4=-(mpm_par.bffe_TR-0.555*T_RFE(1,2))./real(log((slope4-E1)./(slope4.*E1-unity)));

ind1=(slope1-E1)>=0; ind2=(slope2-E1)>=0; ind3=(slope3-E1)>=0; ind4=(slope4-E1)>=0;
T2_1(ind1)=0; T2_2(ind2)=0; T2_3(ind3)=0; T2_4(ind4)=0;
T2_1(mask)=0; T2_2(mask)=0; T2_3(mask)=0; T2_4(mask)=0;
clear ind1 ind2 ind3 ind4 

T2_anal1(mask)=0; T2_anal2(mask)=0;
arg_1(mask)=0; arg_2(mask)=0;
T2_anal(mask)=0;

T2=sqrt(2/3)*sqrt(T2_1.^2+T2_2.^2+T2_3.^2+T2_4.^2);
ind=T2>T2_clip;
T2(ind)=T2_clip; clear ind
clear slope1 slope2 slope3 slope4 E1 E2_1 E2_2 E2_3 E2_4
end
end
%% Calibrate M0 map
if (~strcmpi(calibration,'T1') && ~strcmpi(calculation,'M0&T1') ...
        && ~strcmpi(B1_corr,'no') && exist('M0','var')==1)
if exist('last_slice','var')==0 && strcmpi(mask1_im,'mask_im.mha')
disp('Select a cutoff location below the brain.');
figure; imshow(T1(:,:,87),[0 3000]);
[~,yi]=ginput(1); clf;
last_slice=round(yi);
else
last_slice=size(T1,2);
end
mask_crop=logical(1-mask(1:last_slice,:,:));
Npixels=sum(sum(sum(mask_crop)));
ind=isnan(M0); M0(ind)=0;

M0_vec=(0:0.03:150);
M0_crop=M0(1:last_slice,:,:);
M0_hist=hist(single(M0_crop(:)),M0_vec);
M0_hist_norm=M0_hist/Npixels;
M0_hist_norm(1,1:10)=0; M0_hist_norm(1,4900:end)=0;
M0_max=max(M0_hist_norm);
if strcmpi(cal_mode,'manual')
    disp('For a human brain, click on the WM and the GM peaks; otherwise');
    disp('for a phantom, click to the left and right of the PD peak.');
    figure; plot(M0_vec,M0_hist_norm); ylim([0 M0_max]); xlim([40 120]); 
    title('Normalized M0 histogram');
    [xi,~]=ginput(2); clf;
    WM_GM_centre=(xi(2)+xi(1))/2;
else
    if strcmpi(mask1_im,'mask_im_strip.nii.gz')
    options2=fitoptions('gauss3');
    options2.StartPoint=[1e-3 70 3 1e-3 81 6 0.05e-3 95 20];
    options2.Lower=[0.5e-3 62 2 0.5e-3 72 3 0.01e-3 80 5];
    options2.Upper=[1.5e-3 77 5 1.5e-3 87 8 0.2e-3 120 40];
    option2.Normalize='on'; option2.Robust='on';
    f2=fit(M0_vec',M0_hist_norm','gauss3',options2);
    WM_GM_centre=(f2.b1+f2.b2)/2;
    else
    options2=fitoptions('gauss4');
    options2.StartPoint=[1e-3 70 3 1e-3 81 6 0.5e-3 55 3 0.1e-3 100 20 ];
    options2.Lower=[0.3e-3 62 2 0.3e-3 72 3 0 40 6 0.01e-3 80 5 ];
    options2.Upper=[1.5e-3 77 5 1.5e-3 87 8 1.0e-3 64 12 0.2e-3 120 50 ];
    option2.Normalize='on'; option2.Robust='on';
    f2=fit(M0_vec',M0_hist_norm','gauss4',options2); disp(f2);
    WM_GM_centre=(f2.b1+f2.b2)/2;
    end    
end

    if strcmpi(mpm_par.Anatomy,'Brain')
    M0_corr=76.0/WM_GM_centre;
    else
    M0_corr=100/WM_GM_centre; %For a phantom normalize the PD to 100%
    end
    
M0norm=M0.*M0_corr;
M0=M0norm; clear M0_crop M0_hist mask_crop M0norm
end
    
if strcmpi(New_T1M0_mask,'yes') && (strcmpi(B1_corr,'BiasField'))
    index1=T1>700;
    index2=T1<2000;
    T1_mask2=logical(index1.*index2);
    index1=M0<90;
    index2=M0>50;
    M0_mask2=logical(index1.*index2);
    T1M0_mask=logical(T1_mask2.*M0_mask2);
    
xp=mpm_par.fov(1,1)/x; yp=mpm_par.fov(1,2)/y; zp=mpm_par.fov(1,3)/mpm_par.no_slices;    
metaData=struct('CompressedData','False','ObjectType','Image','NumOfDimensions',3,...
    'BinaryData','True','BinaryDataByteOrderMSB','False','ImageAxesOrientation',[0 0 -1 0 1 0 -1 0 0],...
    'Origin',[mpm_par.fov(1,1)/2 -mpm_par.fov(1,2)/2 mpm_par.fov(1,3)/2],'CenterOfRotation',...
    [mpm_par.fov(1,1)/2 -mpm_par.fov(1,2)/2 mpm_par.fov(1,3)/2],'AnatomicalOrientation','LPS',...
    'Spacing',[xp yp zp],'Size',[x y mpm_par.no_slices],'ElementType','MET_FLOAT','ElementDataFile','Local');
pixelData=single(T1M0_mask);
pixelWithMetaStruct=struct('metaData',metaData,'pixelData',pixelData);
mha_write_volume('T1M0_mask.mha',pixelWithMetaStruct,true);
clear index1 index2 T1mask2 M0mask2 xi yi x y last_slice
end
toc
%% Display sample parametric maps (axial, coronal, sagittal).
sag=100; ax=90; cor=128; 
           if ~strcmpi(cal_mode,'manual') && strcmpi(calibration,'M0&T1') && strcmpi(B1_corr,'BiasField')
y_max=max(T1_max, M0_max);
figure(8); subplot(121); plot(f1,T1_vec,T1_hist_norm); ylim([0 y_max]); xlim([0 4000]); title('Normalized T1 histogram before calibration');
           subplot(122); plot(f2,M0_vec,M0_hist_norm); ylim([0 y_max]); xlim([40 110]); title('Normalized M0 histogram before calibration');
           elseif ~strcmpi(cal_mode,'manual') && strcmpi(calibration,'M0')
figure(8); plot(f2,M0_vec,M0_hist_norm); ylim([0 y_max]); xlim([40 120]); title('Normalized M0 histogram before calibration');
           elseif ~strcmpi(cal_mode,'manual') && strcmpi(calibration,'T1')
figure(8); plot(f1,T1_vec,T1_hist_norm); ylim([0 y_max]); xlim([0 4000]); title('Normalized T1 histogram before calibration');
           end
           
figure(7); subplot(131); imshow(fliplr(squeeze(T1(ax,:,:))),[0 4000]); title('T_{1} axial');
           subplot(132); imshow(fliplr(squeeze(T1(:,cor,:))),[0 4000]); title('T_{1} coronal');
           subplot(133); imshow(squeeze(T1(:,:,sag)),[0 4000]); title('T_{1} sagittal');

           if (exist('T2s','var')==1)
figure(6); subplot(131); imshow(fliplr(squeeze(T2s(ax,:,:))),[0 250]); title('T_{2}^{*} axial');
           subplot(132); imshow(fliplr(squeeze(T2s(:,cor,:))),[0 250]); title('T_{2}^{*} coronal');
           subplot(133); imshow(squeeze(T2s(:,:,sag)),[0 250]); title('T_{2}^{*} sagittal');
           end
           if (exist('M0','var')==1)
figure(5); subplot(131); imshow(fliplr(squeeze(M0(ax,:,:))),[40 120]); title('PD axial');
           subplot(132); imshow(fliplr(squeeze(M0(:,cor,:))),[40 120]); title('PD coronal');
           subplot(133); imshow(squeeze(M0(:,:,sag)),[40 120]); title('PD sagittal');
           end
           if (exist('MTsat','var')==1)
figure(4); subplot(131); imshow(fliplr(squeeze(MTsat(ax,:,:))),[0 0.05]); title('MT_{sat} axial');
           subplot(132); imshow(fliplr(squeeze(MTsat(:,cor,:))),[0 0.05]); title('MT_{sat} coronal');
           subplot(133); imshow(squeeze(MTsat(:,:,sag)),[0 0.05]); title('MT_{sat} sagittal');
           end
           if (exist('T2','var')==1)
figure(3); subplot(131); imshow(fliplr(squeeze(T2(ax,:,:))),[0 350]); title('T_{2} axial');
           subplot(132); imshow(fliplr(squeeze(T2(:,cor,:))),[0 350]); title('T_{2} coronal');
           subplot(133); imshow(squeeze(T2(:,:,sag)),[0 350]); title('T_{2} sagittal');
           end
figure(2); subplot(131); imshow(fliplr(squeeze(c(ax,:,:))),[0.5 1.5]); title('c_{RF}^{+} axial'); colormap(gca, jet(256)); 
           subplot(132); imshow(fliplr(squeeze(c(:,cor,:))),[0.5 1.5]); title('c_{RF}^{+} coronal'); colormap(gca, jet(256));
           subplot(133); imshow(squeeze(c(:,:,sag)),[0.5 1.5]); title('c_{RF}^{+} sagittal'); colormap(gca, jet(256));
           if (exist('c_minus','var')==1)
figure(1); subplot(131); imshow(fliplr(squeeze(c_minus(ax,:,:))),[0.5 1.5]); title('c_{RF}^{-} axial'); colormap(gca, jet(256)); 
           subplot(132); imshow(fliplr(squeeze(c_minus(:,cor,:))),[0.5 1.5]); title('c_{RF}^{-} coronal'); colormap(gca, jet(256)); 
           subplot(133); imshow(squeeze(c_minus(:,:,sag)),[0.5 1.5]); title('c_{RF}^{-} sagittal'); colormap(gca, jet(256)); 
           end
clearvars -except M0 T1 T1_uncorr T1_uncal T2 T2_anal T2_anal1 T2_anal2 T2s MTsat c c_minus mask cal_new             

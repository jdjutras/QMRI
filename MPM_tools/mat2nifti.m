function [output]=mat2nifti(input,option)

if strcmpi(option,'mat')
    s=size(input);
    if length(s)==3
        output=zeros(s(3),s(1),s(2),'single');
        for i=1:s(3)
            output(i,:,:)=imrotate(fliplr(squeeze(input(:,:,i))),-90);
        end
    elseif length(s)==4
        output=zeros(s(3),s(1),s(2),s(4),'single');
        for j=1:s(4)
        for i=1:s(3)
            output(i,:,:,j)=imrotate(fliplr(squeeze(input(:,:,i,j))),-90);
        end
        end
    end
elseif strcmpi(option,'nifti')
        s=size(input);
    if length(s)==3
        output=zeros(s(2),s(3),s(1),'single');
        for i=1:s(1)
            output(:,:,i)=fliplr(imrotate(squeeze(input(i,:,:)),90));
        end
    elseif length(s)==4
        output=zeros(s(2),s(3),s(1),s(4),'single');
        for j=1:s(4)
        for i=1:s(1)
            output(:,:,i,j)=fliplr(imrotate(squeeze(input(i,:,:,j)),90));
        end
        end
    end
end

    
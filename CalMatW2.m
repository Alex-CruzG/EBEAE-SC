function [W,D]=CalMatW2(X,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [W,D]=CalMatW(X,para)
%
% Function to calculate the weight matrix W and vector D for the PISINMF method.
%
% Input Arguments
%   X    = matrix of measurements
%   para = structure type variable
%       \mu =   Variable de regularizacion
%       Tven =  size of window Tven*Tven
%       dimX =  X dimension of the image to be analyzed
%       dimY =  Y dimension of the image to be analyzed
%
% Output Arguments
%   W = weight matrix
%   D = diagonal matrix of W
%
% Ines A. Cruz-Guerrero
% May/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(para.tven,2)==0
    disp('The entered number of window size is even, the default value will be taken.')
    sq=3;
else
    sq=para.tven;
end

X=reshape(X',para.dimX,para.dimY,size(X,1));

sp=floor(sq/2);
W=sparse(para.dimX*para.dimY,para.dimX*para.dimY);
phi=zeros(para.dimX,para.dimY);
H = ones(para.tven);
kH=zeros(para.tven^2,2);
for i=1:para.tven
    kH((i-1)*para.tven+1:i*para.tven,1)=i;
    kH((i-1)*para.tven+1:i*para.tven,2)=1:para.tven;
end
cen=[(sp+1) (sp+1)];
Hdis=reshape(pdist2(cen,kH,'euclidean'),para.tven,para.tven);

for i=1:para.dimX
    for j=1:para.dimY
        C_k=zeros(para.dimX,para.dimY);
        
        if (i-sp>0 && i+sp<para.dimX)&&(j-sp>0 && j+sp<para.dimY)
            C_k(i-sp:i+sp,j-sp:j+sp)=H;
            Uc=Hdis;
            kcom=X(i-sp:i+sp,j-sp:j+sp,:);
        else    
            if i-sp<=0 && j-sp<=0
                C_k(1:i+sp,1:j+sp)=H(sq-(sp+(i-1)):sq,sq-(sp+(j-1)):sq);
                Uc=Hdis(sq-(sp+(i-1)):sq,sq-(sp+(j-1)):sq);
                kcom=X(1:i+sp,1:j+sp,:);
            elseif i+sp>=para.dimX && j+sp>=para.dimY
                C_k(i-sp:para.dimX,j-sp:para.dimY)=H(1:sq-((i+sp)-para.dimX),1:sq-((j+sp)-para.dimY));
                Uc=Hdis(1:sq-((i+sp)-para.dimX),1:sq-((j+sp)-para.dimY));
                kcom=X(i-sp:para.dimX,j-sp:para.dimY,:);
            elseif i-sp<=0 && j+sp>para.dimY
                C_k(1:i+sp,j-sp:para.dimY)=H(sq-(sp+(i-1)):sq,1:sq-((j+sp)-para.dimY));
                Uc=Hdis(sq-(sp+(i-1)):sq,1:sq-((j+sp)-para.dimY));
                kcom=X(1:i+sp,j-sp:para.dimY,:);
            elseif i+sp>para.dimX && j-sp<=0
                C_k(i-sp:para.dimX,1:j+sp)=H(1:sq-((i+sp)-para.dimX),sq-(sp+(j-1)):sq);
                Uc=Hdis(1:sq-((i+sp)-para.dimX),sq-(sp+(j-1)):sq);
                kcom=X(i-sp:para.dimX,1:j+sp,:);
            elseif i-sp<=0
                C_k(1:i+sp,j-sp:j+sp)=H(sq-(sp+(i-1)):sq,:);
                Uc=Hdis(sq-(sp+(i-1)):sq,:);
                kcom=X(1:i+sp,j-sp:j+sp,:);
            elseif j-sp<=0
                C_k(i-sp:i+sp,1:j+sp)=H(:,sq-(sp+(j-1)):sq);
                Uc=Hdis(:,sq-(sp+(j-1)):sq);
                kcom=X(i-sp:i+sp,1:j+sp,:);
            elseif i+sp>=para.dimX
                C_k(i-sp:para.dimX,j-sp:j+sp)=H(1:sq-((i+sp)-para.dimX),:);
                Uc=Hdis(1:sq-((i+sp)-para.dimX),:);
                kcom=X(i-sp:para.dimX,j-sp:j+sp,:);
            elseif j+sp>=para.dimY
                C_k(i-sp:i+sp,j-sp:para.dimY)=H(:,1:sq-((j+sp)-para.dimY));
                Uc=Hdis(:,1:sq-((j+sp)-para.dimY));
                kcom=X(i-sp:i+sp,j-sp:para.dimY,:);
            end
        end
        sigma=0;
        aux=0;
        % C_k(i,j)=0;
        v=abs(acos(sum(X(i,j,:).*kcom,3)./(sqrt(sum(X(i,j,:).^2,3)).*sqrt(sum(kcom.^2,3)))));
        M=1./(sqrt(Uc).*v);
        M(isinf(M))=0;
        aux=sum((X(i,j,:)-kcom).^2,3);
        sigma=sum(aux(:))./(sum(C_k(:))-1);
        Wn=M.*exp(-aux/sigma);
        Wn(isinf(Wn))=0;
        phi(i,j)=sum(Wn(:));
        C_k(logical(C_k))=Wn;
        W(:,(i-1)*para.dimX+j)=C_k(:);
%         C_k(:,:,(i-1)*sc(3)+j)=C_k(:,:,(i-1)*sc(3)+j)/length(C_k(C_k(:,:,(i-1)*sc(3)+j)>0));
    end
end
D=sparse(spdiags(phi(:)));
    
end
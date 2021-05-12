function [Y,Po,A,u,t]=VNIRsynth(N,Npixels,SNR,PSNR)

%%%%%%%%%%%%%%%%%%%%%%%
%
% [Y,P,A,u,t]=VNIRsynth(N,Npixels,SNR,PSNR)
%
% INPUTS
% N --> Order of multi-exponential model
% Npixels --> numbers of pixels in x & y axes
% SNR --> SNR of Gaussian noise (dB)
% PSNR --> PSNR for Shot noise (dB)
%
% OUTPUTS
% Y --> matrix of fluorescence decays of size 186 x (Npixels*Npixels)
% A --> matrix of abundances of N x (Npixels*Npixels)
% P --> matrix of end-members 186 x N
% u --> vector of laser input
% t --> vector of time samples
%
% July/2019
% DUCD
%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthetic FLIM Dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0,
    N=2;
    Npixels=60;
    Ts=0.25e-9;
    SNR=40;
    PSNR=20;
elseif nargin<5,
    PSNR=0;
elseif nargin<4,
    PSNR=0; SNR=0;
elseif nargin<3,
    PSNR=0; SNR=0;Ts=0.25e-9;
elseif nargin<2,
    PSNR=0; SNR=0;Ts=0.25e-9; Npixels=60;
end;

if N>4,
    N=4; disp('The maximum number of components is 4!!');
end;

if SNR ~= 0,
   NoiseGaussian=1;
else
   NoiseGaussian=0;
end;
if PSNR ~= 0,
   NoiseShot=1; 
else
    NoiseShot=0;
end;
if SNR ~= 0 | PSNR ~= 0,
    NoiseMeasurement=1;        
else
     NoiseMeasurement=0;
end;

Nsamp=Npixels;
x=1:Npixels;
y=1:Npixels;
[xx,yy]=meshgrid(x,y);
K=Npixels*Npixels;
load EndMembersVNIR;
L=size(P,1);
P1=P(:,1);
P2=P(:,2);
P3=P(:,3);
P4=P(:,4);


if N==2,

    aa1=7*exp(-0.001*(xx-Nsamp/2).^2-0.001*(yy-Nsamp/2).^2)+0.5;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.001*(yy).^2)+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2)+...
            2.5*exp(-0.001*(xx).^2-0.0001*(yy-Nsamp).^2)+2.5*exp(-0.0001*(xx-Nsamp).^2-0.001*(yy-Nsamp).^2);

     a1=zeros(Nsamp,Nsamp);
     a2=zeros(Nsamp,Nsamp);
    for i=1:Nsamp,
        for l=1:Nsamp,
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l));
        end;
    end;

elseif N==3,

    aa1=2.5*exp(-0.0025*(xx-Nsamp/2).^2-0.0025*(yy-Nsamp/2).^2)+0;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.001*(yy).^2)+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2);
    aa3=2.5*exp(-0.001*(xx).^2-0.0001*(yy-Nsamp).^2)+2.5*exp(-0.0001*(xx-Nsamp).^2-0.001*(yy-Nsamp).^2);

    a1=zeros(Nsamp,Nsamp);
    a2=zeros(Nsamp,Nsamp);
    a3=zeros(Nsamp,Nsamp);

    for i=1:Nsamp,
        for l=1:Nsamp,
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l));
            a3(i,l)=aa3(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l));
        end;
    end;
  
elseif N==4,
    
    aa1=2.5*exp(-0.005*(xx-Nsamp/2).^2-0.0005*(yy-Nsamp/2).^2)+0.5;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.0001*(yy).^2)+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2);
    aa3=2.5*exp(-0.001*(xx).^2-0.0001*(yy-Nsamp).^2);
    aa4=2.5*exp(-0.0001*(xx-Nsamp).^2-0.001*(yy-Nsamp).^2);
    
    a1=zeros(Nsamp,Nsamp);
    a2=zeros(Nsamp,Nsamp);
    a3=zeros(Nsamp,Nsamp);
    a4=zeros(Nsamp,Nsamp);
    for i=1:Nsamp,
        for l=1:Nsamp,
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a3(i,l)=aa3(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a4(i,l)=aa4(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
        end;
    end;
end;

Yy=zeros(Nsamp,Nsamp,L);
for i=1:Nsamp;
    for j=1:Nsamp;
        if N==2,     
            y=a1(i,j)*P1+a2(i,j)*P2;
        elseif N==3,
            y=a1(i,j)*P1+a2(i,j)*P2+a3(i,j)*P3;
        elseif N==4,
            y=a1(i,j)*P1+a2(i,j)*P2+a3(i,j)*P3+a4(i,j)*P4;
        end;
  
        if NoiseMeasurement==1 & NoiseGaussian==1,
             sigmay=sqrt((y'*y)/(10^(SNR/10)));
             yy1=sigmay*randn(L,1);
        else
            yy1=zeros(L,1);
        end;
        if NoiseMeasurement==1 & NoiseShot==1,
            sigmay=sqrt(max(y)/(10^(PSNR/10)));
            yy2=sigmay*randn(L,1).*sqrt(abs(y));
        else
            yy2=zeros(L,1);
        end;     
        Yy(i,j,:)=y+yy1+yy2;
   
    end;
end;

if N==2,     
    Po=[P1 P2];
    A=[reshape(a1,1,K);reshape(a2,1,K)];
elseif N==3,
    Po=[P1 P2 P3];
    A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K)];
elseif N==4,
    Po=[P1 P2 P3 P4];
    A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K); reshape(a4,1,K)];
end;

Y=reshape(Yy,K,L)';





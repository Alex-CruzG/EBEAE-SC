function [Y,P,A,u,t]=mflimsynth(N,Npixels,Ts,SNR,PSNR)

%%%%%%%%%%%%%%%%%%%%%%%
%
% [Y,P,A,u,t]=mflimsynth(N,Npixels,Ts,SNR,PSNR)
%
% INPUTS
% N --> Order of multi-exponential model
% Npixels --> numbers of pixels in x & y axes
% Ts --> sampling time
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
load dataexp1;
L=length(u);
u=u/sum(u);
n=0:(L-1);
t=n*Ts;
U=toeplitz(u',[u(1) zeros(1,L-1)]);

tau11=10;
tau12=20;
tau13=5;
tau14=15;

tau21=25;
tau22=10;
tau23=15;
tau24=5;

tau31=5;
tau32=7;
tau33=35;
tau34=15;

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

    aa1=7*exp(-0.005*(xx-Nsamp/2).^2-0.005*(yy-Nsamp/2).^2)+0.5;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.001*(yy).^2)+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2);
    aa3=3.5*exp(-0.001*(xx).^2-0.0001*(yy-Nsamp).^2)+2.5*exp(-0.0001*(xx-Nsamp).^2-0.001*(yy-Nsamp).^2);

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
    
    aa1=2.5*exp(-0.005*(xx-Nsamp/2).^2-0.0005*(yy-Nsamp/2).^2)+0;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.00025*(yy).^2);%+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2);
    aa3=2.5*exp(-0.001*(xx).^2-0.0002*(yy-Nsamp).^2);
    aa4=2.5*exp(-0.001*(xx-8*Nsamp/9).^2-0.001*(yy-8*Nsamp/9).^2)+2.5*exp(-0.001*(xx-Nsamp/9).^2-0.001*(yy-8*Nsamp/9).^2);
    
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

P11=U*0.6*exp(-n'/tau11);
P21=U*0.2*exp(-n'/tau21);
P31=U*0.2*exp(-n'/tau31);
P1=[P11;P21;P31];

if N>=2,
    P12=U*0.2*exp(-n'/tau12);
    P22=U*0.6*exp(-n'/tau22);
    P32=U*0.2*exp(-n'/tau32);
    P2=[P12;P22;P32];
    if N>=3,
        P13=U*0.15*exp(-n'/tau13);
        P23=U*0.15*exp(-n'/tau23);
        P33=U*0.70*exp(-n'/tau33);
        P3=[P13;P23;P33];
        if N>=4,
            P14=U*0.0*exp(-n'/tau14);
            P24=U*0.4*exp(-n'/tau24);
            P34=U*0.6*exp(-n'/tau34);
            P4=[P14;P24;P34];
        end;
    end;
end;
Yy=zeros(Nsamp,Nsamp,3*L);
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
             yy1=sigmay*randn(3*L,1);
        else
            yy1=zeros(3*L,1);
        end;
        if NoiseMeasurement==1 & NoiseShot==1,
            sigmay=sqrt(max(y)/(10^(PSNR/10)));
            yy2=sigmay*randn(3*L,1).*sqrt(abs(y));
        else
            yy2=zeros(3*L,1);
        end;     
        Yy(i,j,:)=y+yy1+yy2;
   
    end;
end;

if N==2,     
    P=[P1 P2];
    A=[reshape(a1,1,K);reshape(a2,1,K)];
elseif N==3,
    P=[P1 P2 P3];
    A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K)];
elseif N==4,
    P=[P1 P2 P3 P4];
    A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K); reshape(a4,1,K)];
end;

Y=reshape(Yy,K,3*L)';





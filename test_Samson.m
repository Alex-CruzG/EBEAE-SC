%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Evaluation of EBEAE-TV with Real Hyperspectral Dataset
%
% Samson    -->     VNIR HSI captured by the NASA AVIRIS system
%                   available  at  https://rslab.ut.ac.ir/data
%
% Ines A. Cruz-Guerrero
% May/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

addpath('Data_base\Samson');
addpath('EBEAE');
addpath('NMF-QMV');
addpath('GraphL');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3;                % Number of End-members
Nsamples=95;       % Size of the Squared Image Nsamples x Nsamples 
SNR=35;             % Level in dB of Gaussian Noise SNR=45,50,55
PSNR=10;            % Level in dB of Shot Noise PSNR=15,20,25

load samson_1.mat
load end3.mat

Po=normalize(M,'norm',1);
Ao=normalize(A,'norm',1);
Xim = reshape(Ao',nRow,nCol,N);
rows=nRow;
columns=nCol;
L=nBand;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Ground-Truths of Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(5,1,1);
plot(Po,'LineWidth',2); grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('A) Ground-truth');
legend('End-member 1','End-member 2','End-member 3');

figure(2)
for i=1:N
    eval(['subplot(5' num2str(N) num2str(i) ');']);
    eval(['imagesc(Xim(:,:,' num2str(i) '),[0,1]);']);
    title(['End-member ' num2str(i)]);
    colormap jet
end
hold on;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the observed data Y with noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=rows*columns;
Y=Po*Ao;
Y=AddNoiseFLIM(Y,SNR,PSNR);
Y(Y<0) = 0;                             % Forza valores negativos a cero
Y = Y./repmat(sum(Y,1),L,1);            % Condicion de suma a 1

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute GLNMF Methodology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Samson Hyperspecytal Dataset');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('GLNMF Analysis');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

Pini = VCA(Y,'Endmembers',N,'SNR',1,'verbose','no');
Aini = FCLSU(Y, Pini)';
para_nmf.tol = 1e-3;
para_nmf.itermax = 20;
para_nmf.lambda = 1;
para_nmf.mu = 1e-4;

tic;
[iter, P1, A1]= glnmf(Y, N, Pini, Aini, para_nmf);
T_m1 = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1=normalize(P1,'norm',1);
A1=normalize(A1,'norm',1);
[A1, P1, ~, ~] = find_perm(Ao,Po, A1, P1);

figure(1);
subplot(5,1,2);
plot(P1,'LineWidth',2); grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('B) GLNMF');
legend('End-member 1','End-member 2','End-member 3');

figure(2)
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+N) ');']);
    eval(['imagesc(reshape(A1(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    colormap jet
    if i==2, title('B) GLNMF '); end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute PISINMF Methodology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('PISINMF');

para.dimX= Nsamples;
para.dimY= Nsamples;
para.tven= 5 ;
para.tau=1e-3;
para.maxiter=20;
para.delta=10;
para.mu= 1e-04;
para.t=35;
para.alpha=1e-03;

tic
[P3,A3] =  PISINMF(Y,N,para);
T_m3=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P3=normalize(P3,'norm',1);
A3=normalize(A3,'norm',1);
[A3, P3, ~, ~] = find_perm(Ao,Po, A3, P3);

figure(1);
subplot(5,1,3);
plot(P3,'LineWidth',2); grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('C) PISINMF');
legend('End-member 1','End-member 2','End-member 3');

figure(2);
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+2*N) ');']);
    eval(['imagesc(reshape(A3(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    colormap jet
    if i==2
        title('C) PISINMF');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute NMF-QMV with Total variation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('NMF-QMV with Total variation');

term ='totalVar';
beta_candidates =   10.^(-5:5);

tic;
[beta_best, P4, A4, results_save] = NMF_QMV(Y, N, beta_candidates, term);
T_m4=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P4=normalize(P4,'norm',1);
A4=normalize(A4,'norm',1);
[A4, P4, ~, ~] = find_perm(Ao,Po, A4, P4);

figure(1);
subplot(5,1,4);
plot(P4,'LineWidth',2); grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('D) NMF-QMV');
legend('End-member 1','End-member 2','End-member 3');

figure(2);
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+3*N) ');']);
    eval(['imagesc(reshape(A4(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    colormap jet
    if i==2
        title('D) NMF-QMV');
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute EBEAE-TV Methodology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('EBEAE-TV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of EBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initcond=1;
rho=0.9;
lambda=0.0;
epsilon=1e-3;
maxiter=50;
parallel=0;
normalization=1;
disp_iter=0;

mu=1e-04; nu=1e-05; tau=1e-07;
sc=[mu,nu,tau, Nsamples,Nsamples];
paramvec=[initcond,rho,lambda,epsilon,maxiter,parallel,normalization,disp_iter];

tic;
[P5,A,A5,Yh6]=EBEAE_TV(Y,N,paramvec,sc);
T_m5=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A5, P5, ~, ~] = find_perm(Ao,Po, A5, P5);

figure(1);
subplot(5,1,5);
plot(P5,'LineWidth',2); grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('E) EBEAE-TV');
legend('End-member 1','End-member 2','End-member 3');

figure(2);
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+4*N) ');']);
    eval(['imagesc(reshape(A5(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    colormap jet
    if i==2
        title('E) EBEAE-TV');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Estimation Errors on Abundances and End-members, and Execution
% Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Average errors');
Er_p1=[];
Er_a1=[];
Er_p3=[];
Er_a3=[];
Er_p4=[];
Er_a4=[];
Er_p5=[];
Er_a5=[];
for i=1:N
        Er_p1=[Er_p1 norm(Po(:,i)-P1(:,i))];        
        Er_a1=[Er_a1 norm(Ao(i,:)-A1(i,:))];
        Er_p3=[Er_p3 norm(Po(:,i)-P3(:,i))];        
        Er_a3=[Er_a3 norm(Ao(i,:)-A3(i,:))];
        Er_p4=[Er_p4 norm(Po(:,i)-P4(:,i))];        
        Er_a4=[Er_a4 norm(Ao(i,:)-A4(i,:))];
        Er_p5=[Er_p5 norm(Po(:,i)-P5(:,i))];        
        Er_a5=[Er_a5 norm(Ao(i,:)-A5(i,:))];
end

disp('%%%%%%%%%%%%%%%%%%%');
disp('Performance Metrics');
disp(['GLNMF Execution time=' num2str(T_m1) ' s']);
disp(['PISINMF Execution time=' num2str(T_m3) ' s']);
disp(['NMF-QMV Execution time=' num2str(T_m4) ' s']);
disp(['EBEAE-TV Execution time=' num2str(T_m5) ' s']); 
disp('%%%%%%%%%%%%%%%%%%%');
disp(['GLNMF Estimation Error in Measurements=' num2str(norm(Y-P1*A1,'fro')/K)]);
disp(['PISINMF Estimation Error in Measurements=' num2str(norm(Y-P3*A3,'fro')/K)]);
disp(['NMF-QMV Estimation Error in Measurements=' num2str(norm(Y-P4*A4,'fro')/K)]);
disp(['EBEAE-TV Estimation Error in Measurements=' num2str(norm(Y-P5*A5,'fro')/K)]);
disp('%%%%%%%%%%%%%%%%%%%');
disp(['GLNMF Estimation Error in End-members=' num2str(mean(Er_p1))]);
disp(['PISINMF Estimation Error in End-members=' num2str(mean(Er_p3))]);
disp(['NMF-QMV Estimation Error in End-members=' num2str(mean(Er_p4))]);
disp(['EBEAE-TV Estimation Error in End-members=' num2str(mean(Er_p5))]);
disp('%%%%%%%%%%%%%%%%%%%');
disp(['GLNMF Estimation Error in Abundances=' num2str(mean(Er_a1))]);
disp(['PISINMF Estimation Error in Abundances=' num2str(mean(Er_a3))]);
disp(['NMF-QMV Estimation Error in Abundances=' num2str(mean(Er_a4))]);
disp(['EBEAE-TV Estimation Error in Abundances=' num2str(mean(Er_a5))]);
% 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Distance from a point to a set');
Er_p1=[];
Er_a1=[];
Er_p3=[];
Er_a3=[];
Er_p4=[];
Er_a4=[];
Er_p5=[];
Er_a5=[];
for i=1:N
    for j=1:N
        Er_p1=[Er_p1 norm(Po(:,i)-P1(:,j))];        
        Er_a1=[Er_a1 norm(Ao(i,:)-A1(j,:))];
        Er_p3=[Er_p3 norm(Po(:,i)-P3(:,j))];        
        Er_a3=[Er_a3 norm(Ao(i,:)-A3(j,:))];
        Er_p4=[Er_p4 norm(Po(:,i)-P4(:,j))];        
        Er_a4=[Er_a4 norm(Ao(i,:)-A4(j,:))];
        Er_p5=[Er_p5 norm(Po(:,i)-P5(:,j))];        
        Er_a5=[Er_a5 norm(Ao(i,:)-A5(j,:))];
    end
end
disp('%%%%%%%%%%%%%%%%%%%');
disp('Performance Metrics');
disp(['GLNMF Execution time=' num2str(T_m1) ' s']);
disp(['PISINMF Execution time=' num2str(T_m3) ' s']);
disp(['NMF-QMV Execution time=' num2str(T_m4) ' s']);
disp(['EBEAE-TV Execution time=' num2str(T_m5) ' s']); 
disp('%%%%%%%%%%%%%%%%%%%');
disp(['GLNMF Estimation Error in Measurements=' num2str(norm(Y-P1*A1,'fro')/K)]);
disp(['PISINMF Estimation Error in Measurements=' num2str(norm(Y-P3*A3,'fro')/K)]);
disp(['NMF-QMV Estimation Error in Measurements=' num2str(norm(Y-P4*A4,'fro')/K)]);
disp(['EBEAE-TV Estimation Error in Measurements=' num2str(norm(Y-P5*A5,'fro')/K)]);
disp('%%%%%%%%%%%%%%%%%%%');
disp(['GLNMF Estimation Error in End-members=' num2str(min(Er_p1)/(2*N))]);
disp(['PISINMF Estimation Error in End-members=' num2str(min(Er_p3)/(2*N))]);
disp(['NMF-QMV Estimation Error in End-members=' num2str(min(Er_p4)/(2*N))]);
disp(['EBEAE-TV Estimation Error in End-members=' num2str(min(Er_p5)/(2*N))]);
disp('%%%%%%%%%%%%%%%%%%%');
disp(['GLNMF Estimation Error in Abundances=' num2str(min(Er_a1)/(2*N))]);
disp(['PISINMF Estimation Error in Abundances=' num2str(min(Er_a3)/(2*N))]);
disp(['NMF-QMV Estimation Error in Abundances=' num2str(min(Er_a4)/(2*N))]);
disp(['EBEAE-TV Estimation Error in Abundances=' num2str(min(Er_a5)/(2*N))]);

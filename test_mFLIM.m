%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Synthetic Evaluation of EBEAE-TV with m-FLIM Hyperspecytal Dataset
%
% m-FLIM -->  microscopic  Fluorescence  Lifetime  Imagin
%
% Ines A. Cruz-Guerrero
% May/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

addpath('EBEAE');
addpath('NMF-QMV');
addpath('GraphL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Synthetic Dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=4;                % Number of Simulated End-members only n=2,3,4
Nsamples=80;       % Size of the Squared Image Nsamples x Nsamples 
SNR=35;             % Level in dB of Gaussian Noise SNR=45,50,55
PSNR=10;            % Level in dB of Shot Noise PSNR=15,20,25
[Y,Po,Ao]=mflimsynth(N,Nsamples,250e-12,SNR,PSNR);      % Synthetic mFLIM
[L,K]=size(Y);
MHSI=reshape(Y',Nsamples,Nsamples,L);
bands=linspace(450,950,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Ground-Truths of Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(5,1,1);
plot(normalize(Po,'norm',1),'LineWidth',2); grid on;
axis([0 560 0 max(max(normalize(Po,'norm',1)))]);
xlabel('Time Samples');
ylabel('Normalized Intensity');
title('A) Ground-truth');
legend('Endmember 1','Endmember 2','Endmember 3','Endmember 4'); hold on;

figure(2);
for i=1:N
    eval(['subplot(5' num2str(N) num2str(i) ');']);
    eval(['imagesc(reshape(Ao(' num2str(i) ',:),Nsamples,Nsamples),[0,1]);']);
    title(['Endmember ' num2str(i)]);
end
hold on;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute GLNMF Methodology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Synthetic m-flim Hyperspecytal Dataset');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('GLNMF Analysis');

Pini = VCA(Y,'Endmembers',N,'SNR',1,'verbose','no');
Aini = FCLSU(Y, Pini)';
para_nmf.tol = 1e-3;
para_nmf.itermax = 20;
para_nmf.lambda = 0.04;
para_nmf.mu = 0.1;

tic;
[iter, P1, A1]= glnmf(Y, N, Pini, Aini, para_nmf);
T_m1 = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A1, P1, ~, ~] = find_perm(Ao,Po, A1, P1);

figure(2)
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+N) ');']);
    eval(['imagesc(reshape(A1(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    if i==2, title('B) GLNMF Estimation'); end
end

figure(1);
subplot(5,1,2);
plot(P1,'LineWidth',2); grid on;
axis([0 560 0 max(P1(:))]);
xlabel('Time Samples');
ylabel('Normalized Intensity');
title('B) GLNMF Estimation');
legend('Endmember 1','Endmember 2','Endmember 3','Endmember 4'); hold on;

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
para.delta=40;
para.mu= 0.6;
para.t=45;
para.alpha=0.1;

tic
[P3,A3] =  PISINMF(Y,N,para);
T_m3=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated Abundances and End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A3, P3, ~, ~] = find_perm(Ao,Po, A3, P3);

figure(1);
subplot(5,1,3);
plot(P3,'LineWidth',2); grid on;
axis([0 560 0 max(P3(:))]);
xlabel('Time Samples');
ylabel('Normalized Intensity');
title('C) PISINMF Estimation');
legend('Endmember 1','Endmember 2','Endmember 3','Endmember 4'); hold on;

figure(2);
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+2*N) ');']);
    eval(['imagesc(reshape(A3(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    if i==2
        title('C) PISINMF Estimation');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute NMF-QMV Methodology
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

[A4, P4, ~, ~] = find_perm(Ao,Po, A4, P4);

figure(1);
subplot(5,1,4);
plot(P4,'LineWidth',2); grid on;
axis([0 560 0 max(P4(:))]);
xlabel('Time Samples');
ylabel('Normalized Intensity');
title('D) NMF-QMV Estimation');
legend('Endmember 1','Endmember 2','Endmember 3','Endmember 4'); hold on;

figure(2);
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+3*N) ');']);
    eval(['imagesc(reshape(A4(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    if i==2
        title('D) NMF-QMV Estimation');
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
initcond=3;
rho=1;
lambda=0;
epsilon=1e-3;
maxiter=20;
parallel=0;
normalization=1;
disp_iter=0;

mu=0.01; nu=0.1; tau=0.1;
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
axis([0 560 0 max(max(normalize(Po,'norm',1)))]);
xlabel('Time Samples');
ylabel('Normalized Intensity');
title('F) EBEAE-STV Estimation');
legend('Endmember 1','Endmember 2','Endmember 3','Endmember 4'); hold on;

figure(2);
for i=1:N
    eval(['subplot(5,' num2str(N) ',' num2str(i+4*N) ');']);
    eval(['imagesc(reshape(A5(' num2str(i) ',:),Nsamples,Nsamples),[0 1]);']);
    if i==2
        title('E) EBEAE-STV Estimation');
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
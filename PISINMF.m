function [G,A] =  PISINMF(X,m,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [G,A] =  PISINMF(X,m,para)
%
% Spectral Unmixing of Hyperspectral Remote Sensing Imagery via Preserving 
%   the Intrinsic Structure Invariant
%
% Input Arguments
%   X   = matrix of measurements (MxN)
%   n   = order of linear mixture model
%   para= vector of hyper-parameters
%       dimX    = Vertical dimention 
%       dimY    = Horizontal dimention
%       tven    = Size of window
%       tau     = Threshold for convergence
%       maxiter = Maximum number of iterations
%       delta   = Weight value for normalization
%       mu      = Regularization weight in abundance estimation
%       t       = Variance for lambda estimation
%       alpha   = Alpha cero for lambda estimation
%
% Output Argument
%   G   = End-members estimated
%   A   = Abundance estimated
%
%   -   Shao, J. Lan, Y. Zhang, and J. Zou, “Spectral unmixing of hyper-spectral
%       remote sensing imagery via preserving the intrinsic structureinvariant,”
%       Sensors (Switzerland), vol. 18, no. 10, 2018.
%
% Ines A. Cruz-Guerrero
% May/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Go = VCA(X,'Endmembers',m,'SNR',1,'verbose','no');
Ao = max((Go'*Go)\Go'*X,0);


[W,D]=CalMatW2(X,para);
% D=0;
% W=0;

[l,n]=size(X);
iter=0;
while (sqrt(norm(X-Go*Ao)^2/l)/n > para.tau) && (para.maxiter > iter)
    lambda=para.alpha*exp(-iter/para.t);
    
    Go = (Go.*(X*Ao'))./(Go*(Ao*Ao'));
    Xb = X; Xb(l+1,:) = para.delta*ones(1,n);
    Gb = Go; Gb(l+1,:) = para.delta*ones(1,m);
    
    Ao = Ao.*(Gb'*Xb + para.mu*Ao*W)./((Gb'*Gb)*Ao + 0.5*lambda*sqrt(Ao) + para.mu*Ao*D);
    
    Ao(isnan(Ao))=0;
    
    iter=iter+1;
end
G=Go;
A=Ao;
end


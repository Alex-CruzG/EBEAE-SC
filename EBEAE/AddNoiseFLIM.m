function [Y] = AddNoiseFLIM(Y0,SNR,PSNR)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % [Y] = AddNoiseFLIM(Y0,SNR,PSNR)
    %
    % Function to add a certain level of Gaussian noise and Shot noise in decibels
    %
    % Input Arguments
    %   Yo  = matrix of measurements (MxN)
    %   SNR = Gaussian noise level 
    %   PSNR= Shot noise level
    %
    % Output Argument
    %   Y   = matrix of measurements with noise
    %
    % Ines A. Cruz-Guerrero
    % Mayo/2021
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    r = size(Y0,1);
    c = size(Y0,2);
    L = size(Y0,3);
    Y = zeros(r,c,L);
    y = zeros(L,1);
    for i = 1:r
        for j = 1:c
            y(:,1) = Y0(i,j,:);
            if SNR ~= 0
                sigmay = sqrt((y'*y)/(10^(SNR/10)));
                yy1 = sigmay*randn(L,1);
            else
                yy1 = zeros(L,1);
            end
            if PSNR ~= 0
                sigmay = sqrt(max(y)/(10^(PSNR/10)));
                yy2 = sigmay*randn(L,1).*sqrt(abs(y));
            else
                yy2 = zeros(L,1);
            end
            Y(i,j,:) = y + yy1 + yy2;
        end
    end

end
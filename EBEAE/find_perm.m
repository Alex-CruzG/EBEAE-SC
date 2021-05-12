function [A_hat, P_hat, nmse_a, nmse_p] = find_perm(A_true,P_true, A_hat, P_hat)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % [A_hat, P_hat, nmse_a, nmse_p] = find_perm(A_true, P_true, A_hat, P_hat)
    %
    % Function to find the permutation of the estimates that matches the order 
    % of real abundances and real end-members
    %
    % Input Arguments
    %   A_true  = True abundances
    %   P_true  = True end-members
    %   A_hat   = Estimate abundances
    %   P_hat   = Estimate end-members
    % Output Argument
    %   A_true  = Estimated abundances ordered
    %   P_hat   = Estimated end-members ordered
    %   nmse_a  = Error mse of abundance
    %   nmse_p  = Error mse of end-members
    %
    % Ines A. Cruz-Guerrero
    % Mayo/2021
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P = size(A_true,1);
    ords = perms(1:P);
    n = size(ords,1);
    errs = 100*ones(n,1);

    for idx = 1:n
        A = A_hat(ords(idx,:),:);
        errs(idx) = norm(A_true-A, 'fro')/norm(A_true, 'fro');
    end

    [nmse_a, I] = min(errs);
    A_hat = A_hat(ords(I,:), :);
    P_hat = P_hat(:, ords(I,:));
    nmse_p=norm(P_true-P_hat, 'fro')/norm(P_true, 'fro');
end
function [count, A, S]= glnmf(X, P, A_init, S_init, para)

[~, N] = size(X);

delta = 15;
k = 5;
sigma = 5;

try
    load('nmf_graph.mat', 'W', 'd');
    if size(X,2)~=size(d,1) && (size(X,2)~=size(W,1) && size(X,2)~=size(W,2))
        iis = zeros(k*N,1);
        jjs = zeros(k*N,1);
        vvs = zeros(k*N,1);
        for col = 1:N
            w = pdist2(X(:, col)', X','squaredeuclidean');
            [~, I] = sort(w,'ascend'); 
            iis(k*(col-1)+1:k*col) = I(2:k+1); % no self loops 
            jjs(k*(col-1)+1:k*col) = col*ones(k,1);
            vvs(k*(col-1)+1:k*col) = exp(-w(I(2:k+1))/sigma);
        end
        W = sparse(iis, jjs, vvs, N, N);
        d = sum(W,2);
        save('nmf_graph.mat','W', 'd');
    end
catch
    iis = zeros(k*N,1);
    jjs = zeros(k*N,1);
    vvs = zeros(k*N,1);
    for col = 1:N
        w = pdist2(X(:, col)', X','squaredeuclidean');
        [~, I] = sort(w,'ascend'); 
        iis(k*(col-1)+1:k*col) = I(2:k+1); % no self loops 
        jjs(k*(col-1)+1:k*col) = col*ones(k,1);
        vvs(k*(col-1)+1:k*col) = exp(-w(I(2:k+1))/sigma);
    end
    W = sparse(iis, jjs, vvs, N, N);
    d = sum(W,2);
    save('nmf_graph.mat','W', 'd');
end
D = sparse(diag(d));
A = A_init;
S = S_init;
Xbar = [X; delta*ones(1,N)];

count = 0;
while (norm(X-A*S, 'fro')/norm(X, 'fro') > para.tol) && (count < para.itermax)
    
    A = A.*(X*S'./(A*(S*S')));

    Abar = [A; delta*ones(1,P)];
    
    S = max(10^(-8), S);
    mask = S>10^(-4);

    S = S.*(Abar'*Xbar+para.mu*S*W)./(Abar'*Abar*S+0.5*para.lambda*(S.^(-0.5)).*mask+para.mu*S*D);
    count = count +1;
end
end
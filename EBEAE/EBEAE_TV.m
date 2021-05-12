function [P,A,Wm,Yh]=EBEAE_TV(Yo,n,parameters,sc,Po,oae)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   [P,A,Wm,Yh]=EBEAE_TV(Y,n,parameters,Po,oae)
    %
    %   Extended Blind End-member and Abundance Estimation with Total Variation
    %   Estimation of Optimal End-members and Abundances in Linear Mixture Model
    %
    %
    % Input Arguments
    %   Y = matrix of measurements (MxN)
    %   n = order of linear mixture model
    %   parameters = 8x1 vector of hyper-parameters in BEAE methodology
    %              = [initicond rho lambda epsilon maxiter parallel...
    %                   normalization display]
    %       initcond = initialization of end-members matrix {1,2,3}
    %                                 (1) Maximum cosine difference from mean
    %                                      measurement (default)
    %                                 (2) Maximum and minimum energy, and
    %                                      largest distance from them
    %                                 (3) PCA selection + Rectified Linear Unit
    %                                 (4) ICA selection (FOBI) + Rectified
    %                                 Linear Unit
    %       rho = regularization weight in end-member estimation 
    %             (default rho=0.1);
    %       lambda = entropy weight in abundance estimation \in [0,1) 
    %                (default lambda=0);
    %       epsilon = threshold for convergence in ALS method 
    %                 (default epsilon=1e-3); 
    %       maxiter = maximum number of iterations in ALS method
    %                 (default maxiter=20);
    %       parallel = implement parallel computation of abundances (0 -> NO or 1 -> YES)
    %                  (default parallel=0);
    %       normalization = normalization of estimated end-members (0 -> NO or 1 ->YES)
    %                       (default normalization=1);
    %       display = show progress of iterative optimization process (0 -> NO or 1 -> YES)
    %                 (default display=0);
    %
    %   sc  =   vector with the information of spatial coherence
    %       =   [mu, nu, tau, dim(1), dim(2)]
    %           mu  = regularization term of spatial coherence
    %           nu, tau = split Bregman regularization variables
    %           dim = spatial dimention of image to analize
    %
    %   Po = initial end-member matrix (Mxn)
    %   oae = only optimal abundance estimation with Po (0 -> NO or 1 -> YES)
    %         (default oae = 0)
    %
    % Output Arguments
    %
    %   P   = matrix of end-members (Mxn)
    %   A   = abudances matrix normalized (nxN)
    %   Wm  = noise-free abundance matrix (nxN)
    %   Yh = estimated matrix of measurements (MxN)
    %
    % Ines A. Cruz-Guerrero
    % Mayo/2021
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default hyper-parameters of BEAE algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global NUMERROR

    initcond=1;
    rho=0.1;
    lambda=0;
    epsilon=1e-3;
    maxiter=20;
    parallel=0;
    normalization=1;
    display=0;
    NUMERROR=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check concistency of input arguments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin~=6
        oae=0;
    end
    if nargin==0
        disp('The measurement matrix Y has to be used as argument!!');
        return;
    elseif nargin==1
        n=2;
    end
    if nargin==4 || nargin==5 || nargin==6
        if length(parameters)~= 8
            disp('The length of parameters vector is not 8 !!');
            disp('Default values of hyper-parameters are used instead');
        else
            initcond=round(parameters(1));
            rho=parameters(2);
            lambda=parameters(3);
            epsilon=parameters(4);
            maxiter=parameters(5);
%             downsampling=parameters(6);
            parallel=parameters(6);
            normalization=parameters(7);
            display=parameters(8);
            if initcond~=1 && initcond~=2 && initcond~=3 && initcond~=4
                disp('The initialization procedure of end-members matrix is 1,2,3 or 4!');
                disp('The default value is considered!');
                initcond=1;
            end
            if rho<0
                disp('The regularization weight rho cannot be negative');
                disp('The default value is considered!');
                rho=0.1;
            end
            if lambda<0 || lambda>=1
                disp('The entropy weight lambda is limited to [0,1)');
                disp('The default value is considered!');
                lambda=0;
            end
            if epsilon<0 || epsilon>0.5
                disp('The threshold epsilon cannot be negative or >0.5');
                disp('The default value is considered!');
                epsilon=1e-3;
            end
            if maxiter<0 && maxiter<100
                disp('The upper bound maxiter cannot be negative or >100');
                disp('The default value is considered!');
                maxiter=20;
            end
            if parallel~=0 && parallel~=1
                disp('The parallelization parameter is 0 or 1');
                disp('The default value is considered!');
                parallel=0;
            end
            if normalization~=0 && normalization~=1
                disp('The normalization parameter is 0 or 1');
                disp('The default value is considered!');
                normalization=1;
            end
            if display~=0 && display~=1
                disp('The display parameter is 0 or 1');
                disp('The default value is considered!');
                display=0;
            end
        end
        if n<2
            disp('The order of the linear mixture model has to greater than 2!');
            disp('The default value n=2 is considered!');
            n=2;
        end
    end
    if ((nargin==4||nargin==5 || nargin==6) && ~isempty(sc))
        if length(sc)~=5
            disp('The parameters of spatial coherence are: lambda = regularization weith, nu, tau, and dim = dimentions of the image to analize')
            disp('The evaluation will be carried out without considering spatial coherence!');
            sc=[];
        else
            if ~(sc(1)>=0 && sc(1)<=1)
                disp('The parameter mu must be <=1');
                sc(1)=1;
            end
            if ~(sc(2)>=0 && sc(2)<=1)
                disp('The parameter nu must be <=1');
                sc(2)=1;
            end
            if ~(sc(3)>=0 && sc(3)<=1)
                disp('The parameter tau must be <=1');
                sc(3)=1;
            end
            if (sc(4)*sc(5) ~= size(Yo,2))
                disp('The spatial dimensions do not match the number of measurements in the Y matrix');
                disp('The evaluation will be carried out without considering spatial coherence!');
                sc=[];
            end
        end
    else
        disp('The evaluation will be carried out without considering spatial coherence!');
    end
    if nargin==5 || nargin==6
        if ~ismatrix(Po)
            disp('The initial end-members Po must be a matrix !!');
            disp('The initialization is considered by the maximum cosine difference from mean measurement');
            initcond=1;
        else
            if size(Po,1)==size(Yo,1) && size(Po,2)==n
                initcond=0;
            else
                disp('The size of Po must be Mxn!!');
                disp('The initialization is considered based on the input dataset');
                initcond=1;
            end
        end
    end
    if nargin==6
        if oae~=0 && oae~=1
            disp('The assignment of oae is incorrect!!');
            disp('The initial end-members Po will be improved iteratively from a selected sample');
            oae=0;
        elseif oae==1 && initcond~=0
            disp('The initial end-members Po is not defined properly!');
            disp('Po will be improved iteratively from a selected sample');
            oae=0;
        end
    end
    if nargin>7
        disp('The number of input arguments is 6 maximum');
        disp('Please check the help documentation');
        return;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Random downsampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~ismatrix(Yo)
        disp('The measurements matrix Y has to be a matrix');
        return;
    end
    [M,No]=size(Yo);
    if M>No
        disp('The number of spatial measurements has to be larger to the number of time samples!');
        return;
    end
    
    N=No;
    % I=1:No;
    % N=round(No*(1-downsampling));
    % Is=randperm(No,N);
    Y=Yo;%(:,Is);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if normalization==1
        mYm=sum(Y,1);
        mYmo=sum(Yo,1);
    else
        mYm=ones(1,N);
        mYmo=ones(1,No);
    end
    Ym=Y./repmat(mYm,[M 1]);
    Ymo=Yo./repmat(mYmo,[M 1]);
    NYm=norm(Ym,'fro');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Selection of Initial End-members Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if initcond==1 || initcond==2
        if initcond==1
            Po=zeros(M,n);
            index=1;
            pmax=mean(Yo,2);
            Yt=Yo;
            Po(:,index)=pmax;
        elseif initcond==2
            index=1;
            Y1m=sum(abs(Yo),1);
            [~,Imax]=max(Y1m);
            [~,Imin]=min(Y1m);
            pmax=Yo(:,Imax);
            pmin=Yo(:,Imin);
            K=size(Yo,2);
            II=1:K;
            Yt=Yo(:,setdiff(II,[Imax Imin]));
            Po(:,index)=pmax;
            index=index+1;
            Po(:,index)=pmin;
        end
        while index<n
            ymax=zeros(1,index);
            Imax=zeros(1,index);
            for i=1:index
                e1m=sum(Yt.*repmat(Po(:,i),1,size(Yt,2)),1)./sqrt(sum(Yt.^2,1))./sqrt(sum(Po(:,i).^2,1));
                [ymax(i),Imax(i)]=min(abs(e1m));
            end
            [~,Immax]=min(ymax);
            IImax=Imax(Immax);
            pmax=Yt(:,IImax);
            index=index+1;
            Po(:,index)=pmax;
            II=1:size(Yt,2);
            Yt=Yt(:,setdiff(II,IImax));
        end
    elseif initcond==3
        [~,~,VV]=svd(Ym',0);
        W=VV(:,1:n);
        Po=W.*repmat(sign(W'*ones(M,1))',M,1); 
    elseif initcond==4
        Yom=mean(Ym,2);
        Yon = Ym - repmat(Yom,1,N);
        [~,S,VV]=svd(Yon',0);
        Yo_w= pinv(sqrtm(S))*VV'*Ym; 
        [V,~,~] = svd((repmat(sum(Yo_w.*Yo_w,1),M,1).*Yo_w)*Yo_w');
        W=VV*sqrtm(S)*V(1:n,:)'; 
        Po=W.*repmat(sign(W'*ones(M,1))',M,1);
    end    

    Po(Po<0)=0;
    Po(isnan(Po))=0;
    Po(isinf(Po))=0;

    if normalization==1
        mPo=sum(Po,1);
        P=Po./repmat(mPo,[M 1]);
    else
        P=Po;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Alternated Least Squares Procedure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iter=1;
    J=1e5;
    Jp=1e6;
    W=zeros(n,size(Ym,2));
    Wm=[];
    tic;
    if display==1
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp('EBEAE Linear Unmixing');
            disp(['Model Order =' num2str(n)]);
            if oae==1
                disp('Only the abundances are estimated from Po');
            elseif oae==0 && initcond==0
                disp('The end-members matrix is initialized externally by matrix Po');
            elseif oae==0 && initcond==1
                disp('Po is constructed based on the maximum cosine difference from mean measurement'); 
            elseif oae==0 && initcond==2
                disp('Po is constructed based on the maximum and minimum energy, and largest difference from them');
            elseif oae==0 && initcond==3
                disp('Po is constructed based on the PCA selection + Rectified Linear Unit');
            elseif oae==0 && initcond==4
                disp('Po is constructed based on the ICA selection (FOBI) + Rectified Linear Unit');
            end
    end
    if ~isempty(sc)
        while (Jp-J)/Jp >= epsilon && iter < maxiter && oae==0 && NUMERROR==0

            Am = abundanceSC(Ym,P,W,parallel,sc(1));    % abundance estimation modification
            Pp=P;
            if NUMERROR==0
                P = endmember(Ym,Am,rho,normalization);
            end
            Jp=J;
            W=denoise_abundance(Am,sc(1),sc(2),sc(3),sc(4),sc(5),10);   W=normalize(W,'norm',1);
            J=norm(Ym-P*Am,'fro');
            if J > Jp
                P=Pp; break;
            end
            if display ==1
                disp(['Number of iteration =' num2str(iter)]);
                disp(['Percentage Estimation Error =' num2str(100*J/NYm) '%']);
            end
            iter=iter+1;
        end
        if NUMERROR==0
            Am = abundanceSC(Ymo,P,W,parallel,sc(1));
            W=denoise_abundance(Am,sc(1),sc(2),sc(3),sc(4),sc(5),10);   W=normalize(W,'norm',1);
            ElapTime=toc;
            if display ==1
                disp(['Elapsep Time =' num2str(ElapTime)]);
            end
            A=normalize(Am,'norm',1);
            Wm=W;
            Yh=P*Wm;
        else
            disp('Please revise the problem formulation, not reliable results');
            P=[];
            A=[];
            Yh=[];
            Wm=[];
        end
    else
        while (Jp-J)/Jp >= epsilon && iter < maxiter && oae==0 && NUMERROR==0
            Am = abundance(Ym,P,lambda,parallel);
            Pp=P;
            if NUMERROR==0
                P = endmember(Ym,Am,rho,normalization); 
            end
            Jp=J;
            J=norm(Ym-P*Am,'fro');
            if J > Jp
                P=Pp; break;
            end
            if display ==1
                disp(['Number of iteration =' num2str(iter)]);
                disp(['Percentage Estimation Error =' num2str(100*J/NYm) '%']);
            end
            iter=iter+1;
        end
        if NUMERROR==0
            Am = abundance(Ymo,P,lambda,parallel); 
            ElapTime=toc;
            if display ==1
                disp(['Elapsep Time =' num2str(ElapTime)]);
            end
            A=normalize(Am,'norm',1);
            Yh=P*A;
        else
            disp('Please revise the problem formulation, not reliable results');
            P=[];
            A=[];
            Yh=[];
        end
    end
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
end

function A = abundanceSC(Y,P,W,parallel, mu)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % A = abundanceSC(Y,P,W,parallel,mu)
    %
    % Estimation of Optimal Abundances with spatial coherence in Linear Mixture Model
    %
    % Input Arguments
    %   Y           = matrix of measurements
    %   P           = matrix of end-members
    %   parallel    = implementation in parallel of the estimation
    %   mu          = regularization term of spatial coherence \in [0,1]
    %
    % Output Argument
    %   A           = abundance matrix (nxN)
    %
    % Ines A. Cruz-Guerrero
    % Mayo/2021
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check arguments dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global NUMERROR

    [M,N]=size(Y);
    n=size(P,2);
    A=zeros(n,N);

    if size(P,1) ~= M
        disp('ERROR: the number of rows in Y and P does not match');
        NUMERROR=1;
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute fixed vectors and matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    em=eye(n);
    c = ones(n,1);
    d = 1;
    Go=P'*P;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start Computation of Abundances
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parallel==1
        parfor k=1:N
            yk=Y(:,k);
            byk=yk'*yk;
            bk=P'*yk+byk*mu*W(:,k);
            G=Go+em*mu*byk;
            Gi=em/G;
            T1=Gi*c;
            T2=c'*T1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute Optimal Unconstrained Solution
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dk=(bk'*T1-1)/T2;
            ak = Gi*(bk-dk*c);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check for Negative Elements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(sum(ak>=0) ~=n)
                Iset = zeros(1,n);
                while(sum(ak<0) ~= 0)    
                    Iset(ak<0) = 1;
                    L = length(find(Iset));

                    Q = n+1+L;
                    Gamma = zeros(Q);
                    Beta = zeros(Q,1);

                    Gamma(1:n,1:n) = G/byk;
                    Gamma(1:n,n+1) = c;
                    Gamma(n+1,1:n) = c';

                    cont = 0;
                    for i = 1:n
                        if(Iset(i)~= 0)
                            cont = cont + 1;
                            ind = i; 
                            Gamma(ind,n+1+cont) = 1;
                            Gamma(n+1+cont,ind) = 1;   
                        end
                    end
                    Beta(1:n) = bk/byk;
                    Beta(n+1) = d;
                    delta = Gamma\Beta;
                    ak = delta(1:n);
                    ak(abs(ak)<1e-9) = 0;
                end    
            end
            A(:,k) = single(ak); 
        end
    else
        for k=1:N
            yk=Y(:,k);
            byk=yk'*yk;
            bk=P'*yk+byk*mu*W(:,k);
            G=Go+em*mu*byk;
            Gi=em/G;
            T1=Gi*c;
            T2=c'*T1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute Optimal Unconstrained Solution
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dk=(bk'*T1-1)/T2;
            ak = Gi*(bk-dk*c);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check for Negative Elements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(sum(ak>=0) ~=n)
                Iset = zeros(1,n);
                while(sum(ak<0) ~= 0)    
                    Iset(ak<0) = 1;
                    L = length(find(Iset));

                    Q = n+1+L;
                    Gamma = zeros(Q);
                    Beta = zeros(Q,1);

                    Gamma(1:n,1:n) = G/byk;
                    Gamma(1:n,n+1) = c;
                    Gamma(n+1,1:n) = c';

                    cont = 0;
                    for i = 1:n
                        if(Iset(i)~= 0)
                            cont = cont + 1;
                            ind = i; 
                            Gamma(ind,n+1+cont) = 1;
                            Gamma(n+1+cont,ind) = 1;   
                        end
                    end
                    Beta(1:n) = bk/byk;
                    Beta(n+1) = d;
                    delta = Gamma\Beta;
                    ak = delta(1:n);
                    ak(abs(ak)<1e-9) = 0;
                end    
            end
            A(:,k) = single(ak); 
        end
    end
end

function A = abundance(Y,P,lambda,parallel)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % A = abundance(Y,P,lambda,parallel)
    %
    % Estimation of Optimal Abundances in Linear Mixture Model
    %
    % Input Arguments
    %   Y       = matrix of measurements
    %   P       = matrix of end-members
    %   lambda  = entropy weight in abundance estimation \in (0,1)
    %   parallel= implementation in parallel of the estimation
    %
    % Output Argument
    %   A = abundances matrix 
    %
    % Daniel U. Campos-Delgado
    % Oct/2020
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check arguments dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global NUMERROR

    [M,N]=size(Y);
    n=size(P,2);
    A=zeros(n,N);
    if size(P,1) ~= M
        disp('ERROR: the number of rows in Y and P does not match');
        NUMERROR=1;
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute fixed vectors and matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c = ones(n,1);
    d = 1;
    Go=P'*P;
    lmin=min(eig(Go));
    G=Go-eye(n)*lmin*lambda;
    while rcond(G)<1e-6
        lambda=lambda/2;
        G=Go-eye(n)*lmin*lambda;
        if lambda<1e-6
           disp('Unstable numerical results in abundances estimation, update lambda!!');
           NUMERROR=1;
           return;
        end
    end
    Gi=eye(n)/G;
    T1=Gi*c;
    T2=c'*T1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start Computation of Abundances
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parallel==1
        parfor k=1:N
            yk=Y(:,k);
            bk=P'*yk;
            byk=yk'*yk;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute Optimal Unconstrained Solution
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dk=(bk'*T1-1)/T2;
            ak = Gi*(bk-dk*c);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check for Negative Elements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(sum(ak>=0) ~=n)
                Iset = zeros(1,n);
                while(sum(ak<0) ~= 0)    
                    Iset(ak<0) = 1;
                    L = length(find(Iset));

                    Q = n+1+L;
                    Gamma = zeros(Q);
                    Beta = zeros(Q,1);

                    Gamma(1:n,1:n) = G/byk;
                    Gamma(1:n,n+1) = c;
                    Gamma(n+1,1:n) = c';

                    cont = 0;
                    for i = 1:n
                        if(Iset(i)~= 0)
                            cont = cont + 1;
                            ind = i; 
                            Gamma(ind,n+1+cont) = 1;
                            Gamma(n+1+cont,ind) = 1;   
                        end
                    end
                    Beta(1:n) = bk/byk;
                    Beta(n+1) = d;
                    delta = Gamma\Beta;
                    ak = delta(1:n);
                    ak(abs(ak)<1e-9) = 0;
                end    
            end
            A(:,k) = single(ak); 
        end
    else
        for k=1:N
            yk=Y(:,k);
            byk=yk'*yk;
            bk=P'*yk;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute Optimal Unconstrained Solution
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dk=(bk'*T1-1)/T2;
            ak = Gi*(bk-dk*c);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check for Negative Elements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(sum(ak>=0) ~=n)
                Iset = zeros(1,n);
                while(sum(ak<0) ~= 0)    
                    Iset(ak<0) = 1;
                    L = length(find(Iset));
                    
                    Q = n+1+L;
                    Gamma = zeros(Q);
                    Beta = zeros(Q,1);

                    Gamma(1:n,1:n) = G/byk;
                    Gamma(1:n,n+1) = c;
                    Gamma(n+1,1:n) = c';

                    cont = 0;
                    for i = 1:n
                        if(Iset(i)~= 0)
                            cont = cont + 1;
                            ind = i; 
                            Gamma(ind,n+1+cont) = 1;
                            Gamma(n+1+cont,ind) = 1;   
                        end
                    end

                    Beta(1:n) = bk/byk;
                    Beta(n+1) = d;
                    delta = Gamma\Beta;
                    ak = delta(1:n);
                    ak(abs(ak)<1e-9) = 0;
                end    
            end
            A(:,k) = single(ak); 
        end
    end
end

function P = endmember(Y,A,rho,normalization)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %  P = endmember(Y,A,rho,normalization)
    %
    % Estimation of Optimal End-members in Linear Mixture Model
    %
    % Input Arguments
    % Y = Matrix of measurements
    % A =  Matrix of abundances
    % rho = Weighting factor of regularization term
    % normalization = normalization of estimated profiles (0=NO or 1=YES)
    % 
    % Output Arguments
    % P = Matrix of end-members
    %
    % Daniel U. Campos-Delgado
    % Oct/2020
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check arguments dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global NUMERROR

    [n,N]=size(A);
    [M,K]=size(Y);
    P=zeros(M,n);
    R=sum(n - (1:(n-1)));
    W=repmat((1./K./sum(Y.^2,1))',1,n);

    if size(Y,2) ~= N
        disp('ERROR: the number of columns in Y and A does not match');
        NUMERROR=1;
        return;
    end
    O = single(n*eye(n) - ones(n,n));   
    n1 = single(ones(n,1));
    m1 = single(ones(M,1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct Optimal End-members Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T0=A*(W.*A') + rho*O/R;
    while rcond(T0)<1e-6
       rho=rho/10;
       T0=A*(W.*A') + rho*O/R;
       if rho<1e-6
           disp('Unstable numerical results in end-members estimation, update rho!!');
           NUMERROR=1;
           return;
       end
    end
    V = eye(n)/T0;
    T2 = Y*(W.*A')*V;
    if normalization == 1
        T1 = single(eye(M) - (1/M)*(m1*m1'));
        T3 = (1/M)*m1*n1';
        P_est = T1*T2 + T3;
    else
        P_est=T2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate and Project Negative Elements 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_est(P_est<0) = 0;
    P_est(isnan(P_est))=0;
    P_est(isinf(P_est))=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalize Optimal Solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if normalization==1
        Psum=sum(P_est,1);
        P=P_est./repmat(Psum,M,1);
    else 
        P=P_est;
    end
end

function dW=denoise_abundance(A,mu,nu,tau,m,n,maxiter)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % dW = denoise_abundance(A,mu,nu,tau,m,n,maxiter)
    %
    % Estimation of Optimal noise-free abundance with total variation theory in Linear Mixture Model
    %
    % Input Arguments
    %   A       = matrix of abundances
    %   mu      = regularization term of spatial coherence \in [0,1]
    %   nu, tau = regularization term of split Bragman \in [0,1]
    %   m, n    = vertical and horizontal spatial dimensions
    %   maxiter = regularization term of spatial coherence \in [0,1]
    %
    % Output Argument
    %   dW      = noise-free abundance matrix (nxN)
    %
    % Ines A. Cruz-Guerrero
    % Mayo/2021
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition of the soft-thresholding function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SoftTh=@(B,lambda)  sign(B).*max(0,abs(B)-(lambda));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization of variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dim,mn]=size(A);
    b1=zeros(mn,1);
    W=A';    b2=b1;  p=b1;   q=b1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition of derivatives
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dh=spdiags([-ones(n,1) ones(n,1)],[0 1],n,n);   Dh(n,:) = 0;    Dh = kron(Dh,speye(m));
    Dv=spdiags([-ones(m,1) ones(m,1)],[0 1],m,m);   Dv(m,:) = 0;    Dv = kron(speye(n),Dv);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start computation of noise-free abundance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:maxiter
        for j=1:dim
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Least squares stage
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Ay=mu*A(j,:)' + nu*Dh'*(p-b1)+nu*Dv'*(q-b2);
            [W(:,j),~]=lsqr(@afun,Ay,1e-15,10,[],[],W(:,j));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update soft-thresholding function
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p=SoftTh(Dh*W(:,j)+b1,tau/nu);
            q=SoftTh(Dv*W(:,j)+b2,tau/nu);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update of Bregman variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            b1=b1+Dh*W(:,j)-p;
            b2=b2+Dv*W(:,j)-q;
        end
    end
    W(W<0)=0;
    dW=W';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Handle for least squares
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y = afun(W,str)
        tempval= nu*((Dh'*(Dh*W))+(Dv'*(Dv*W)))+ mu*W;
        switch str
            case 'transp'
                y = tempval;
            case 'notransp'
                y = tempval;
        end
    end
end

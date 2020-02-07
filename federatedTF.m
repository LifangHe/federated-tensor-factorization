function [rmses, tElapsed, communication, tElapsedIter, A, u]= federatedTF(maxiter, Xs, XsSum, dim, rank, K, guide, cutoffs, omega, gamma, paraUpdateRatio,lambda_q )
%% TRIP until max iter
% Yejin Kim

%% Server: set algorithm parameters
rmses=zeros(1, maxiter);
%rmses_w_ori=zeros(1, maxiter);
communication=0;
tElapsedIter=zeros(1, maxiter);

lambda_a = 0; % weight for guide
%lambda_q = 1e-2; % weight for pairwise 
sparsity = [0.01, 0.1, 0.1]; % sparse threshold

mu = 1e-10;
eta = 1e-10;
%omega=10; 
%gamma=1e-13;

mumax = 1e10;
etamax=1e10;
omegamax = 1e10;
gammamax = 1e10;

%paraUpdateRatio = 1.15; 


%% Set clients workspace

for k=1:K
    client(k).Hi={}; 
    client(k).gi={};
    
    %set initial parameters

    client(k).omega=omega;
    client(k).gamma=gamma;
    client(k).X=Xs{k};
    client(k).sparsity= [0.01, 0.1, 0.1];

    client(k).omegamax = 1e10;
    client(k).gammamax = 1e10;
    client(k).paraUpdateRatio = paraUpdateRatio;

end


%% Server: compute statistics for X and guide

%T_ori=ktensor(A_ori);
%C_ori=ktensor(u_ori);

numGuide = size(guide, 2);

guide = [guide, zeros(dim(2), rank - numGuide)];

W = [eye(numGuide), zeros(numGuide, rank - numGuide); zeros(rank - numGuide,rank)]; % K-by-K symetric matrix


%% Server: initialization

% interaction tensor 

A = cell(1,3);
for n = 2:3 % n=1 initialize at institutions
    A{n} = rand(dim(n),rank);
    A{n}(A{n} <= sparsity(n)) = 0; 
end


 B = A;

% bias tensor

u = cell(1,3);

for n = 2:3 %n=1 initialize at institutions
    u{n} = rand(dim(n),1)/10;
end

v = u;

% Lagrange multipliers
Y = cell(1,3);
for n = 1:3
    Y{n} = zeros(dim(n),rank);
end

p = cell(1,3);
for n = 1:3
    p{n} = zeros(dim(n),1);
end

Hi=cell(1, K);
for k=1:K
    Hi{k}=cell(1,3);
    for n=2:3
        Hi{k}{n}=zeros(dim(n), rank);
    end
end


gi=cell(1, K);
for k=1:K
    gi{k}=cell(1,3);
    for n=2:3
        gi{k}{n}=zeros(dim(n), 1);
    end
end


%% Hospitals: compute statistics for X / initialization
% client= a structure contains client i's information, dim, A, u, normX.

for k=1:K
    client(k).dim=size(client(k).X); %data
    client(k).normX=norm(client(k).X);
    
    % dist interaction tensor
    
    client(k).Ai=cell(1,3); 
    % n = 1
    client(k).Ai{1}=zeros(dim(1),rank);
    client(k).Ai{1}(cutoffs(k):(cutoffs(k+1)-1),:)=rand(cutoffs(k+1)-cutoffs(k), rank);
    client(k).Ai{1}(client(k).Ai{1}<=sparsity(1)) = 0;

    for n = 2:3       
        client(k).Ai{n} = A{n};
    end   
    
    % dist bias tensor
    
    client(k).ui= cell(1,3);
    % n =1
    client(k).ui{1}=zeros(dim(1),1);
    client(k).ui{1}(cutoffs(k):(cutoffs(k+1)-1),:)=rand(cutoffs(k+1)-cutoffs(k), 1);
    
    for n = 2:3
        client(k).ui{n} = u{n};
    end
    
    client(k).Hi=Hi{k};
    client(k).gi=gi{k};
    
    client(k).A=A; %for n=2 3
    client(k).u=u; %for n=2 3
    
end

% communication A u
Abytes=whos('A');
ubytes=whos('u');
communication=communication+Abytes.bytes*K;
communication=communication+ubytes.bytes*K;

%% main loop

dist_time=zeros(maxiter*3*3, K); % per iteration, 3 modes, 3 distirubted computation (initial, A, u)
dt=0;
tStart=tic;
for iter = 1: maxiter

    %client iter 1
    dt=dt+1;
    for k=1:K
        tStartD=tic;
        
        client(k).Ri=plusKtensor(client(k).X, -ktensor(client(k).ui));
        client(k).Ei=plusKtensor(client(k).X, -ktensor(client(k).Ai));  
        %client(k).Ri = client(k).X - tensor(ktensor(client(k).ui));
        %client(k).Ei = client(k).X - tensor(ktensor(client(k).Ai));

        
        id1 = find(client(k).Ri<0);
        id2 = find(client(k).Ei<0);
        if ~isempty(id1)
            client(k).Ri(id1) = 0;
        end
        if ~isempty(id2)
            client(k).Ei(id2) = 0;
        end
        dist_time(dt, k)=toc(tStartD);
    end
   
    
    for n = 1: 3
        if (n ==1)
            % update Ai
            
            Ai=cell(1,K);
            
            %client iter 2
            dt=dt+1;
            for k=1:K
                tStartD=tic;
                
                piitpii=ones(rank,rank);
                for nn=[1:n-1, n+1:3]
                    piitpii=piitpii .*( client(k).Ai{nn}' * client(k).Ai{nn} );% compute \Pi^t\Pi in Eq.?
                end
          
                term1 = 2 * mttkrp(client(k).Ri, client(k).Ai,n);
                term2 = 2 * (piitpii) ;
                client(k).Ai{n}= term1/term2;
                              
                %JUST FOR PERFORMANCE TEST send client(i).Ai{n} to server Ai.
                Ai{k}=client(k).Ai{n};
                dist_time(dt, k)=toc(tStartD);
            end
            
            % JUST FOR PERFORMANCE TEST 
            %update sum(Ai_n)
            sumAi_n=zeros(size(Ai{1}));
            for k=1:K
                sumAi_n=sumAi_n+Ai{k};
            end               
            A{n}=  sumAi_n ;
            
           
            % update ui
            
            ui=cell(1,K);
            
            %client iter 3
            dt=dt+1;
            for k= 1:K
                tStartD=tic;
                LambdaitLambdai = 1;
                for nn = [1:n-1,n+1:3]
                    LambdaitLambdai = LambdaitLambdai .* (client(k).ui{nn}' * client(k).ui{nn});
                end
                
                term1 = 2 * mttkrp(client(k).Ei, client(k).ui, n);
                term2 = 2 * LambdaitLambdai ; 
                client(k).ui{n}= term1/term2;

                %JUST FOR PERFORMANCE TEST. send to server
                ui{k}=client(k).ui{n};               
                dist_time(dt, k)=toc(tStartD);
            end
            
            %JUST FOR PERFORMANCE TEST
            % update sum(ui_n)
            sumui_n=zeros(size(ui{1}));
            for k=1:K
                sumui_n=sumui_n+ui{k};
            end
            u{n}=sumui_n;
            
            
        else %n = 2 , 3


            Ai=cell(1,K);
            Hi=cell(1,K);
            
            %client iter 1
            dt=dt+1;
            for k=1:K
                tStartD=tic;            
                
                % update Ai
                
                piitpii=ones(rank,rank);
                for nn=[1:n-1, n+1:3]
                    piitpii=piitpii .*( client(k).Ai{nn}' * client(k).Ai{nn} );% compute \Pi^t\Pi in Eq.?
                end
                
                term1 = 2 * mttkrp(client(k).Ri, client(k).Ai,n) + client(k).omega * client(k).A{n} + client(k).Hi{n} ;               
                term2 = 2 * (piitpii) + client(k).omega * eye(rank);
                client(k).Ai{n}= term1/term2;
                
                
                % update Hi
                client(k).Hi{n} = client(k).Hi{n} + omega * (client(k).A{n} - client(k).Ai{n});              
                
                % Send client(i).Ai Hi to server
                Ai{k}=client(k).Ai{n};
                Hi{k}=client(k).Hi{n};
                
                dist_time(dt, k)=toc(tStartD);
            end
            
            % Communication
            
            AiBytes=whos('Ai');
            HiBytes=whos('Hi');
            communication=communication+AiBytes.bytes;
            communication=communication+HiBytes.bytes;
            
            % update sum(Ai) sum(Hi)
            
            sumAi_n=zeros(size(Ai{1}));
            for k=1:K
                sumAi_n=sumAi_n+Ai{k};
            end
                 
            sumHi_n=zeros(size(Hi{1}));
            for k=1:K
                sumHi_n=sumHi_n+Hi{k};
            end
            
            
            % update A{n}
            % only mode=2 has guide information and disjoint regularizer.

            if n == 2
                
                a = lambda_q * B{n} * B{n}';
                b = lambda_a * W + mu * eye(rank) + omega * K* eye(rank);                
                c = lambda_a * guide * W  + lambda_q * B{n} + mu * B{n} + Y{n} + omega * sumAi_n - sumHi_n;
                
                try 
                A{n} = sylvester(a,b,c);
                
                %A{n}=lyap(a,b,-c);
                catch
                    %fileName=strcat('UCSD', dataset, 'omega', num2str(omegas(o)), 'gamma', num2str(gammas(g)), 'updateRatio', num2str(updateRatios(u)), 'sylvesterError', '.mat');
                    %save(fileName);
                    %return;
                    disp(strcat('sylvester error:', num2str(iter), ' at ', num2str(omega), ' ', num2str(gamma), ' ', num2str(paraUpdateRatio), '\n'));
                end
                
            else

                A{n}= (mu * B{n} + Y{n} + omega * sumAi_n - sumHi_n )/(mu + K * omega);

            end
            
            
            % update B{n}
            B{n} = A{n} + Y{n}/mu;
            
            % update Y{n}
            Y{n} = Y{n} + mu * (B{n} - A{n});
       
           
            % update ui gi
            
            ui=cell(1,K);
            gi=cell(1,K);
            
            %client iter 2
            dt=dt+1;
            for k= 1:K

                tStartD=tic;
                LambdaitLambdai = 1;
                for nn = [1:n-1,n+1:3]
                    LambdaitLambdai = LambdaitLambdai .* (client(k).ui{nn}' * client(k).ui{nn});
                end
                
                % update ui
                
                term1 = 2 * mttkrp(client(k).Ei, client(k).ui, n) + client(k).gamma * client(k).u{n} + client(k).gi{n};
                term2 = 2 * LambdaitLambdai + client(k).gamma;
                client(k).ui{n}= term1/term2;
                
                % update gi{n}        
                client(k).gi{n}=client(k).gi{n} + gamma * (client(k).u{n} - client(k).ui{n});
                
                % communicate send ui gi to server
                ui{k}=client(k).ui{n};
                gi{k}=client(k).gi{n};
                
                dist_time(dt, k)=toc(tStartD);
            end
            
            % communication
            uiBytes=whos('ui');
            giBytes=whos('gi');
            communication=communication+uiBytes.bytes;
            communication=communication+giBytes.bytes;
            
            % update sum(ui_n) sum(gi_n)
            sumui_n=zeros(size(ui{1}));
            for k=1:K
                sumui_n=sumui_n+ui{k};
            end
            
            sumgi_n=zeros(size(gi{1}));
            for k=1:K
                sumgi_n=sumgi_n+gi{k};
            end
            
            
            % update u{n}
            
            u{n}= (eta * v{n} + p{n} + gamma * sumui_n - sumgi_n)/(eta + K * gamma);

            % update v{n}
            v{n} = u{n} + p{n}/eta;
            
            % update p{n}
            p{n} = p{n} + eta * (v{n} - u{n});
            
            
            % send A{n}, u{n} to clients
            for k=1:K
                client(k).A{n}=A{n};
                client(k).u{n}=u{n};
            end
            Abytes=whos('A');
            ubytes=whos('u');
            communication=communication+Abytes.bytes*K;
            communication=communication+ubytes.bytes*K;
 
            
        end
    end
        
        
        
    
    % update mu, eta
    mu = min(paraUpdateRatio * mu, mumax);
    
    eta = min(paraUpdateRatio * eta, etamax);
   
    omega = min(paraUpdateRatio * omega, omegamax);
    
    gamma = min(paraUpdateRatio * gamma, gammamax);
    
    for k=1:K
        
        
        client(k).omega = min(client(k).paraUpdateRatio * client(k).omega, client(k).omegamax);
        
        client(k).gamma = min(client(k).paraUpdateRatio * client(k).gamma, client(k).gammamax);
    end

    
   
    
    %% Hospitals: compute the fit
    
    T = ktensor(A);
    
    C = ktensor(u);
    
    normresidual=norm(plusKtensor(XsSum, -(T+C)));
    %normresidual = sqrt( normX^2 + norm(T+C)^2 - 2 * innerprod(XsSum,T+C) );
    
    rmses(iter)=normresidual/(sqrt(nnz(XsSum)));
   % rmses_w_ori(iter)= norm(T_ori+C_ori - (T + C))/(sqrt(double(dim(1)*dim(2)*dim(3))));
    %rmses_w_ori(iter)= sqrt(normAuOri^2 + norm(T+C)^2 - 2 * innerprod(T_ori+ C_ori,T+C))/(sqrt(double(dim(1)*dim(2)*dim(3))));
    
    ignoreTime=sum(sum(dist_time))- sum(max(dist_time, [], 2));
    tElapsedIter(iter)=toc(tStart) - ignoreTime;
    
end


%Computing time
tElapsed=toc(tStart);
ignoreTime=sum(sum(dist_time))- sum(max(dist_time, [], 2));
tElapsed=tElapsed-ignoreTime;

communication=communication*2/3; % remove n=1
communication=(communication/(1024*1024))/(15);

%{
%% clean up final results
  T = arrange(ktensor(A));  % columns are normalized        
  C = arrange(ktensor(u));


for n = 1: 3
    T{n}(T{n} <= sparsity(n)) = 0; % sparse factor  
end


if printitn>0
    %Just for measuring performance
    Xsum=tensor(zeros(size(Xs{1})));
    for k=1:K
        Xsum=Xsum+Xs{k};
    end
    
    normresidual = sqrt( norm(Xsum)^2 + norm(T+C)^2 - 2 * innerprod(Xsum,T+C) );
    
   
    fit = 1 - (normresidual / norm(Xsum)^2);
    
    fprintf(' Final fit = %e \n', fit);
end
%}

end
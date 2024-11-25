function [G,e] = APTR_SGD(Y,P_Omega, miu, theta,R,Z_true)
%% fix theta[t] = 0.001 (step size), miu[t] = 0.01 for t = 1,...
% Initialize
G = Initialize(Y,R);
I = size(P_Omega);
d = length(size(Y));
e = zeros(I(d),1);
T = {}; % stores the newly adjusted core tensors 
  

%% Main Algorithm
for t=1:I(d)
    fprintf('\n t = %d ============= \n \n',t)

    Y_t = subsref(Y,substruct('()',[repelem({':'},d-1) t]));
    P_t = subsref(P_Omega,substruct('()',[repelem({':'},d-1) t]));


%% Estimate the temporal TR-core
    A  = TCP(G(1:d-1));
    Yt = Y_t(:);
    A  = modek_unfolding(A,2);
    P  = modek_unfolding(P_Omega,d);
    Pt = P_t(:);
 
    R_N = A'*diag(Pt)*A + 1e-4*eye(R(d-1)*R(d));
    s_N = A'*diag(Pt)*Yt;
    Q   = pinv(R_N) * s_N; 

    er  = norm(A*Q - Yt);
    Gt  = reshape(Q, [R(d-1),1,R(d)]);
    Zt  = Z_true{1,d}(:,t,:);
    fprintf(' er-Yt(%d) = %f \n',t,er)

    G{1,d}(:,t,:) = reshape(Q, [R(d-1), 1, R(d)]); % OK
    Z_tem = G;
    Z_tem{1,d} =  G{1,d}(:,t,:);
%% Calculate and adjust the permanent (non-temporal) core tensors 
  
    for k = 1:d-1
        %disp(k); 
        Y_t_k = modek_unfolding(Y_t,k);
        P_t_k = modek_unfolding(P_t,k);
        
        G_k = ten2mat(G{k},2);
        B   = TCP(Z_tem([k+1:d 1:k-1])); 
        B   = modek_unfolding(B,2);
        ER_Y_t_k =  norm(Y_t_k - G_k*B','fro');
        fprintf(' er-Y_t_k(%d) = %d \n',k,ER_Y_t_k)

        grad = -(P_t_k.*(Y_t_k - G_k*B')*B) + (miu)/t* G_k;
        temp = ten2mat(G{k}, 2) - (theta)* grad;
        T{k} = temp;
    end
%% 
    for k = 1:d-1
         if k == 1
             G{k} = mat2ten(T{k},[R(d), I(1), R(1)],2); 
             % reshape(T{k}, [R(d), I(1), R(1)]);
         else
             G{k} = mat2ten(T{k},[R(k-1), I(k), R(k)],2); 
             %reshape(T{k}, [R(k-1), I(k), R(k)]);
         end
    end
    Z_tem = G;
    Z_tem{1,d} = G{1,d}(:,t,:);
%% Performance evaluation   
    
    % comment: The error should not be evaluated on the entire tensor, but only on the new tensor slice data
    Y_t_re = ten2mat(G{1},2)* modek_unfolding(TCP(Z_tem([2:d])),2)';
    e(t) = norm(ten2mat(Y_t,1) - Y_t_re,'fro')/norm(ten2mat(Y_t,1),'fro') ; 
    fprintf(' reconstructed error = %d\n',e(t))

    % Y_re = ten2mat(G{1},2)* modek_unfolding(TCP(G([2:d])),2)';
    % e(t) = norm(ten2mat(Y,1) - Y_re,'fro')/norm(ten2mat(Y,1),'fro') ; 
     
end
% retrieve Y 
end

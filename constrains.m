function [c,ceq] = constrains(x, sum_I_Mn, v_sk0, C_k, R_sk, Delta, T, K)



cnt = 1;
I_b = x(cnt:cnt+T-1);
cnt = cnt + T;

I_sk = x(cnt:cnt+K*T-1);
cnt = cnt + K*T;

V_sk = x(cnt:cnt+K*T-1);

I_sk = reshape(I_sk,[T,K]);
V_sk = reshape(V_sk,[T,K]);


%% constraint 1
ceq1_t = I_b + sum(I_sk, 2) - sum_I_Mn;


ceq1 = norm(ceq1_t, 1);

%% constraint 2

ceq2_tk = zeros(T, K);

for k=1:K
    sum_val = cumsum( Delta / C_k(k) * I_sk(:,k) + R_sk(k) * abs(I_sk(:,k)) );
    ceq2_tk(:,k) = -V_sk(:,k) - sum_val;
end

ceq2_t = sum(ceq2_tk,2);
ceq2 = norm(ceq2_t, 1);

%% constraint 3
ceq3_k = v_sk0 - V_sk(T, :);
ceq3 = sum(ceq3_k);

% scale ceq3 by T do be within the same decade as the other constrains
ceq3 = ceq3 * T;

%% combine

% ceq1
% ceq2
% ceq3

% TODO: Check for scales of the values maybe one dominates the other
ceq = ceq1 + ceq2 + ceq3;
c  = [];


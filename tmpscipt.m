% first constraint
vals_c1 = I_b + sum(I_sk_out - I_sk_in, 2) - sum(I_Mn, 2);

figure(4)
plot(vals_c1)

% second constraint
vals_c2 = zeros(T,1);
for k=1:K
    sumval = Delta/C_k(k) * I_sk(:,k) + R_sk_max(k) * abs(I_sk(:,k));
    
    vals_c2 = vals_c2 + V_sk(:,k);
    vals_c2 = vals_c2 - cumsum(sumval);
end

figure(5)
plot(vals_c1)

% third consstraint
disp(V_sk(1) - V_sk(end))


% fourth constraint - lb

% (I_sk_out-I_sk_in) <= -L_k
vals_c31 = (I_sk_out-I_sk_in)-L_k;
for k=1:K
    plot(vals_c31(:,k))
    hold on
end
plot(vals_c1 .* 0, 'k')

% (I_sk_out-I_sk_in) <= L_k

vals_c32 = L_k -(I_sk_out-I_sk_in);
for k=1:K
    plot(vals_c32(:,k))
    hold on
end
plot(vals_c1 .* 0, 'k')

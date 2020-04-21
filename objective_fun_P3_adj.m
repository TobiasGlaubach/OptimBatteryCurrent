function fval = objective_fun_P3_adj(x, R_sk, gamma, epsilon, delta, T, K)

I_b = x(1:T);

I_sk = x(T+1:T+K*T);
I_sk = reshape(I_sk,[T,K]);

sum_I_s = sum(R_sk .* norm(I_sk, 1));

fval = gamma * ( epsilon *  mean(I_b) + (1 - epsilon) * var(I_b) ) + delta * sum_I_s;

fval = 100000 * fval;

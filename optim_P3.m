clear;
close all;

debug_lvl = 0
iter_req = 1;
count_failed = 1;

T = 200;
K = 4;

TolCon = 1e-12;
MaxIter = 10000;
MaxFunEvals = 1000000;

%% generate the input data

gen_test_data;

% HACK! this is given nowhere in the paper besides the info 
% "We assume that the voltage and current are measured using discrete 
% signals with a sufficiently small sampling period"

% just set some value which is reasonable for sampling
Delta = 0.01;
 
%% initial conditions

I_b0        = ones(T, 1) * 0;
I_sk_out0   = ones(T, K) * 1;
I_sk_in0    = ones(T, K) * 2;
V_sk0       = ones(T+1, K) * 3;
L_k0        = ones(T, K) * 4;

I_b0        = reshape(I_b0,1,[]).';
I_sk_out0   = reshape(I_sk_out0,1,[]).';
I_sk_in0 	= reshape(I_sk_in0,1,[]).';
V_sk0       = reshape(V_sk0,1,[]).';
L_k0        = reshape(L_k0,1,[]).';

x0 = [  I_b0;
        I_sk_out0;
        I_sk_in0;
        V_sk0;
        L_k0];
    

%% construct the problem

prepare_P3;

%% parameters

% from fig 6 in the paper
epsilon =  0.7;
sigma_1 = 26.0;
sigma_2 =  0.8;
gamma = 0.001;
delta = 1.;


%% nonlinear constraint
nonlcon = [];


%% input for MIAD

% initial values for testrun from fig 7
sigma_1 = 2.0;
sigma_2 = 0.5;
alpha = 2;
beta_1 = 2.0;
beta_2 = 0.5;

if exist('x_start_orig_0.mat','file')
    load('x_start_orig_0.mat', 'x');
    x0 = x;
else
    rng('shuffle')
    x0 = rand(size(x));
end



%% parameters

% for debugging from fig 6 in the paper
epsilon =  0.7;
sigma_1 = 26.0;
sigma_2 =  0.8;
gamma = 0.001;
delta = 1.;


%% solver options

% according to matab documentation the interior-point Algorithm of fmincon
% implements the Barrier Function approach to the minization
% https://de.mathworks.com/help/optim/ug/constrained-nonlinear-optimization-algorithms.html
% options = optimoptions(@fmincon,'Algorithm','interior-point', ...
%             'Display','off', ...
%             'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt} ...
%             );

options = optimoptions(@fmincon,'Algorithm','interior-point', ...
    'PlotFcn',{@optimplotconstrviolation,@optimplotfval,@optimplotfirstorderopt}, ...
    'Display','iter', 'MaxIter', MaxIter, 'MaxFunEvals', MaxFunEvals, ...
    'TolCon', TolCon);

%% MIAD

if alpha < 1
   error('alpha must be bigger than 1') 
end
if beta_1 < 0 && beta_1 >= sigma_1
   error('beta_1 must be bigger than 0 and smaller sigma_1') 
end
if beta_2 < 0 && beta_2 >= sigma_2
   error('beta_2 must be bigger than 0 and smaller sigma_2') 
end


n = 1;
m = 1;
run = 1;
n_iter = 0;
% HACK: this is not in the paper but we need to initialize the values and
% it seems valid, ssince they will be overwritten anyways
temp_1 = sigma_1;
temp_2 = sigma_2;
while n_iter < iter_req
    
    % I_Mn is constant for us
    if debug_lvl > 0
        disp(['sigma_1: ', num2str(sigma_1)])
        disp(['sigma_2: ', num2str(sigma_2)])
        disp(['n_iter: ', num2str(n_iter)])
        
    end
    
    fun = @(x)objective_fun_P3( x, sigma_1, sigma_2, epsilon, gamma, delta, T, F, f);
    
    x0 = x;
    [x, fval, exitflag, output] = ...
        fmincon(fun,x,A,b,Aeq,beq,lb,ub,nonlcon,options);

    output
    
    is_feasible = exitflag ~= -2;
    if(is_feasible)
        disp(['Solution is OK (feasible) exitflag=', num2str(exitflag)])
        disp(['fval: ', num2str(fval)])
    else
        disp('Solution is NOT OK (NOT feasible)')
        output.message
    end
    
    if (is_feasible + count_failed) > 0
        n_iter = n_iter + 1;
    end
    
    if (is_feasible)
        if sigma_1 >= beta_1 && sigma_2 >= beta_2
            if n == 1
                temp_1 = sigma_1;
                sigma_1 = sigma_1 - beta_1;
                n = 2;
            else
                temp_2 = sigma_2;
                sigma_2 = sigma_2 - beta_2;
                n = 1;
            end
            % I_Mn is constant for us
        end
    else
        if temp_1 - sigma_1 == beta_1
            sigma_1 = sigma_1 + beta_1;
        end
        if temp_2 - sigma_2 == beta_2
            sigma_2 = sigma_2 + beta_2;
        end
        
        if temp_1 - sigma_1 ~= beta_1 && temp_2 - sigma_2 ~= beta_2
            if m == 1
                sigma_1 = alpha * sigma_1;
                m = 2;
            else
                sigma_2 = alpha * sigma_2;
                m = 1;
            end
        end
    end
end

%% get results

cnt = 1;
I_b = x(cnt:cnt+T-1);
cnt = cnt + T;

I_sk_out = x(cnt:cnt + K*T-1);
cnt = cnt + K*T;

I_sk_in = x(cnt:cnt + K*T-1);
cnt = cnt + K*T;

V_sk = x(cnt:cnt + K*(T+1)-1);
cnt = cnt + K*(T+1);

L_k = x(cnt:cnt + K*T-1);

I_sk = I_sk_out - I_sk_in;

I_sk = reshape(I_sk,[T,K]);
I_sk_out = reshape(I_sk_out,[T,K]);
I_sk_in = reshape(I_sk_in,[T,K]);
V_sk = reshape(V_sk,[T+1,K]);
L_k = reshape(L_k,[T,K]);
t = 1:T;

%% plot results
figure();

subplot(3,1,1);
plot(t, I_b);
ylabel('I_b')
legend('I_b')

subplot(3,1,2);
l = {};
for k=1:K
    plot(t, I_sk_in(:,k), '-');
    l{k} = ['I_{s', num2str(k), '}'];
    hold on;
end
ylabel('I_{s}')
legend(l)


subplot(3,1,3);
l = {};
for k=1:K
    plot([0, t], V_sk(:,k));
    l{k} = ['V_{s', num2str(k), '}'];
    hold on;
end
ylabel('V_{s}')
legend(l)
xlabel('t')

%% test constrains as given in equations

% first constraint
vals_c1 = I_b + sum(I_sk_out - I_sk_in, 2) - sum(I_Mn, 2);

figure()
plot(vals_c1)
title('I_b + sum(I_{sk}^{out} - I_{sk}^{in}) - sum(I_{M_n}) | should all be zero')

% second constraint
vals_c2 = zeros(T,1);
for k=1:K
    sumval = Delta/C_k(k) * I_sk(:,k) + R_sk_max(k) * abs(I_sk(:,k));
    
    vals_c2 = vals_c2 + V_sk(2:end,k);
    vals_c2 = vals_c2 - cumsum(sumval);
end

figure()
plot(vals_c1)
title('second sconstraint | should all be zero')

% third constraint
disp(['should be zero:', num2str(V_sk(1) - V_sk(end))])

% fourth constraint - lb

% (I_sk_out-I_sk_in) <= -L_k
vals_c31 = (I_sk_out-I_sk_in)-L_k;

figure()
for k=1:K
    plot(vals_c31(:,k))
    hold on
end
plot(vals_c1 .* 0, 'k')
title('( I_{sk}_{out}-I_{sk}^{in} )+L_k  | should all be higher than zero')

% (I_sk_out-I_sk_in) <= L_k
vals_c32 = L_k -(I_sk_out-I_sk_in);

figure()
for k=1:K
    plot(vals_c32(:,k))
    hold on
end
plot(vals_c1 .* 0, 'k')
title('L_k -( I_{sk}^{ou} - I_{sk}^{in} ) | should all be higher than zero')

%% test constraints as implemented in matrix form

figure()
v = Aeq * x - beq;
plot(v < 0)
title('Aeq * x - beq <= 0 | should all be true')

figure()
v = A * x - b;
plot(abs(v) < options.TolCon );
title('A * x - b = 0 | should all be true')

%% test boundaries as implemented in matrix form

figure()
plot(x >= lb)
title('x >= lb | should all be true')

figure()
plot(x <= ub)
title('x <= ub | should all be true')


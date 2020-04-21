clear;
close all;
debug_lvl = 1

v_sk0 = 0;
T = 200;
K = 4;

TolCon = 1e-16;
MaxIter = 10000;
MaxFunEvals = 1000000;


%% generate the inpput data

gen_test_data;


% HACK! this is given nowhere in the paper besides the info 
% "We assume that the voltage and current are measured using discrete 
% signals with a sufficiently small sampling period"

% just set some value which is reasonable for sampling
Delta = 0.01;
 
%% initial conditions

I_b0        = ones(T, 1) * 0;
I_sk0       = ones(T, K) * 1;
V_sk0       = ones(T, K) * 3;

I_b0        = reshape(I_b0,1,[]).';
I_sk0       = reshape(I_sk0,1,[]).';
V_sk0       = reshape(V_sk0,1,[]).';

x0 = [  I_b0;
        I_sk0;
        V_sk0];
    


%% input for MIAD

% initial values for testrun from fig 7
sigma_1 = 2.0;
sigma_2 = 0.5;
alpha = 2;
beta_1 = 2.0;
beta_2 = 0.5;


%% parameters

% for debugging from fig 6 in the paper
epsilon =  0.7;
sigma_1 = 26.0;
sigma_2 =  0.8;
gamma = 0.001;
delta = 1.;


%% set the boundaries

% from mathworks documentation:
% x(i) >= lb(i) for all i.
% x(i) <= ub(i) for all i.

% no constraint on I_b
lb_I_b = ones(T, 1)* -inf;
ub_I_b = ones(T, 1) * inf;

% no constraint on I_sk
lb_I_sk = ones(T*K, 1)* -inf;
ub_I_sk = ones(T*K, 1) * inf;

% 0 <= V_sk <= V_sk_max | for each k
lb_V_sk = zeros(T*K, 1);
ub_V_sk = [];
for k=1:K
    ub_V_sk = [ub_V_sk; ones(T, 1) * V_sk_max(k)];
end

lb = [lb_I_b; lb_I_sk; lb_V_sk];
ub = [ub_I_b; ub_I_sk; ub_V_sk];


%% solver options

% according to MATLAB documentation the interior-point Algorithm of fmincon
% implements the Barrier Function approach to the minization problem
% https://de.mathworks.com/help/optim/ug/constrained-nonlinear-optimization-algorithms.html

options = optimoptions(@fmincon,'Algorithm','interior-point', ...
    'PlotFcn',{@optimplotconstrviolation,@optimplotfval,@optimplotfirstorderopt}, ...
    'Display','off', 'MaxIter', MaxIter, 'MaxFunEvals', MaxFunEvals, ...
    'TolCon', TolCon);

%% nonlinear constraint
nonlcon = @(x)constrains(x, sum(I_Mn,2), v_sk0, C_k, R_sk_max, Delta, T, K);


%% objective function
method = 'logbarrier';

%% initial conditions

% x0 = zeros(size(x0));

% rng('shuffle')
% x0 = rand(size(x0));

load('x_start0.mat', 'x')
x0 = x;

%% MIAD

x = x0;

fun = @(x) objective_fun_P3_adj(x, R_sk_max, gamma, epsilon, delta, T, K);

while 1
    
    x0 = x;
    [x, fval, exitflag, output] = ...
        fmincon(fun,x,[],[],[],[],lb,ub,nonlcon,options);
    
    output
    disp(output.message)
    
    save('x_start0', 'x', 'output');
    
    if exitflag == -1
        break
    end
    
    if exitflag ~= 0
        if exitflag ~= -2
            filename = ['minima_', datestr(now,'yyyy-mm-dd_HHMMSS') '.mat'];
            save(filename)
        end
    
        rng('shuffle');
        x = rand(size(x));
    end
    
end

is_feasible = exitflag ~= -2;
if(is_feasible)
    disp(['Solution is OK (feasible) exitflag=', num2str(exitflag)])
    disp(['fval: ', num2str(fval)])
else
    disp(['Solution is NOT OK (NOT feasible) exitflag=', num2str(exitflag)])
    output.message
end

%% get results

cnt = 1;
I_b = x(cnt:cnt+T-1);
cnt = cnt + T;

I_sk = x(cnt:cnt + K*T-1);
cnt = cnt + K*T;

V_sk = x(cnt:cnt + K*T-1);
cnt = cnt + K*T;

I_sk = reshape(I_sk,[T,K]);
V_sk = reshape(V_sk,[T,K]);
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
    plot(t, I_sk(:,k), '-');
    l{k} = ['I_{s', num2str(k), '}'];
    hold on;
%     plot(t, I_sk(:,k), '--');
%     plot(t, I_sk_out(:,k), ':');
end
ylabel('I_{s}')
legend(l)


subplot(3,1,3);
l = {};
for k=1:K
    plot(t, V_sk(:,k));
    l{k} = ['V_{s', num2str(k), '}'];
    hold on;
end
ylabel('V_{s}')
legend(l)
xlabel('t')

%% test constrains as given in equations

% first constraint
vals_c1 = I_b + sum(I_sk, 2) - sum(I_Mn, 2);

figure()
plot(vals_c1)
title('$$I_b + \sum_{k \in K}(I_{s_k}) - \sum_{n \in N}(I_{M_n}) = 0$$ ??','interpreter','latex')
ylabel('should all be zero')

% second constraint
vals_c2 = zeros(T,1);
for k=1:K
    sumval = cumsum(Delta/C_k(k) * I_sk(:,k) + R_sk_max(k) * abs(I_sk(:,k)));
    vals_c2 = vals_c2 + V_sk(:,k) - sumval;
end

figure()
plot(vals_c2)
title('second sconstraint | should all be zero')

% third constraint
disp(['should be zero:', num2str(v_sk0 - V_sk(end,:))])




clear;
close all;
debug_lvl = 1

%% generate the inpput data
gen_test_data;

% HACK! this is given nowhere in the paper besides the info 
% "We assume that the voltage and current are measured using discrete 
% signals with a sufficiently small sampling period"

% just set some value which is reasonable for sampling
Delta = 200e-6;
 
%% initial conditions

I_b0        = ones(T, 1) * 0;
I_sk_out0   = ones(T, K) * 1;
I_sk_in0    = ones(T, K) * 2;
V_sk0       = ones(T, K) * 3;
L_k0       = ones(T, K) * 4;

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
    


%% constraint 1

% constraint : I_B + sum_k(I_sk_out - I_sk_in) - sum_n(I_Mn) = 0

% eq_lh = I_B + sum(I_sk_out - I_sk_in);
% eq_rh = sum(I_Mn, 2);

% TODO: check math against paper again, not sure if this is right

E = eye(T, T);
A1 = E;
A4 = zeros(T, T*K);
A5 = zeros(T, T*K);

A2 = [];
A3 = [];
for k=1:K
    A2 = [A2,  E];
    A3 = [A3, -E];
end

%        I_b,      I_sk_out, I_sk_in,  V_sk, L_k
eq1_A = [A1,       A2,       -A3,      A4,   A5];
eq1_b = sum(I_Mn, 2);

if debug_lvl > 0
    % test for debugging --> this must not fail
    eq1_A * x0 - eq1_b
end

%% constraint 2

% constraint : (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in = 0

% eq_lh = (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in;
% eq_rh = 0;

% TODO: check math against paper again, not sure if this is right

Z = zeros(T, T);
Z_k = zeros(T, K*T);

% part 1: (A-I) * V_sk

% HACK: The matrice construction in the paper seems to be wrong and will 
% result in wrong dimensions since it will result in K*(T+1) rows, which can 
% not be concat with the rest of the constraints. By my calculations it 
% should be a simple unity matrix for each k. Gonna implement those here.
tmp = eye(T);
A = [];
for t=1:K
    A = [A, tmp];
end

% part 2: D_k_out * I_sk_out - D_k_in * I_sk_in


% HACK: The matrice construction in the paper seems to be wrong and will 
% result in wrong dimensions since it will result in K*(T+1) rows, which can 
% not be concat with the rest of the constraints. By my calculations they 
% should be lower triangular matrices. Gonna implement those here.

%construct lower triangular matrix
tmp = zeros(T);
for t=1:T
    tmp(t,1:t) = 1;
end

% HACK: The equivalence transformations from the paper might be wrong, my
% calculations arrived at R_sk - Delta / C_k for both i_in and i_out
%  --> not 100% sure about this, so I will leave the values from the paper.
D_k_out = [];
D_k_in = [];
for k=1:K
    D_k_out = [D_k_out, R_sk_max(k) + Delta / C_k(k) * tmp];
    D_k_in = [D_k_in, R_sk_max(k) - Delta / C_k(k) * tmp];
end

%         I_b,  I_sk_out, I_sk_in, V_sk, L_k
eq2_A = [ Z,    -D_k_out, -D_k_in, A,    Z_k];

% HACK: there seems to be an error in the paper, since rhs is set to be
% length T in the paper, but needs to be T+1 for the math to work
eq2_b = zeros(T, 1);

if debug_lvl > 0
    % test for debugging --> this must not fail
    D_k_out * I_sk_out0
    D_k_in * I_sk_in0
    A * V_sk0

    eq2_A * x0 - eq2_b
end

%% constraint 3

% constraint : E * V_sk = 0

% eq_lh = E * V_sk
% eq_rh = 0;

Z_k = zeros(1, K*T);

% HACK: there seems to be an error in the paper, since E is set to be
% length T+1 in the paper, but needs to be T for the math to work
E_sub = zeros(1, T);
E_sub(1) = 1;
E_sub(end) = -1;

E = [];
for k=1:K
    E = [E, E_sub];
end

if debug_lvl > 0
    % test for debugging --> this must not fail
    E * V_sk0
end

Z1 = zeros(1, T);
Z2 = zeros(1, K*T);

%         I_b,  I_sk_out,   I_sk_in, V_sk, L_k
eq3_A = [ Z1,   Z2,         Z2,       E,   Z_k ];
eq3_b = 0;

if debug_lvl > 0
    eq3_A * x0 - eq3_b
end

%% combine all constraints to one
Aeq = [
    eq1_A;
    eq2_A;
    eq3_A];

beq = [
    eq1_b;
    eq2_b;
    eq3_b];

if debug_lvl > 0
    % test for debugging --> this must not fail
    Aeq * x0 - beq
end

%% input for MIAD

% initial values for testrun from fig 7
sigma_1 = 2.0;
sigma_2 = 0.5;
alpha = 2;
beta_1 = 2.0;
beta_2 = 0.5;

x0 = x0 .* 0;



%% parameters

% for debugging from fig 6 in the paper
epsilon =  0.7;
sigma_1 = 26.0;
sigma_2 =  0.8;
gamma = 0.001;
delta = 1.;


%% inequality constraint 1

% constraint : -L_k <= I_sk_out - I_sk_in <= L_k
% can be rewritten into: 
%   (1): -L_k <= I_sk_out - I_sk_in
%   (2): -L_k >= I_sk_out - I_sk_in
%  which can be rewritten to:
%   (1): -I_sk_out + I_sk_in - L_k <= 0
%   (2):  I_sk_out - I_sk_in - L_k <= 0
% which in matrix algebra is:

Z          = zeros(T);
Z_k        = zeros(T, K*T);

E_I_sk_out = [];
E_I_sk_in = [];
E_L_k = [];
for k=1:K
    E_I_sk_out = [E_I_sk_out, eye(T)];
    E_I_sk_in  = [E_I_sk_in,  eye(T)];
    E_L_k      = [E_L_k,      eye(T)];
end

%      I_b,  I_sk_out,    I_sk_in,   V_sk,  L_k
A_11 = [Z,   -E_I_sk_out,  E_I_sk_in, Z_k,  -E_L_k];
b_11 = zeros(T, 1);

%      I_b,  I_sk_out,    I_sk_in,   V_sk,  L_k
A_12 = [Z,    E_I_sk_out, -E_I_sk_in, Z_k,  -E_L_k];
b_12 = zeros(T, 1);


%%  combine all inequality constraints
A = [   A_11;
        A_12];

b = [   b_11;
        b_12];

if debug_lvl > 0
    % test for debugging --> this must not fail
   A * x0 - b 
end
%% set the boundaries

% from mathworks documentation:
% x(i) >= lb(i) for all i.
% x(i) <= ub(i) for all i.

% no constraint on I_b
lb_I_b = ones(T, 1)* -inf;
ub_I_b = ones(T, 1) * inf;

% 0 <= I_sk_out
lb_I_sk_out = zeros(T*K, 1);
ub_I_sk_out = ones(T*K, 1) * inf;

% 0 <= I_sk_in
lb_I_sk_in = zeros(T*K, 1);
ub_I_sk_in = ones(T*K, 1) * inf;

% 0 <= V_sk <= V_sk_max | for each k
lb_V_sk = zeros(T*K, 1);
ub_V_sk = [];
for k=1:K
    ub_V_sk = [ub_V_sk; ones(T, 1) * V_sk_max(k)];
end

% no constraint on L_k
lb_L_k = ones(T*K, 1) * -inf;
ub_L_k = ones(T*K, 1) *  inf;

lb = [lb_I_b; lb_I_sk_out; lb_I_sk_in; lb_V_sk; lb_L_k];
ub = [ub_I_b; ub_I_sk_out; ub_I_sk_in; ub_V_sk; ub_L_k];


%% nonlinear constraint
nonlcon = [];

%% solver options

% according to matab documentation the interior-point Algorithm of fmincon
% implements the Barrier Function approach to the minization
% https://de.mathworks.com/help/optim/ug/constrained-nonlinear-optimization-algorithms.html
% options = optimoptions(@fmincon,'Algorithm','interior-point', ...
%             'Display','off', ...
%             'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt} ...
%             );
options = optimoptions(@fmincon,'Algorithm','interior-point', ...
    'Display','off', 'MaxIter', 10000, 'MaxFunEvals', 10000);
iter_req = 1;
count_failed = 0;

%% build objective function for P2

% minimize : sum_k(R_sk * f_sub * L_k)
%   where:
%       R_sk is 1x1
%       L_k is Tx1
%       f_sub is 1xT
%     for each k

%    I_b         I_sk_out,     I_sk_in,      V_sk
f = [zeros(1,T), zeros(1,T*K), zeros(1,T*K), zeros(1,T*K)];
for k=1:K
    f = [f,  R_sk_max(k) * ones(1,T)];
end

if debug_lvl > 0
    % test for debugging --> this must not fail
    f * x0
end


%% objective function
method = 'logbarrier';

fun = @(x)objective_fun_P3(x, sigma_1, sigma_2, epsilon, gamma, delta, T, method, f);

        

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

    [x, fval, exitflag, output] = ...
        fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

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

V_sk = x(cnt:cnt + K*T-1);
cnt = cnt + K*T;

L_k = x(cnt:cnt + K*T-1);

I_sk = I_sk_out - I_sk_in;

I_sk = reshape(I_sk,[T,K]);
I_sk_out = reshape(I_sk_out,[T,K]);
I_sk_in = reshape(I_sk_in,[T,K]);
V_sk = reshape(V_sk,[T,K]);
L_k = reshape(L_k,[T,K]);
t = 1:T;

%% plot results
figure(3);

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

%% test constrains

% first constraint
vals_c1 = I_b + sum(I_sk_out - I_sk_in, 2) - sum(I_Mn, 2);

figure(4)
plot(vals_c1)
title('should all be zero')

% second constraint
vals_c2 = zeros(T,1);
for k=1:K
    sumval = Delta/C_k(k) * I_sk(:,k) + R_sk_max(k) * abs(I_sk(:,k));
    
    vals_c2 = vals_c2 + V_sk(:,k);
    vals_c2 = vals_c2 - cumsum(sumval);
end

figure(5)

plot(vals_c1)

% third constraint
disp(['should be zero:', num2str(V_sk(1) - V_sk(end))])
title('should all be zero')

% fourth constraint - lb

% (I_sk_out-I_sk_in) <= -L_k
vals_c31 = (I_sk_out-I_sk_in)-L_k;

figure(6)

for k=1:K
    plot(vals_c31(:,k))
    hold on
end
plot(vals_c1 .* 0, 'k')
title('should all be higher than zero')

% (I_sk_out-I_sk_in) <= L_k

vals_c32 = L_k -(I_sk_out-I_sk_in);

figure(7)
title('should all be higher than zero')
for k=1:K
    plot(vals_c32(:,k))
    hold on
end
plot(vals_c1 .* 0, 'k')





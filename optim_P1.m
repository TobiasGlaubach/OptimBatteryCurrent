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

I_b0        = reshape(I_b0,1,[]).';
I_sk_out0   = reshape(I_sk_out0,1,[]).';
I_sk_in0 	= reshape(I_sk_in0,1,[]).';
V_sk0       = reshape(V_sk0,1,[]).';

x0 = [  I_b0;
        I_sk_out0;
        I_sk_in0;
        V_sk0];
    
%% for debugging

x = x0;


%% constraint 1

% constraint : I_B + sum_k(I_sk_out - I_sk_in) - sum_n(I_Mn) = 0

% eq_lh = I_B + sum(I_sk_out - I_sk_in);
% eq_rh = sum(I_Mn, 2);

% x = [I_b0, I_sk_out0, I_sk_in0, V_sk0]


% len_I_b0 = length(I_b0);
% len_I_sk_out0 = length(I_sk_out0);
% len_I_sk_in0 = length(I_sk_in0);
% len_V_sk0 = length(V_sk0);

A1 = eye(T, T);
A2 = ones(T, T*K);
A3 = ones(T, T*K);
A4 = zeros(T, T*K);

%        I_b,      I_sk_out, I_sk_in,  V_sk
eq1_A = [A1,       A2,       -A3,      A4];
eq1_b =  sum(I_Mn, 2);

if debug_lvl > 0
    % test for debugging --> this must not fail
    eq1_A * x0 - eq1_b
end

%% constraint 2

% constraint : (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in = 0

% eq_lh = (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in;
% eq_rh = 0;
% Z = zeros(T+1, T);
Z = zeros(T, T);

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

%         I_b,  I_sk_out, I_sk_in, V_sk
eq2_A = [ Z,    -D_k_out, -D_k_in, A ];

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

%         I_b,  I_sk_out,   I_sk_in, V_sk
eq3_A = [ Z1,   Z2,         Z2,       E ];
eq3_b = 0;

if debug_lvl > 0
    eq3_A * x0 - eq3_b
end

%% combine all contraints to one
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

lb = [lb_I_b; lb_I_sk_out; lb_I_sk_in; lb_V_sk];
ub = [ub_I_b; ub_I_sk_out; ub_I_sk_in; ub_V_sk];

%% inequality constraint

% no inequality constraint here, only equality constraint
A = [];
b = [];
% 
% disp('WARNING!!!! SETTING UNEQUALITY TO ZERO! FOR DEBUGGING')
% Aeq = [];
% beq = [];

%% nonlinear constraint
nonlcon = [];

%% solver options

% according to matab documentation the interior-point Algorithm of fmincon
% implements the Barrier Function approach to the minization
% https://de.mathworks.com/help/optim/ug/constrained-nonlinear-optimization-algorithms.html
options = optimoptions(@fmincon,'Algorithm','interior-point', ...
            'Display','off', ...
            'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt} ...
            );
        

%% parameters

% for debugging from fig 6 in the paper
epsilon =  0.7;
sigma_1 = 26.0;
sigma_2 =  0.8;
gamma = 0.001;
delta = 1.;

method = 'logbarrier';

fun = @(x)objective_fun_P1(x, sigma_1, sigma_2, epsilon, T, method);

%% input for MIAD

n = 1;
m = 1;

alpha = rand() + 1 + eps;
beta_1 = rand() * (sigma_1 + eps);
beta_2 = rand() * (sigma_2 + eps);

%% MIAD
run = 1;
while run == 1
    
    % I_Mn is constant for us
    
    [x, fval, exitflag, output] = ...
        fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    output
    
    is_feasible = exitflag == -2;
    if(is_feasible)
        disp('Solution is OK (feasible)')
    else
        disp('Solution is NOT OK (NOT feasible)')
        output.message
    end
    
    
    

    
% HACK: needed to uncomment this, since: 
%       A .) we are only simulating 1 step
%       B .) temp_1 and temp_2 are potentially undefined

%     if is_feasible == 1
%         if sigma_1 >= beta_1 && sigma_2 >= beta_2
%             if n == 1
%                 temp_1 = sigma_1;
%                 sigma_1 = sigma_1 - beta_1;
%                 n = 2;
%             else
%                 temp_2 = sigma_2;
%                 sigma_2 = sigma_2 - beta_2;
%                 n = 1;
%             end
%             % I_Mn is constant for us
%         end
%     else
%         if temp_1 - sigma_1 == beta_1
%             sigma_1 = sigma_1 + beta_1;
%         end
%         if temp_2 - sigma_2 == beta_2
%             sigma_2 = sigma_2 + beta_2;
%         end
%         
%         if temp_1 - sigma_1 ~= beta_1 && temp_p2 - sigma_2 ~= beta_2
%             if m == 1
%                 sigma_1 = alpha * sigma_1;
%                 m = 2;
%             else
%                 sigma_2 = alpha * sigma_2;
%                 m = 1;
%             end
%         end
%     end
    
    run = 0;
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
% cnt = cnt + K*T;

I_sk = I_sk_out - I_sk_in;

I_sk = reshape(I_sk,[T,K]);
I_sk_out = reshape(I_sk_out,[T,K]);
I_sk_in = reshape(I_sk_in,[T,K]);
V_sk = reshape(V_sk,[T,K]);

t = 1:T;

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



%% constraint 1

% constraint : I_B + sum_k(I_sk_out - I_sk_in) - sum_n(I_Mn) = 0

% eq_lh = I_B + sum(I_sk_out - I_sk_in);
% eq_rh = sum(I_Mn, 2);

% in matrix form this should be:
% [I,   I,    I,    I,    I,    -I,   -I,   -I,   -I] * 
% [I_b, I_1o, I_2o, I_3o, I_4o, I_1i, I_2i, I_3i, I_4i]^T

I = eye(T, T);
A1 = I;
A4 = zeros(T, (T+1)*K);
A5 = zeros(T, T*K);

A2 = [];
A3 = [];
for k=1:K
    A2 = [A2,  I];
    A3 = [A3,  I];
end

%        I_b, I_sk_out, I_sk_in,  V_sk, L_k
eq1_A = [A1,  A2,       -A3,      A4,   A5];
eq1_b = sum(I_Mn, 2);

if debug_lvl > 0
    % test for debugging --> this must not fail
    disp('eq1_A * x0 - eq1_b');
    eq1_A * x0 - eq1_b
end

%% constraint 2 original implementation

% constraint : (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in = 0

% eq_lh = (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in;
% eq_rh = 0;

% should in matrix form be 
% [0, A-I, -D_k_out, -D_k_in, 0] *
% [I_b, V_sk, I_sk_out, I_sk_in, L_k]^T
% where the non zero terms must actually be concat k times

% part 1: (A-I) * V_sk
A11 = zeros(1, T);
A11(1,1) = 1;
A21 = eye(T);
A22 = zeros(T,1);

A_k = [ A11, 0;
        A21, A22];

% build block diagonalmatrix with K blocks on diagonal
args = cell(K);
for k=1:K
    args{k} = A_k;
end
A = blkdiag(args{:});
Ia = eye(size(A));

% part 2: D_k_out * I_sk_out - D_k_in * I_sk_in
tmp = [ zeros(1,T); 
        eye(T)];
    
args1 = cell(K);
args2 = cell(K);
for k=1:K
    c1 = R_sk_max(k) + Delta / C_k(k);
    c2 = R_sk_max(k) - Delta / C_k(k);
    args1{k} = c1 * tmp;
    args2{k} = c2 * tmp;
end
D_out = blkdiag(args1{:});
D_in = blkdiag(args2{:});

Z = zeros(size(D_out, 1), T);
Z_k = zeros(size(D_out, 1), K*T);

%         I_b,  I_sk_out, I_sk_in,  V_sk,    L_k
eq2_A = [ Z,    -D_out,   -D_in,    (A-Ia), Z_k];
eq2_b = zeros(size(D_out, 1), 1);


% get all non zero rows
idx = sum(abs(eq2_A), 2) > 0;
% exclude all zero rows
eq2_A = eq2_A(idx,:);
eq2_b = eq2_b(idx,:);


% test for debugging --> this must not fail
if debug_lvl > 0
    disp('D_out * I_sk_out0');
    D_out * I_sk_out0
    disp('D_in * I_sk_in0');
    D_in * I_sk_in0
    disp('A * V_sk0');
    A * V_sk0
    disp('eq2_A * x0 - eq2_b');
    eq2_A * x0 - eq2_b
end


%% constraint 3

% constraint : E * V_sk = 0

% eq_lh = E * V_sk
% eq_rh = 0;

E_sub = zeros(1, T+1);
E_sub(1) = 1;
E_sub(end) = -1;

% stack K times
E = [];
for k=1:K
    E = [E, E_sub];
end

Z1 = zeros(1, T);
Z2 = zeros(1, K*T);

%         I_b,  I_sk_out,   I_sk_in, V_sk, L_k
eq3_A = [ Z1,   Z2,         Z2,       E,   Z2 ];
eq3_b = 0;

if debug_lvl > 0
    % test for debugging --> this must not fail
    disp('E * V_sk0');
    E * V_sk0
    disp('eq3_A * x0 - eq3_b');
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
    disp('Aeq * x0 - beq');
    Aeq * x0 - beq
end


%% inequality constraint 1

% constraint : -L_k <= I_sk_out - I_sk_in <= L_k
% can be rewritten into: 
%   (1): -L_k <= I_sk_out - I_sk_in
%   (2): I_sk_out - I_sk_in <= L_k      | switch sides
% can be rewritten into: 
%   (1): -L_k <= I_sk_out - I_sk_in
%   (2):  L_k >= I_sk_out - I_sk_in     | * -1
% can be rewritten into: 
%   (1): -L_k <=  I_sk_out - I_sk_in    | -I_sk_out | +I_sk_in
%   (2): -L_k <= -I_sk_out + I_sk_in    | +I_sk_out | -I_sk_in
%  which can be rewritten to:
%   (1): -I_sk_out + I_sk_in - L_k <= 0
%   (2):  I_sk_out - I_sk_in - L_k <= 0
% which in matrix algebra is:

Z          = zeros(T*K, T);
Z_k        = zeros(T*K, K*(T+1));

E_I_sk_out = eye(T*K);
E_I_sk_in = eye(T*K);
E_L_k = eye(T*K);

%      I_b,  I_sk_out,    I_sk_in,   V_sk,  L_k
A_11 = [Z,   -E_I_sk_out,  E_I_sk_in, Z_k,  -E_L_k];
b_11 = zeros(T*K, 1);

%      I_b,  I_sk_out,    I_sk_in,   V_sk,  L_k
A_12 = [Z,    E_I_sk_out, -E_I_sk_in, Z_k,  -E_L_k];
b_12 = zeros(T*K, 1);


%%  combine all inequality constraints
A = [   A_11;
        A_12];

b = [   b_11;
        b_12];

if debug_lvl > 0
    % test for debugging --> this must not fail
   disp('A * x0 - b');
   A * x0 - b 
end

%% construct F vector

F = sparse(eye(T-1, T));
for t=1:T-1
    F(t, t+1) = -1;
end

if debug_lvl > 0
    % test for debugging --> this must not fail
   disp('F * I_b0');
   F * I_b0
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
lb_V_sk = zeros((T+1)*K, 1);
ub_V_sk = [];
for k=1:K
    ub_V_sk = [ub_V_sk; ones(T+1, 1) * V_sk_max(k)];
end

% no constraint on L_k
lb_L_k = ones(T*K, 1) * -inf;
ub_L_k = ones(T*K, 1) *  inf;

lb = [lb_I_b; lb_I_sk_out; lb_I_sk_in; lb_V_sk; lb_L_k];
ub = [ub_I_b; ub_I_sk_out; ub_I_sk_in; ub_V_sk; ub_L_k];


%% nonlinear constraint
nonlcon = [];


%% build objective function for P2

% minimize : sum_k(R_sk * f_sub * L_k)
%   where:
%       R_sk is 1x1
%       L_k is Tx1
%       f_sub is 1xT
%     for each k

%    I_b         I_sk_out,     I_sk_in,      V_sk
f = [zeros(1,T), zeros(1,T*K), zeros(1,T*K), zeros(1,(T+1)*K)];
for k=1:K
    f = [f,  R_sk_max(k) * ones(1,T)];
end

if debug_lvl > 0
    % test for debugging --> this must not fail
    disp('f * x0');
    f * x0
end

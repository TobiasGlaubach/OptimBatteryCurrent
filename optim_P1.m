
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


len_I_b0 = length(I_b0);
len_I_sk_out0 = length(I_sk_out0);
len_I_sk_in0 = length(I_sk_in0);
len_V_sk0 = length(V_sk0);

%        I_b,              I_sk_out,               I_sk_in,                  V_sk
eq1_A = [eye(T, len_I_b0), ones(T, len_I_sk_out0), -1*ones(T, len_I_sk_in0), zeros(T, len_V_sk0)];
eq1_b =  sum(I_Mn, 2);


% test for debugging --> this must not fail
eq1_A * x0 - eq1_b

%% constraint 2

% constraint : (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in = 0

% eq_lh = (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in;
% eq_rh = 0;
Z = zeros(T+1, T);

% part 1: (A-I) * V_sk
A11 = zeros(1, T);
A11(1) = 1;
A21 = eye(T);
A22 = zeros(T,1);

A = [   A11, 0;
        A21, A22];
I = eye(size(A));

% part 2: D_k_out * I_sk_out - D_k_in * I_sk_in
tmp = [zeros(1,T); eye(T)];
D_k_out = [];
D_k_in = [];


for k=1:K
    D_k_out = [D_k_out, R_sk_max(k) + Delta / C_k(k) * tmp];
    D_k_in = [D_k_in, R_sk_max(k) - Delta / C_k(k) * tmp];
end

%         I_b,  I_sk_out, I_sk_in, V_sk
eq2_A = [ Z,    -D_k_out, -D_k_in, (A-I) ];
% HACK: there seems to be an error in the paper, since rhs is set to be
% length T in the paper, but needs to be T+1 for the math
eq2_b = zeros(T+1, 1);


% test for debugging --> this must not fail
D_k_out * I_sk_out0
D_k_in * I_sk_in0
(A-I) * V_sk0

eq2_A * x0 - eq2_b

%%

E_sub = zeros(1, T+1);
E_sub(1) = 1;
E_sub(end) = -1;

E = [E, E_sub];
E * V_sk - Z

%%

    I_Mn = M_n(:, n);
    
    
    % todo: optimize here
    x0 = [I_b0, I_sk_out0, I_sk_in0, V_sk0];
    

% options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% x = fmincon(@(x)x,1,[],[],[],[],0,[],[],options)

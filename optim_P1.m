
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

A1 = eye(T, len_I_b0);
A2 = ones(T, len_I_sk_out0);
A3 = -1*ones(T, len_I_sk_in0);
A4 = zeros(T, len_V_sk0);

%        I_b,      I_sk_out, I_sk_in,  V_sk
eq1_A = [A1,       A2,       A3,       A4];
eq1_b =  sum(I_Mn, 2);


% test for debugging --> this must not fail
eq1_A * x0 - eq1_b

%% constraint 2

% constraint : (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in = 0

% eq_lh = (A-I) * V_sk - D_k_out * I_sk_out - D_k_in * I_sk_in;
% eq_rh = 0;
% Z = zeros(T+1, T);
Z = zeros(T, T);

% HACK: The matrice construction in the paper seems to be wrong and will 
% result in wrong dimensions since it will result in K*(T+1) rows, which can 
% not be concat with the rest of the constraints. By my calculations it 
% should be a simple unity matrix. Gonna implement those here.

% part 1: (A-I) * V_sk
% A11 = zeros(1, T);
% A11(1) = 1;
% A21 = eye(T);
% A22 = zeros(T,1);
% 
% A = [   A11, 0;
%         A21, A22];
% I = eye(size(A));

A = eye(size(V_sk0));


% part 2: D_k_out * I_sk_out - D_k_in * I_sk_in


% HACK: The matrice construction in the paper seems to be wrong and will 
% result in wrong dimensions since it will result in K*(T+1) rows, which can 
% not be concat with the rest of the constraints. By my calculations they 
% should be lower triangular matrices. Gonna implement those here.

%consstruct lower triangular matrix
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


% test for debugging --> this must not fail
D_k_out * I_sk_out0
D_k_in * I_sk_in0
A * V_sk0

eq2_A * x0 - eq2_b


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

% test for debugging --> this must not fail
E * V_sk0

Z1 = zeros(1, T);
Z2 = zeros(size(I_sk_out0));

%         I_b,  I_sk_out,   I_sk_in, V_sk
eq3_A = [ Z1,   Z2,         Z,       E ];
eq3_b = 0;

eq3_A * x0 - eq3_b

%% combine all contraints to one
size(eq1_A * x0)
size(eq2_A * x0)
size(eq3_A * x0)

size(eq1_A)
size(eq2_A)
size(eq3_A)


%%

    I_Mn = M_n(:, n);
    
    
    % todo: optimize here
    x0 = [I_b0, I_sk_out0, I_sk_in0, V_sk0];
    

% options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% x = fmincon(@(x)x,1,[],[],[],[],0,[],[],options)

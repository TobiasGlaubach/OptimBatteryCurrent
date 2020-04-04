close all;
clear;
clc;

T = 4
K = 2
N = 2;
debug_lvl = 1;

%%
I_b = sym('I_b', [1, T]);
I_sk_in = [];
I_sk_out = [];
V_sk = [];
L_k = [];
I_Mn = [];

for k=1:K
    I_sk_in  = [I_sk_in,  sym(['I_s', num2str(k), '_in_'], [1, T])];
    I_sk_out = [I_sk_out, sym(['I_s', num2str(k), '_out_'], [1, T])];
    V_sk     = [V_sk,     sym(['V_s', num2str(k), '_0']), sym(['V_s', num2str(k), '_'], [1, T])];
    L_k      = [L_k,      sym(['L_k', num2str(k), '_'], [1, T])];
end

for n=1:N
    I_Mn = [ I_Mn, sym(['I_m', num2str(n), '_'], [1, T]).' ];
end

R_sk_max = sym('R_s', [1, K]);
C_k = sym('C_', [1, K]);
Delta = sym('Delta');
V_sk_max = sym('V_sk_max', [1, K]);

%%

I_b = transpose(I_b);
I_sk_in = I_sk_in.';
I_sk_out = I_sk_out.';
V_sk = V_sk.';
L_k = L_k.';

%%

I_b0 = I_b;
I_sk_out0 = I_sk_out;
I_sk_in0 = I_sk_in;
V_sk0 = V_sk;
L_k0 = L_k;

x0 = [  I_b0;
        I_sk_out0;
        I_sk_in0;
        V_sk0;
        L_k0]
    
%%


prepare_P3;
%% generate the test data as well as M1-M3
K = 4;
N = 3;  % usually 6 but we only have cases 1 to 3

C_k = [50, 150, 310, 350];          % in [F]
V_sk_max = [2.7, 2.7, 2.7, 2.7];    % in [V]
R_sk_max = [20, 14, 2.2, 3.2];    % in [mOhm] --> TODO: check for units

% T = 5;     % for debugging
T = 200;
M1 = 5 * ones(T, 1);
M2 = zeros(T,1);
M3 = ones(T, 1);
t = 0:T-1;

for i=1:10:T-10
    M2(i+8) = 50;
    M2(i+9) = 50;
end

for i=7:35:T
    M3(i) = 100;
end

rng(10000)
% hack, this is not in the paper, I just came up with something 
M4 = rand(T, 1) * 10 + 10 * sin(linspace(0, 150 * 2 * pi, T)');
% M4 = smooth(M4, 4);
M4 = M4 - mean(M4);

M5 = smooth(rand(T, 1) * 20, 10);
M5 = M5 - mean(M5);

M6 = smooth(rand(T, 1) * 40, 10) + 2 * sin(linspace(0, 20 * 2 * pi, T)');
M6 = M6 - mean(M6);
M6(19:40) = 0;
M6(80:90) = 0;
M6(125:135) = 0;
M6(160:175) = 0;

I_Mn = [M1, M2, M3, M4, M5, M6];


%% plot test data M1-M3
figure(1);
subplot(3,1,1);
plot(t, M1);
ylabel('M1')

subplot(3,1,2);
plot(t, M2);
ylabel('M2')
ylim([-5, 60])

subplot(3,1,3);
plot(t, M3);
ylabel('M3')
xlabel('t')
ylim([-5, 110])

%% plot test data M4-M6
figure(1);
subplot(3,1,1);
plot(t, M4);
ylabel('M4')

subplot(3,1,2);
plot(t, M5);
ylabel('M5')

subplot(3,1,3);
plot(t, M6);
ylabel('M6')
xlabel('t')


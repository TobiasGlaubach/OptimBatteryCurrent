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

I_Mn = [M1, M2, M3];


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
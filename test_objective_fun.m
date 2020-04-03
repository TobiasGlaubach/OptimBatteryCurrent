
method = 'logbarrier';

epsilon =  30;
sigma_1 = 25;
sigma_2 =  0.8;
gamma = 0.001;
delta = 1.;

%%
T = 200;

a = rand(T,1) * 20 + 1;
f = zeros(1,T);

figure()
fun = @(x)objective_fun_P3(x, sigma_1, sigma_2, epsilon, gamma, delta, T, method, f);
for i=5:20
    aa = smooth(a, i);
    plot(aa);
    hold on;
    val = fun(aa);
    disp(val)
end

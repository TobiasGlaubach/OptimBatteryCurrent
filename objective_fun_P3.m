function [ fval ] = objective_fun_P3( x, sigma_1, sigma_2, epsilon, gamma, delta, T, F, f)
    
    % get I_b the only needed values from x
    I_b = x(1:T, 1);
    
    val_1 = epsilon * sum(penalty_function(I_b, sigma_1));
    val_2 = (1-epsilon) * sum(penalty_function(F * I_b, sigma_2));
    val_3 = (f * x);
    
    fval = gamma * (val_1 + val_2) + delta * val_3;
end

function theta = penalty_function(u, sigma)
    theta = ones(size(u)) * inf;
    idx = abs(u) < sigma;
    theta(idx) = -sigma^2 * log(1-(u(idx)/sigma).^2);
end
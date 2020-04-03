function [ fval ] = objective_fun_P3( x, sigma_1, sigma_2, epsilon, gamma, delta, T, method, f)
    
    % get I_b the only needed values from x
    I_b = x(1:T, 1);

    
    val_1 = epsilon * sum(penalty_function(I_b, sigma_1, method));
    % HACK: there is a mistake in the paper where the sum starts at 1
    % but if the diff is used, the sum must start at 2
    val_2 = (1-epsilon) * sum(penalty_function(diff(I_b), sigma_2, method));
    
    fval = gamma * (val_1 + val_2) + delta * f * x;
    
    return 
    

end


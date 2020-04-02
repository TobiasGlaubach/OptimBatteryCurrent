function [ fval ] = objective_fun_P1( I_B, sigma_1, sigma_2, epsilon)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    
%     persistent F
%     if isempty(F)
%         T = len(i_B);
%         F = [diag(ones(T,1)), zeros(T,1)] - [zeros(T,1), diag(ones(T,1))];
%         disp(F)
%     end
%     
    
    val_1 = epsilon * sum(penalty_function(I_B, sigma_1));
    % there is a mistake in the paper where the sum starts at 1
    % but if the diff is used, the sum must start at 2
    val_2 = (1-epsilon) * sum(penalty_function(diff(I_b), sigma_2));
    
    fval = val_1 + val_2;
    
    return 
    

end


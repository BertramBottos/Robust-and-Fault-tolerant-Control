function [g,mu_1,idx] = GLR(r,M,mu_0,sigma)
%GLR 
%M =initial guess of window size
%mu_0 = mean value when no fault occured 
%sigma = standard deviation of residual

g = zeros(size(r));     % Test statistics
idx = ones(size(r));    % Index of fault occurence sample estimation
mu_1 = g;               % Estimated parameter
for k = M:length(r)     % For each new sample
    S = zeros(M,1);     % Define the log-likelihood ratio
    z = r(k-M+1:k);     % Define the residual samples inside the window M
    for j = 1:M         % Iterate for all the time instants of the window
        sum_sq = 0;
        for i = j:M
            sum_sq = sum_sq + (z(i) - mu_0);
        end
        S(j) = sum_sq^2/((M - j + 1)*2*sigma^2);
    end
    [g(k), idx(k)] = max(S); % Get the value of g(k) and the sample index
    mu_1(k) = sum(r(k-M+idx(k)+1:k))/(M - idx(k) + 1);
end

end


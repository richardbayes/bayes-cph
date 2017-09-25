function out = rejectsamp(n,C,C1)
    smallval = 1e-6;
    lambda = fzero(@(x) lambda_opt(x,smallval,C,C1),[1e-6,10^6]);
    if (1 - exp(-lambda*smallval)) > .01
        warning('Probability of getting a bad rejection candidate may be greater than 1%');
    end
    M = exp(log(fu(smallval,C,C1,0)) - log(exppdf(smallval,1/lambda))) * 1.05;
    % Rejection sampling
    ii = 1;
    out = zeros(n,1);
    while ii <= n
        cand = exprnd(1/lambda,n,1);
        ind = log(rand(n,1)) < (log(fu_vec(cand,C,C1,0)) - (log(M) + log(exppdf(cand,1/lambda))));
        cand2 = cand(ind);
        ngood = sum(ind);
        topind = min(n, ii + ngood - 1);
        out(ii:topind) = cand2(1:min(ngood,n-ii+1));
        ii = topind + 1;        
    end
end

function out = lambda_opt(lambda,x,C,C1)
    out = exppdf(x,1/lambda) - fu(x,C,C1,0) + .01;
end
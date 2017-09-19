
function out = log_ptheta_t(theta,y,a,a0,b0,aa,bb)
    n = size(y,1);
    ind = (1:n)';
    astar = sum(y(:,2)) + a0;
    B = (a + theta*(n-ind)) ./ (a + theta*(n-ind+1));
    bstar = -a*sum(y(:,1).*log(B)) + b0;
    theinside = gammaln(astar) - gammaln(a0) + a0*log(b0) - astar*log(bstar) + ...
        sum(y(:,2) .* log(-a*log(B)));
    out = -(theinside + log(gampdf(theta,aa,1/bb)));
end
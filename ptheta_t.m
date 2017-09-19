function out = ptheta_t(theta,y,a,a0,b0,aa,bb,C)
    n = size(y,1);
    ntheta = length(theta);
    ind = (1:n)';
    astar = sum(y(:,2)) + a0;
    if(size(theta,1) == 1)
        out = zeros(1,ntheta);
    else
        out = zeros(ntheta,1);
    end
    for ii=1:ntheta
        B = (a + theta(ii)*(n-ind)) ./ (a + theta(ii)*(n-ind+1));
        bstar = -a*sum(y(:,1).*log(B)) + b0;
%         out(ii) = gamma(astar)/gamma(a0) * b0^a0 / bstar^astar * ...
%             exp(sum(y(:,2) .* log(-a*log(B)))) * gampdf(theta(ii),aa,1/bb);
        theinside = gammaln(astar) - gammaln(a0) + a0*log(b0) - astar*log(bstar) + ...
            sum(y(:,2) .* log(-a*log(B))) + C;
        out(ii) = exp(theinside) * gampdf(theta(ii),aa,1/bb);
    end
end


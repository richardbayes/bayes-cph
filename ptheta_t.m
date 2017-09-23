function out = ptheta_t(theta,y,a,omega,aa,bb,C)
    n = size(y,1);
    ntheta = length(theta);
    ind = (1:n)';
    if(size(theta,1) == 1)
        out = zeros(1,ntheta);
    else
        out = zeros(ntheta,1);
    end
    for ii=1:ntheta
        B = (a + theta(ii)*(n-ind)) ./ (a + theta(ii)*(n-ind+1));
        theinside = sum( a*omega*y(:,1) .* log(B) + y(:,2) .* log(-a*omega*log(B)) ) + C;
        out(ii) = exp(theinside) * gampdf(theta(ii),aa,1/bb);
    end
end


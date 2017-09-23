function out = log_ptheta_t(theta,y,a,omega,aa,bb)
    n = size(y,1);
    ind = (1:n)';
    B = (a + theta*(n-ind)) ./ (a + theta*(n-ind+1));
    theinside = sum( a*omega*y(:,1) .* log(B) + y(:,2) .* log(-a*omega*log(B)) );
    out = -(theinside + log(gampdf(theta,aa,1/bb)));
end
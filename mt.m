function [out,C] = mt(y,a,a0,b0,aa,bb)
    out = integral(@(theta)ptheta_t(theta,y,a,a0,b0,aa,bb,0),0,Inf);
    C = 0;
    if out == 0
        options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
        [~,lmax] = fminunc(@(theta) log_ptheta_t(theta,y,a,a0,b0,aa,bb),1,options);
        C = -lmax - log(1000);
        out = integral(@(theta)ptheta_t(theta,y,a,a0,b0,aa,bb,-C),0,Inf);
    end
end


function out = slicesampler(n,C,C1,x,thin)
    out = zeros(n,1);
    for ii=1:(n*thin)
        y = rand*fu(x,C,C1,0);
        xmax = fzero(@(x) fu(x,C,C1,y),x);
        x = rand*xmax;
        if mod(ii,thin) == 0
            out(ii/thin) = x;
        end
    end
end




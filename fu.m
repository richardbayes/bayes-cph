function out = fu(u,C,C1,y)
    out = 1/u * (exp(-C*u) - exp(-C1*u)) / log(C1/C) - y;
end
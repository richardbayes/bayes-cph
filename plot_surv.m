function [GAMPRT,grid,THETAS] = plot_surv(y,X,tree,n,grid,thetas,a,omega,ntheta)
    tree = fatten_tree(tree,X);
    if isempty(omega)
        omega = tree.omega;
    end
    if isempty(a)
        a = tree.a;
    end
    if isempty(ntheta)
        ntheta=1;
    end
    nprt = tree.Ntermnodes;
    % find the MAP estimate for theta in each partition element
    [tind,~] = termnodes(tree);
    prt = zeros(nprt,1);
    for ii=1:nprt
        ind = tree.Allnodes{tind(ii)}.Xind;
        prt(ind) = ii;
    end
    
    
    if isempty(thetas) && (ntheta == 1)
        % find the MAP estimate for theta in each partition element
        options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
        %[tind,~] = termnodes(tree);
        thetas = zeros(nprt,1);
        %prt = zeros(nprt,1);
        for ii=1:nprt
            %ind = tree.Allnodes{tind(ii)}.Xind;
            %prt(ind) = ii;
            %ypart = y(ind,:);
            ypart = y(prt == ii,:);
            [thetatmp,~] = fminunc(@(theta) log_ptheta_t(theta,ypart,a,omega,...
                tree.theta_shape,tree.theta_rate),1,options);
            thetas(ii) = thetatmp;
        end
        THETAS = zeros(ntheta,nprt);
        THETAS(1,:) = thetas;
    elseif isempty(thetas) % draw ntheta draws of theta
        THETAS = zeros(ntheta,nprt);
        for ii=1:nprt
            ypart = y(prt == ii,:);
            thetapart = @(x) ptheta_t(x,ypart,a,omega,...
                tree.theta_shape,tree.theta_rate,0);
            slicesamp = slicesample(1,ntheta,'pdf',thetapart,'burnin',50,'thin',10);
            THETAS(:,ii) = slicesamp;
        end
    end

    if isempty(grid)
        tt = linspace(.01,max(y(:,1)),25);
        ngrid = 25;
    elseif length(grid) == 1
        tt = linspace(.01,max(y(:,1)),grid);
        ngrid = grid;
    else
        tt = grid;
        ngrid = length(grid);
    end
    grid = tt;
    GAM = zeros(n*ntheta,ngrid);
    GAMPRT = zeros(ntheta,n*ntheta,ngrid);
    % For every set of thetas
    for kk=1:ntheta
        theind = (n*(kk-1)+1):(n*kk);
        disp(strcat([num2str(kk),'/',num2str(ntheta)]));
        thetas = THETAS(kk,:);
        A = zeros(size(y,1),1);
        for jj=1:size(y,1)
            tmp = 0;
            for ii=1:length(thetas)
                tmp = tmp + sum(y(prt==ii,1) >= y(jj,1)) * thetas(ii);
            end
            A(jj) = tmp;
        end
        for ii=1:length(tt)
            %disp(strcat([num2str(kk),'/',num2str(nprt),' ',num2str(ii),'/',num2str(length(tt))]));
            thegam = zeros(n,1);
            U = zeros(n,1);
            jj = 1;
            while tt(ii) >= y(jj,1) && jj < size(y,1);
                if jj == 1
                    thediff = y(jj,1);
                else
                    thediff = y(jj,1) - y(jj-1,1);
                end
                thegam = thegam + gamrnd(a*omega*thediff,1/(a + A(jj)),n,1);
                if y(jj,2) % only do for non-censored observations
                    %U = U + slicesampler(nsamp,a + A(jj),a + A(jj+1),.5,10);
                    U = U + rejectsamp(n,a + A(jj),a + A(jj+1));
                    %U = U + slicesample(.5,nsamp,'pdf',@(x) fu(x,a + A(jj),a + A(jj+1),0));
                end
                jj = jj + 1;
            end
            if jj == 1
                thediff = tt(ii);
            else
                thediff = tt(ii) - y(jj-1,1);
            end
            delta = gamrnd(a*omega*thediff,1/(a + A(jj)),n,1);
            %GAM(kk,:,ii) = thegam + U + delta;
            GAM(theind,ii) = thegam + U + delta;
        end
        % Now obtain the GAM function in each partition using the thetas
        tgam = GAM(theind,:);
        for ii=1:nprt
            GAMPRT(ii,theind,:) = tgam*thetas(ii);
        end
        %T(kk,:) = tt;
    end
        
    %end
end
    
        
        
        
% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

sdte_shocks = sigma*(dt^0.5)*e_irf;

for t=1:nval_IRF
    
    if t>1
        
        g0=g1;

        myA=A1{BirfposD(t-1,1),NirfposD(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{BirfposD(t-1,1),NirfposU(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{BirfposU(t-1,1),NirfposD(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{BirfposU(t-1,1),NirfposU(t-1,1)};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;
        
        Nirf(t) = Nirf(t-1,1) + ( alpha*Zeta*((Birf(t-1,1)+Nirf(t-1,1)).^alpha) - delta*(Birf(t-1,1)+Nirf(t-1,1)) - rirf(t-1,1)*Birf(t-1,1) - rhohat*Nirf(t-1,1) )*dt + (Birf(t-1,1)+Nirf(t-1,1))*sdte_shocks(t-1);

    end

    girf(:,:,t)=g1;
    
    Birf(t)=sum(sum(g1.*a2*da));
    Birf(t)=max([Birf(t) Bmin]);
    Birf(t)=min([Birf(t) Bmax]);
    BirfposD(t) =floor((Birf(t)-Bmin)/dB)+1;
    BirfposU(t) = ceil((Birf(t)-Bmin)/dB)+1;
    wB=(B_grid(BirfposU(t))-Birf(t))/dB;

    Nirf(t)=max([Nirf(t) Nmin]);
    Nirf(t)=min([Nirf(t) Nmax]);
    NirfposD(t) =floor((Nirf(t)-Nmin)/dN)+1;
    NirfposU(t) = ceil((Nirf(t)-Nmin)/dN)+1;
    wN=(N_grid(NirfposU(t))-Nirf(t))/dN;

    cirft      = wB*wN*squeeze(c(:,:,BirfposD(t,1),NirfposD(t,1))) + wB*(1-wN)*squeeze(c(:,:,BirfposD(t,1),NirfposU(t,1))) + (1-wB)*wN*squeeze(c(:,:,BirfposU(t,1),NirfposD(t,1))) + (1-wB)*(1-wN)*squeeze(c(:,:,BirfposU(t,1),NirfposU(t,1)));
    cirf(:,:,t)=cirft;
    Cirf(t)    =sum(sum(g1.*cirft*da));
    
    rirf(t,1) = alpha * Zeta * ((Birf(t,1)+Nirf(t,1)).^(alpha-1)) - delta - sigma2*((Birf(t,1)+Nirf(t,1))./Nirf(t,1));
    
end

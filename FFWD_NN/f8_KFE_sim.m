% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This function simulates the dynamics of the distribution given a path of shocks, interpolating over a finer grid that the one employed in the HJB solution

function [KFE_sim] = f8_KFE_sim(KFE_inputs)

    Bsim    = zeros(KFE_inputs.each_sim,1);  % Actual value of B in the simulation 
    BposU   = zeros(KFE_inputs.each_sim,1);  % Closest upper grid point index position (grid used for HJB) 
    BposD   = zeros(KFE_inputs.each_sim,1);  % Closest lower grid point index position (grid used for HJB)

    Nsim    = zeros(KFE_inputs.each_sim,1);
    NposU   = zeros(KFE_inputs.each_sim,1);
    NposD   = zeros(KFE_inputs.each_sim,1);

    BBposU  = zeros(KFE_inputs.each_sim,1);  % Closest upper fine grid point index position (fine grid used for PLM converggence criteria)  
    BBposD  = zeros(KFE_inputs.each_sim,1);
    NNposU  = zeros(KFE_inputs.each_sim,1);
    NNposD  = zeros(KFE_inputs.each_sim,1);

    rsim    = zeros(KFE_inputs.each_sim,1);

    g_big   = zeros(KFE_inputs.nval_a,KFE_inputs.nval_z,KFE_inputs.each_sim);

    for t=1:KFE_inputs.each_sim

        if t==1                        % every simulation starts from the deterministric steady state
            g1       =KFE_inputs.g_ss;
            Nsim(1,1)=KFE_inputs.N_ss;
        else

            g0=g1; 
            % Interpolation step to compute next period's distribution
            % The following block computes the transition g(t+1) = (I - A*dt)^(-1)*g(t) for the four closest grid points

            myA        = KFE_inputs.A1{BposD(t-1,1),NposD(t-1,1)};
            myA        = (speye(KFE_inputs.nval_a*KFE_inputs.nval_z)-myA'*KFE_inputs.dt);
            g1_stacked = myA\g0(:);
            g1_stacked = g1_stacked/sum(g1_stacked*KFE_inputs.da);
            g1DD       = [g1_stacked(1:KFE_inputs.nval_a) g1_stacked(KFE_inputs.nval_a+1:KFE_inputs.nval_a*2)];

            myA        = KFE_inputs.A1{BposD(t-1,1),NposU(t-1,1)};
            myA        = (speye(KFE_inputs.nval_a*KFE_inputs.nval_z)-myA'*KFE_inputs.dt);
            g1_stacked = myA\g0(:);
            g1_stacked = g1_stacked/sum(g1_stacked*KFE_inputs.da);
            g1DU       = [g1_stacked(1:KFE_inputs.nval_a) g1_stacked(KFE_inputs.nval_a+1:KFE_inputs.nval_a*2)];

            myA        = KFE_inputs.A1{BposU(t-1,1),NposD(t-1,1)};
            myA        = (speye(KFE_inputs.nval_a*KFE_inputs.nval_z)-myA'*KFE_inputs.dt);
            g1_stacked = myA\g0(:);
            g1_stacked = g1_stacked/sum(g1_stacked*KFE_inputs.da);
            g1UD       = [g1_stacked(1:KFE_inputs.nval_a) g1_stacked(KFE_inputs.nval_a+1:KFE_inputs.nval_a*2)];

            myA        = KFE_inputs.A1{BposU(t-1,1),NposU(t-1,1)};
            myA        = (speye(KFE_inputs.nval_a*KFE_inputs.nval_z)-myA'*KFE_inputs.dt);
            g1_stacked = myA\g0(:);
            g1_stacked = g1_stacked/sum(g1_stacked*KFE_inputs.da);
            g1UU       = [g1_stacked(1:KFE_inputs.nval_a) g1_stacked(KFE_inputs.nval_a+1:KFE_inputs.nval_a*2)];

            % The we compute distribution next period, averaging those four results
            g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;

            Nsim(t) = Nsim(t-1,1) ...
                + ( KFE_inputs.alpha*KFE_inputs.Zeta*((Bsim(t-1,1)+Nsim(t-1,1)).^KFE_inputs.alpha) ...
                    - KFE_inputs.delta*(Bsim(t-1,1)+Nsim(t-1,1)) ...
                    - rsim(t-1,1)*Bsim(t-1,1) ...
                    - KFE_inputs.rhohat*Nsim(t-1,1) )*KFE_inputs.dt ...
                + (Bsim(t-1,1)+Nsim(t-1,1))*KFE_inputs.sdte_shocks(t);

        end

        g_big(:,:,t)=g1;

        % Before calculating the weights, ensure that the simulation is bounded within the limits

        Bsim(t)=sum(sum(g1.*KFE_inputs.a2*KFE_inputs.da));

        Bsim(t)  =   max([Bsim(t) KFE_inputs.Bmin+0.000001]);
        Bsim(t)  =   min([Bsim(t) KFE_inputs.Bmax-0.000001]);
        BposD(t) = floor((Bsim(t)-KFE_inputs.Bmin)/KFE_inputs.dB)+1;
        BposU(t) =  ceil((Bsim(t)-KFE_inputs.Bmin)/KFE_inputs.dB)+1;
        wB=(KFE_inputs.B_grid(BposU(t))-Bsim(t))/KFE_inputs.dB;       % weight of posD

        Nsim(t)  =   max([Nsim(t) KFE_inputs.Nmin+0.000001]);
        Nsim(t)  =   min([Nsim(t) KFE_inputs.Nmax-0.000001]);
        NposD(t) = floor((Nsim(t)-KFE_inputs.Nmin)/KFE_inputs.dN)+1;
        NposU(t) =  ceil((Nsim(t)-KFE_inputs.Nmin)/KFE_inputs.dN)+1;
        wN=(KFE_inputs.N_grid(NposU(t))-Nsim(t))/KFE_inputs.dN;       % weight of posD

        BBposD(t) = floor((Bsim(t)-KFE_inputs.Bmin)/KFE_inputs.dBB)+1;
        BBposU(t) =  ceil((Bsim(t)-KFE_inputs.Bmin)/KFE_inputs.dBB)+1;
        NNposD(t) = floor((Nsim(t)-KFE_inputs.Nmin)/KFE_inputs.dNN)+1;
        NNposU(t) =  ceil((Nsim(t)-KFE_inputs.Nmin)/KFE_inputs.dNN)+1;

        rsim(t,1) = KFE_inputs.alpha * KFE_inputs.Zeta * ((Bsim(t,1)+Nsim(t,1)).^(KFE_inputs.alpha-1)) - KFE_inputs.delta - KFE_inputs.sigma2*((Bsim(t,1)+Nsim(t,1))./Nsim(t,1));

    end

    KFE_sim.Bsim    = Bsim;
    KFE_sim.BposU   = BposU;
    KFE_sim.BposD   = BposD;

    KFE_sim.Nsim    = Nsim;
    KFE_sim.NposU   = NposU;
    KFE_sim.NposD   = NposD;

    KFE_sim.BBposU  = BBposU;
    KFE_sim.BBposD  = BBposD;
    KFE_sim.NNposU  = NNposU;
    KFE_sim.NNposD  = NNposD;

    KFE_sim.rsim    = rsim;

    KFE_sim.g_big   = g_big;

end

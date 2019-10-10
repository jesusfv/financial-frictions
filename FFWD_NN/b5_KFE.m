% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This code simulates the complete path of the distribution given the
% shock

% INPUTS
KFE_inputs.each_sim = each_sim;
KFE_inputs.nval_a   = nval_a;
KFE_inputs.nval_z   = nval_z;
KFE_inputs.g_ss     = g_ss;
KFE_inputs.N_ss     = N_ss;
KFE_inputs.A1       = A1;
KFE_inputs.a2       = a2;
KFE_inputs.dt       = dt;
KFE_inputs.da       = da;
KFE_inputs.dB       = dB;
KFE_inputs.dN       = dN;
KFE_inputs.dBB      = dBB;
KFE_inputs.dNN      = dNN;
KFE_inputs.Bmin     = Bmin;
KFE_inputs.Bmax     = Bmax;
KFE_inputs.Nmin     = Nmin;
KFE_inputs.Nmax     = Nmax;
KFE_inputs.alpha    = alpha;
KFE_inputs.Zeta     = Zeta;
KFE_inputs.delta    = delta;
KFE_inputs.sigma2   = sigma2;
KFE_inputs.rhohat   = rhohat;
KFE_inputs.B_grid   = B_grid;
KFE_inputs.N_grid   = N_grid;

Bsim    = [];
BposU   = [];
BposD   = [];

Nsim    = [];
NposU   = [];
NposD   = [];

BBposU  = [];
BBposD  = [];
NNposU  = [];
NNposD  = [];

rsim    = [];

g_big   = [];

% MAIN LOOP

b5_use_parallel

if use_parallel == 0

    disp('KFE computations will NOT use the parallel toolbox');

    for itp=1:multi_sim   % Number of threads     

        KFE_inputs.sdte_shocks = sdte_shocks((itp-1)*each_sim+1:itp*each_sim);
        KFE_sim = f8_KFE_sim(KFE_inputs);

        Bsim    = [ Bsim   ; KFE_sim.Bsim   ];
        BposU   = [ BposU  ; KFE_sim.BposU  ];
        BposD   = [ BposD  ; KFE_sim.BposD  ];

        Nsim    = [ Nsim   ; KFE_sim.Nsim   ];
        NposU   = [ NposU  ; KFE_sim.NposU  ];
        NposD   = [ NposD  ; KFE_sim.NposD  ];

        BBposU  = [ BBposU ; KFE_sim.BBposU ];
        BBposD  = [ BBposD ; KFE_sim.BBposD ];
        NNposU  = [ NNposU ; KFE_sim.NNposU ];
        NNposD  = [ NNposD ; KFE_sim.NNposD ];

        rsim    = [ rsim   ; KFE_sim.rsim   ];

        g_big   = cat(3, g_big, KFE_sim.g_big);

    end
    
    clear KFE_sim KFE_inputs; 
    
else
    
    disp('KFE computations will use the parallel toolbox');

    myparpool = gcp();

    % first do calculations using myparpool
    for itp=1:multi_sim   % Number of threads     

        KFE_inputs.sdte_shocks = sdte_shocks((itp-1)*each_sim+1:itp*each_sim);
        f(itp) = parfeval(myparpool,@f8_KFE_sim,1,KFE_inputs);

    end

    % then fetch results
    for itp=1:multi_sim

        KFE_sim = fetchOutputs(f(itp));

        Bsim    = [ Bsim   ; KFE_sim.Bsim   ];
        BposU   = [ BposU  ; KFE_sim.BposU  ];
        BposD   = [ BposD  ; KFE_sim.BposD  ];

        Nsim    = [ Nsim   ; KFE_sim.Nsim   ];
        NposU   = [ NposU  ; KFE_sim.NposU  ];
        NposD   = [ NposD  ; KFE_sim.NposD  ];

        BBposU  = [ BBposU ; KFE_sim.BBposU ];
        BBposD  = [ BBposD ; KFE_sim.BBposD ];
        NNposU  = [ NNposU ; KFE_sim.NNposU ];
        NNposD  = [ NNposD ; KFE_sim.NNposD ];

        rsim    = [ rsim   ; KFE_sim.rsim   ];

        g_big   = cat(3, g_big, KFE_sim.g_big);

    end

    clear KFE_sim KFE_inputs; 
    delete(gcp('nocreate'));
    
end




% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

% some preliminary things I may need later

Bsim_used=[];
Nsim_used=[];
nsim_used=[];
for it_sim=1:multi_sim
    Bsim_used =[Bsim_used ; Bsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)];
    Nsim_used =[Nsim_used ; Nsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)];
    nsim_used =[nsim_used  (it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim];
end
nsim_used=nsim_used';
lvrg_sim= (Bsim_used+Nsim_used)./Nsim_used;

[lvrg_sort,lvrg_index]=sort(lvrg_sim);

sim05_=lvrg_index(round(0.05*size(lvrg_index,1)));
sim05 =nsim_used(sim05_);
g05 = squeeze(g_big(:,:,sim05));
B05 = Bsim(sim05);
N05 = Nsim(sim05);

for it1=2:19
    eval(['sim' num2str(5*it1) '_=lvrg_index(round(' num2str(0.05*it1) '*size(lvrg_index,1)));']);
    eval(['sim' num2str(5*it1) ' =nsim_used(sim' num2str(5*it1) '_);']);
    eval(['g'   num2str(5*it1) ' =squeeze(g_big(:,:,sim' num2str(5*it1) '));']);
    eval(['B'   num2str(5*it1) ' =Bsim(sim' num2str(5*it1) ');']);
    eval(['N'   num2str(5*it1) ' =Nsim(sim' num2str(5*it1) ');']);
end

% now plot everything

c2_PLM_NPLM_phase
c3_sss
c4_errors_zone
c5_IRF

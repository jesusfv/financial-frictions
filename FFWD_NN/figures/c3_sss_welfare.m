% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

V_DD=squeeze(V1(:,:,sss_BposD,sss_NposD));
V_DU=squeeze(V1(:,:,sss_BposD,sss_NposU));
V_UD=squeeze(V1(:,:,sss_BposU,sss_NposD));
V_UU=squeeze(V1(:,:,sss_BposU,sss_NposU));
V1_hlsss = wB*wN*V_DD + wB*(1-wN)*V_DU + (1-wB)*wN*V_UD + (1-wB)*(1-wN)*V_UU;

V_DD=squeeze(V1(:,:,sss_BposD,sss_NposD));
V_DU=squeeze(V1(:,:,sss_BposD,sss_NposU));
V_UD=squeeze(V1(:,:,sss_BposU,sss_NposD));
V_UU=squeeze(V1(:,:,sss_BposU,sss_NposU));
V1_llsss = wB*wN*V_DD + wB*(1-wN)*V_DU + (1-wB)*wN*V_UD + (1-wB)*(1-wN)*V_UU;

% report welfare

ss_BposD =floor((B_ss-Bmin)/dB)+1;
ss_BposU = ceil((B_ss-Bmin)/dB)+1;
wB=(B_grid(ss_BposU)-B_ss)/dB;

ss_NposD =floor((N_ss-Nmin)/dN)+1;
ss_NposU = ceil((N_ss-Nmin)/dN)+1;
wN=(N_grid(ss_NposU)-N_ss)/dN;

V_DD=squeeze(V1(:,:,ss_BposD,ss_NposD));
V_DU=squeeze(V1(:,:,ss_BposD,ss_NposU));
V_UD=squeeze(V1(:,:,ss_BposU,ss_NposD));
V_UU=squeeze(V1(:,:,ss_BposU,ss_NposU));
V1_ss = wB*wN*V_DD + wB*(1-wN)*V_DU + (1-wB)*wN*V_UD + (1-wB)*(1-wN)*V_UU;

disp('Average welfare at each SSS')
disp([ 'Average welfare at DSS:    ' num2str(sum(sum(V1_ss.*g_ss*da))) ])
disp([ 'Average welfare at LL-SSS: ' num2str(sum(sum(V1_llsss.*g_llsss*da))) ])
disp([ 'Average welfare at HL-SSS: ' num2str(sum(sum(V1_hlsss.*g_hlsss*da))) ])
disp(' ')

disp('Consumption equivalent of shift from DSS to each SSS')
disp([ 'LL-SSS: ' num2str( (sum(sum(V1_llsss.*g_llsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 ) ])
disp([ 'HL-SSS: ' num2str( (sum(sum(V1_hlsss.*g_hlsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 ) ])
disp(' ')

disp('Same, as percentage of aggregate consumption at each SSS')
disp([ 'LL-SSS: ' num2str( ((sum(sum(V1_llsss.*g_llsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 )/sum(sum(c_llsss.*g_llsss*da))*100) '%'])
disp([ 'HL-SSS: ' num2str( ((sum(sum(V1_hlsss.*g_hlsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 )/sum(sum(c_hlsss.*g_hlsss*da))*100) '%'])
disp(' ')



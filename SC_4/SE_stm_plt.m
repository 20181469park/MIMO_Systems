
close all
clear
clc




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR_dB=-10;
Mrx=200;
Nrx=50;
T=200;
Ms=1:1:Nrx;
Ntx=Ms;
p_tau_dB_org=5 ;
K_sbs=2;
Mtx=K_sbs;
tau_p=(Ms+Ntx)*K_sbs;

LP='MRC';

f_ora=MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ms,Ms,tau_p,p_tau_dB_org,LP);
f_or1=((T-tau_p)/(T)).*(f_ora);
h=plot(Ms,f_or1,'b-');


hold on
grid on



p_tau_dB=p_tau_dB_org;

[Ms_1,Ms_2]=Opt_stm_search(Mrx,Nrx,K_sbs,SNR_dB,p_tau_dB,T,Mtx,LP);


%calculate final ans

Ms_max=Ms_1;
tau_p_max=(Ms_max+Ms_max)*K_sbs;
fmaxa=MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ms_max,Ms_max,tau_p_max,p_tau_dB,LP);
fmax=((T-tau_p_max)/(T)).*(fmaxa);

hb= plot(Ms_max,fmax,'ks');

hold on




xlim([1 24])


legend('MRC/MRT','Optimal number of streams','Interpreter','latex');



xlabel('Number of streams','Interpreter','latex')
ylabel ('Sum SE (bits/s/Hz)','Interpreter','latex')

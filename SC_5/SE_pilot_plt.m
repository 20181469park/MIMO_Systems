



close all
clear
clc




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR_dB=0;

Mrx=200;
Nrx=50;
T=100;
Ntx=2;
Ms=2;

p_tau_dB_org=-3;
K_sbs=1;
Mtx=K_sbs;

tau_p=(Ms+Ntx)*K_sbs:0.1:25;


    LP='MRC';
    f_ora= MRC_ZF_app_5(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,tau_p,p_tau_dB_org,LP);
    f_or=((T-tau_p)/(T)).*(f_ora);
    plot(tau_p,f_or,'b');

    hold on
    grid on

    p_tau_dB=p_tau_dB_org;
    [x_1,x_2]=Opt_pilot_search(Mrx,Nrx,K_sbs,SNR_dB,p_tau_dB,T,Mtx,Ntx,Ms,LP);
   
    %calculate final ans

    x_max=(x_1+x_2)/2;
    fmaxa= MRC_ZF_app_5(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,x_max,p_tau_dB,LP);
    fmax=((T-x_max)/(T)).*(fmaxa);
    ha=plot(x_max,fmax,'ko');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



legend('MRC analytical','Optimal point','Interpreter','latex')
xlabel('Pilot length $\tau_p$','Interpreter','latex')
ylabel ('Sum SE (bit/s/Hz)','Interpreter','latex')
xlim([4 25])
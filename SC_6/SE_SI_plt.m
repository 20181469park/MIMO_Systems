

close all
clear 
clc

SNR_dB=-5;
Ms=3;
Ntx=3;
T=200;
Mrx=70;
Nrx=70;
K_sbs=3;
tau_p=(Ms+Ntx)*K_sbs;

Mtx=K_sbs;
Mtx2=Mrx;
Ntx2=Nrx;
Nrx2=Ntx;
simul=10^1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LP='MRC';
duplex='FD';
sigma_sk_sq_dB=-20:5:25;
sigma_ck_sq_dB=sigma_sk_sq_dB;
sigma_m_sq_dB=sigma_sk_sq_dB;



[USR1g] =  MRC_ZF_app_6(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,sigma_sk_sq_dB,LP,duplex);
USR1gha=((T-tau_p)/(T)).*(USR1g);
plot(sigma_sk_sq_dB,USR1gha,'b-');
grid on
hold on

[USR1gg] = MRC_ZF_simul_6(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,sigma_sk_sq_dB,sigma_ck_sq_dB,sigma_m_sq_dB,LP,duplex,simul);
USR1gh=((T-tau_p)/(T)).*(USR1gg);
plot(sigma_sk_sq_dB,USR1gh,'bo');




xlabel('SI $$\sigma_l^2$$(dB)', 'interpreter','latex' )
ylabel ('Sum SE (bits/s/Hz)' ,'interpreter','latex')
disp(['Mrx=Nrx=' num2str(Mrx)]);
disp([' Number_small_cells=' num2str(K_sbs)]);
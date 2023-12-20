
close all
clear 
clc

SNR_dB=-20:5:30;
T=200;
p_tau_dB=10;
Mrx=100;
Nrx=50;
Ms=2;
Ntx=2;
K_sbs=4;
Mtx=K_sbs;
tau_p=((Ms+Ntx)*K_sbs);
Nrx2=Ntx;
simul=10^2;




LP='ZF';
[USR1] = MRC_ZF_app_3(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,tau_p,p_tau_dB,LP);
Spec1=((T-tau_p)/(T)).*(USR1);
plot(SNR_dB,Spec1,'r-.');

grid on
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ZF_sim%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[USR11] = MRC_ZF_simul_3(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,tau_p,p_tau_dB,LP,simul);
Spec11=((T-tau_p)/(T)).*(USR11);
plot(SNR_dB,Spec11,'r*');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







legend('ZF analytical','ZF simulation','Interpreter','latex')
xlabel('SNR(dB)','Interpreter','latex')
ylabel ('Sum SE (bits/s/Hz)','Interpreter','latex')

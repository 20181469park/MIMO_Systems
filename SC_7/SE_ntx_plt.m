



close all
clear 
clc

SNR_dB=10;
Ms=2;
Ntx=2;
T=200;
Mrx=10:5:150;
Nrx=10:5:150;
simul=10^1;
K_sbs=4;
tau_p=(Ms+Ntx)*K_sbs;
Mtx=K_sbs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LP='ZF';
duplex='HD';
[USR3a] =  MRC_ZF_app_7(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,LP,duplex);
USR3=((T-tau_p)/(2*T)).*(USR3a);
h1=plot(Mrx,USR3,'r--');
grid on
hold on





xlim([10 150])
legend('ZF analytical','Interpreter','latex')

xlabel('Number of antennas','Interpreter','latex')
ylabel ('Sum SE(bits/s/Hz)','Interpreter','latex')

disp([' Number_small_cells=' num2str(K_sbs)]);
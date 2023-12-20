
close all
clear
clc

SNR_dB=-20:5:30;

simul=10^1;
Ms=2;
Ntx=2;
Mrx=300;
Nrx=150;

K_sbs=2;
Mtx=K_sbs;
LP='ZF';

[USR1] = MRC_ZF_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,LP);

h1=plot(SNR_dB,USR1,'r-.');
grid on
hold on



[USR2] = MRC_ZF_simul(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,simul,LP);
h2=plot(SNR_dB,USR2,'r*');

hold on



legend( 'ZF-D analytical','ZF-D simulation','Interpreter','Latex')

xlim([-20 30])
ylim([0 60])







xlabel('SNR (dB)','Interpreter','Latex')
ylabel ('Sum-rate (bits/s/Hz)','Interpreter','Latex')

%--------------------------------------------------------------------------
function [x_1,x_2]=Opt_pilot_search(Mrx,Nrx,numb_sbs,SNR_dB,p_tau_dB,T,Mtx,Ntx,Ms,LP)



x_low=(Ms+Ntx)*numb_sbs;                            % start of interval
x_high=T;                            % end of interval
R=0.5*(sqrt(5)-1);
d=R*(x_high-x_low);
x_1=x_high-d;
x_2=x_low+d;

f1a= MRC_ZF_app_5(Mrx,Nrx,numb_sbs,SNR_dB,Mtx,Ntx,Ms,x_1,p_tau_dB,LP);
f1=((T-x_1)/(T)).*(f1a);


f2a= MRC_ZF_app_5(Mrx,Nrx,numb_sbs,SNR_dB,Mtx,Ntx,Ms,x_2,p_tau_dB,LP);
f2=((T-x_2)/(T)).*(f2a);

%main loop
tol=1e-4;
err=inf;

while err>tol

    if f1>f2

        x_high=x_2;
        x_2=x_1;
        f2=f1;
        d=R*(x_high-x_low);
        x_1=x_high-d;

        f1a= MRC_ZF_app_5(Mrx,Nrx,numb_sbs,SNR_dB,Mtx,Ntx,Ms,x_1,p_tau_dB,LP);
        f1=((T-x_1)/(T)).*(f1a);

    elseif f1<f2

        x_low=x_1;
        x_1=x_2;
        f1=f2;
        d=R*(x_high-x_low);
        x_2=x_low+d;

        f2a= MRC_ZF_app_5(Mrx,Nrx,numb_sbs,SNR_dB,Mtx,Ntx,Ms,x_2,p_tau_dB,LP);
        f2=((T-x_2)/(T)).*(f2a);

    else
        x_low=(x_1+x_2)/2;
        x_high=x_low;

    end

    err=2*abs(x_high-x_low)/(x_high-x_low);

end

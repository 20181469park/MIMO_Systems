function [Ms_1,Ms_2]=Opt_stm_search(Mrx,Nrx,K_sbs,SNR_dB,p_tau_dB,T,Mtx,LP)



Ms_low=1;                               % start of interval
Ms_high=Nrx;                            % end of interval


R=0.5*(sqrt(5)-1);
d=R*(Ms_high-Ms_low);
Ms_1=round(Ms_high-d);
tau_p_1=(Ms_1+Ms_1)*K_sbs;
Ms_2=round(Ms_low+d);
tau_p_2=(Ms_2+Ms_2)*K_sbs;

f1a=MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ms_1,Ms_1,tau_p_1,p_tau_dB,LP);
f1=((T-tau_p_1)/(T)).*(f1a);


f2a=MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ms_2,Ms_2,tau_p_2,p_tau_dB,LP);
f2=((T-tau_p_2)/(T)).*(f2a);



tol=1e-4;
err=inf;

while err>tol
    
    if f1>f2
        
        Ms_high=Ms_2;
        Ms_2=Ms_1;
        f2=f1;
        d=R*(Ms_high-Ms_low);
        Ms_1=round(Ms_high-d);
        tau_p_1=(Ms_1+Ms_1)*K_sbs;

        f1a=MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ms_1,Ms_1,tau_p_1,p_tau_dB,LP);
        f1=((T-tau_p_1)/(T)).*(f1a);

    elseif f1<f2
       
        Ms_low=Ms_1;
        Ms_1=Ms_2;
        f1=f2;
        d=R*(Ms_high-Ms_low);
        Ms_2=round(Ms_low+d);
        tau_p_2=(Ms_2+Ms_2)*K_sbs;

        f2a=MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ms_2,Ms_2,tau_p_2,p_tau_dB,LP);
        f2=((T-tau_p_2)/(T)).*(f2a);

    else
        Ms_low=(Ms_1+Ms_2)/2;
        Ms_high=Ms_low;
   end
   
    err=2*abs(Ms_high-Ms_low)/(Ms_high-Ms_low);
end

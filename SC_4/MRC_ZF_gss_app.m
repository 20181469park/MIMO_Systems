function [USR] = MRC_ZF_gss_app(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,tau_p,p_tau_dB,LP)





sigma_sk_sq=0.3;
sigma_ck_sq=0.3;
sigma_m_sq=0.3;


p_tau=10^(p_tau_dB/10);

USR=zeros(1,length(Ms));

for ind_m=1:length(Ms)
    N_one=Ntx(ind_m);
    Ms_one=Ms(ind_m);
    beta_ulk=0.1*ones(N_one,K_sbs);
    beta_dlk=0.1*ones(Ms_one,K_sbs);
    SNR=10^(SNR_dB/10);
    p1=SNR;
    p2= SNR;
    hat_beta_ulk_sq=(tau_p(ind_m)*p_tau.*beta_ulk.^2)./(1+(tau_p(ind_m)*p_tau.*beta_ulk));
    hat_beta_dlk_sq=(tau_p(ind_m)*p_tau.*beta_dlk.^2)./(1+(tau_p(ind_m)*p_tau.*beta_dlk));


    Num1=zeros(N_one,1);
    sum_Its=zeros(N_one,1);
    It=zeros(N_one,1);
    Denum1=zeros(N_one,1);
    R_mka=zeros(N_one,1);

    Num2=zeros(Ms_one,1);
    sum_Int_dl=zeros(Ms_one,1);
    Ls_sc_N_dl=zeros(Ms_one,1);
    Denum2=zeros(Ms_one,1);
    R_ska=zeros(Ms_one,1);


    for kk=1:K_sbs

        if strcmp('MRC',LP)

            for nn=1:N_one
                Num1(nn,kk)=p1*hat_beta_ulk_sq(nn,kk)*Mrx;


                It(nn,kk)=p1*sum(beta_ulk(:));


                Denum1(nn,kk)=(p2*sigma_m_sq*Mtx*Ms_one)+1;


                R_mka(nn,kk)=log2(1+(Num1(nn,kk)/(Denum1(nn,kk)+It(nn,kk))));


            end
        else
            for nn=1:N_one
                til_beta_ulk=beta_ulk-hat_beta_ulk_sq;


                Num1(nn,kk)=p1*(Mrx-(K_sbs*N_one))*hat_beta_ulk_sq(nn,kk);

                sum_Its(nn,kk)=p1*sum(til_beta_ulk(:));

                Denum1(nn,kk)=(p2*sigma_m_sq*Mtx*Ms_one)+1;

                R_mka(nn,kk)=log2(1+(Num1(nn,kk)/(Denum1(nn,kk)+sum_Its(nn,kk))));

            end
        end
        R_mk=sum(R_mka(:,kk));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp('MRC',LP)

            for ms=1:Ms_one
                Num2(ms,kk)=p2*hat_beta_dlk_sq(ms,kk)*Nrx;

                sum_Int_dl(ms,kk)= p1*sum(beta_dlk(:));

                Denum2(ms,kk)=(p1*(sigma_sk_sq)*N_one)+(p1*((K_sbs -1)*sigma_ck_sq*N_one))+1;

                R_ska(ms,kk)=log2(1+(Num2(ms,kk)/(Denum2(ms,kk)+sum_Int_dl(ms,kk))));

            end
        else
            for ms=1:Ms_one

                Num2(ms,kk)=p2*(Nrx-(K_sbs*Ms_one))*hat_beta_dlk_sq(ms,kk);

                til_beta_dlk=beta_dlk-hat_beta_dlk_sq;

                sum_Int_dl(ms,kk)=p2*sum(til_beta_dlk(:));

                Ls_sc_N_dl(ms,kk)=(p1*(sigma_sk_sq)*N_one)+(p1*((K_sbs-1)*sigma_ck_sq)*N_one)+1;

                R_ska(ms,kk)=log2(1+(Num2(ms,kk)/(sum_Int_dl(ms,kk)+Ls_sc_N_dl(ms,kk))));

            end
        end

        R_sk=sum(R_ska(:,kk));


        rate_sum=(R_mk+R_sk);

        USR(1,ind_m)=USR(1,ind_m)+rate_sum;
    end


end





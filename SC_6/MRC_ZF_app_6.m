function [USR] = MRC_ZF_app_6(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,sigma_sk_sq_dB,LP,duplex)




tau_p=(Ms+Ntx)*K_sbs;
beta_ulk=1*ones(Ntx,K_sbs);
beta_dlk=1*ones(Ms,K_sbs);



USR=zeros(1,length(sigma_sk_sq_dB));

for ind_lp=1:length(sigma_sk_sq_dB)

    if strcmp('HD',duplex)
       
        p1=2*10^(SNR_dB/10);

        sigma_sk_sq=0;                                   % SBS self interference
        sigma_ck_sq=sigma_sk_sq;                                   % jth SBS- kth SBS interference
        sigma_m_sq=sigma_sk_sq;                                    % MBS self interference


       
    else
        p1=10^(SNR_dB/10);
        sigma_sk_sq=10.^(sigma_sk_sq_dB(ind_lp)/10);      % SBS self interference
        sigma_ck_sq= sigma_sk_sq;      % jth SBS- kth SBS interference
        sigma_m_sq= sigma_sk_sq;        % MBS self interference

    end
    p2=p1;
    p_tau=p1;

    hat_beta_ulk_sq=(tau_p*p_tau.*beta_ulk.^2)./(1+(tau_p*p_tau.*beta_ulk));
    hat_beta_dlk_sq=(tau_p*p_tau.*beta_dlk.^2)./(1+(tau_p*p_tau.*beta_dlk));

    Num1=zeros(Ntx,1);
    It=zeros(Ntx,1);
    Denum1=zeros(Ntx,1);
    sum_Its=zeros(Ntx,1);
    R_mka=zeros(Ntx,1);

    Num2=zeros(Ms,1);
    sum_Int_dl=zeros(Ms,1);
    Denum2=zeros(Ms,1);
    Ls_sc_N_dl=zeros(Ms,1);
    R_ska=zeros(Ms,1);

    for kk=1:K_sbs

       if strcmp('MRC',LP)

                for nn=1:Ntx

                    Num1(nn,kk)=p1*hat_beta_ulk_sq(nn,kk)*Mrx;


                    It(nn,kk)=p1*sum(beta_ulk(:));


                    Denum1(nn,kk)=(p2*sigma_m_sq*Mtx*Ms)+1;


                    R_mka(nn,kk)=log2(1+(Num1(nn,kk)/(Denum1(nn,kk)+It(nn,kk))));



                end
            

       else

                for nn=1:Ntx

                    til_beta_ulk=beta_ulk-hat_beta_ulk_sq;


                    Num1(nn,kk)=p1*(Mrx-(K_sbs*Ntx))*hat_beta_ulk_sq(nn,kk);

                    sum_Its(nn,kk)=p1*sum(til_beta_ulk(:));

                    Denum1(nn,kk)=(p2*sigma_m_sq*Mtx*Ms)+1;

                    R_mka(nn,kk)=log2(1+(Num1(nn,kk)/(Denum1(nn,kk)+sum_Its(nn,kk))));


                end
        end

        R_mk=sum(R_mka(:,kk));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp('MRC',LP)
          
            for ms=1:Ms

                Num2(ms,kk)=p2*hat_beta_dlk_sq(ms,kk)*Nrx;

                sum_Int_dl(ms,kk)= p1*sum(beta_dlk(:));

                Denum2(ms,kk)=(p1*(sigma_sk_sq)*Ntx)+(p1*((K_sbs -1)*sigma_ck_sq*Ntx))+1;

                R_ska(ms,kk)=log2(1+(Num2(ms,kk)/(Denum2(ms,kk)+sum_Int_dl(ms,kk))));

            end


        else

            for ms=1:Ms


                Num2(ms,kk)=p2*(Nrx-(K_sbs*Ms))*hat_beta_dlk_sq(ms,kk);

                til_beta_dlk=beta_dlk-hat_beta_dlk_sq;

                sum_Int_dl(ms,kk)=p2*sum(til_beta_dlk(:));

                Ls_sc_N_dl(ms,kk)=(p1*(sigma_sk_sq)*Ntx)+(p1*((K_sbs -1)*sigma_ck_sq)*Ntx)+1;

                R_ska(ms,kk)=log2(1+(Num2(ms,kk)/(sum_Int_dl(ms,kk)+Ls_sc_N_dl(ms,kk))));

            end
        end
        R_sk=sum(R_ska(:,kk));
        rate_sum=R_mk+R_sk;



        USR(1,ind_lp)=USR(1,ind_lp)+rate_sum;
    end


end
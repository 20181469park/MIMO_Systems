

function [USR] =MRC_ZF_simul(Mrx,Nrx,K_sbs,SNR_dB,Mtx,Ntx,Ms,simul,LP)




tau_p=(Ms+Ntx)*K_sbs;

beta_ulk=0.1*ones(Ntx,K_sbs);
beta_dlk=0.1*ones(Ms,K_sbs);

sigma_sk_sq=0.3;                % SBS self interference
sigma_ck_sq=0.3;      % jth SBS- kth SBS interference
sigma_m_sq=0.3;                 % MBS self interference

p_tau_dB=10;
p_tau=10^(p_tau_dB/10);
USR=zeros(1,length(SNR_dB));



for ind=1:length(SNR_dB)

    SNR=10^(SNR_dB(ind)/10);
    p1=SNR;
    p2= SNR;
    hat_beta_ulk_sq=(tau_p*p_tau.*beta_ulk.^2)./(1+(tau_p*p_tau.*beta_ulk));
    hat_beta_dlk_sq=(tau_p*p_tau.*beta_dlk.^2)./(1+(tau_p*p_tau.*beta_dlk));
    til_beta_ulk_sq=(beta_ulk)./(1+(tau_p*p_tau.*beta_ulk));
    til_beta_dlk_sq=(beta_dlk)./(1+(tau_p*p_tau.*beta_dlk));


    hat_Beta_ulk_sq = diag((hat_beta_ulk_sq(:)));
    til_Beta_ulk_sq = diag((til_beta_ulk_sq(:)));
    hat_Beta_dlk_sq = diag((hat_beta_dlk_sq(:)));
    til_Beta_dlk_sq = diag((til_beta_dlk_sq(:)));

    Av_ul = zeros(Ntx,K_sbs);
    Aav_ul = zeros(Ntx,K_sbs);
    Bv_ul = zeros(Ntx,K_sbs);
    Cv_ul = zeros(Ntx,K_sbs);
    Dv_ul = zeros(Ntx,K_sbs);
    Ev_ul = zeros(Ntx,K_sbs);


    Av_dl = zeros(Ms,K_sbs);
    Aav_dl = zeros(Ms,K_sbs);
    Bv_dl = zeros(Ms,K_sbs);
    Cv_dl = zeros(Ms,K_sbs);
    Dv_dl = zeros(Ms,K_sbs);
    Ev_dl = zeros(Ms,K_sbs);
    Fv_dl = zeros(Ms,K_sbs);

    for  expn=1:simul


        hat_H_ul =sqrt(1/2)*(randn(Mrx,K_sbs*Ntx)+1i*randn(Mrx,K_sbs*Ntx))*(sqrt(hat_Beta_ulk_sq));
        til_H_ul=sqrt(1/2)*(randn(Mrx,K_sbs*Ntx)+1i*randn(Mrx,K_sbs*Ntx))*(sqrt(til_Beta_ulk_sq));
        H_ul=hat_H_ul+til_H_ul;


        hat_H_dl =sqrt(1/2)*(randn(K_sbs*Nrx,Mtx*Ms)+1i*randn(K_sbs*Nrx,Mtx*Ms))*(sqrt(hat_Beta_dlk_sq));
        til_H_dl=sqrt(1/2)*(randn(K_sbs*Nrx,Mtx*Ms)+1i*randn(K_sbs*Nrx,Mtx*Ms))*(sqrt(til_Beta_dlk_sq));
        H_dl=hat_H_dl+til_H_dl;


        G_m =sqrt(1/2)*(randn(Mrx,Mtx*Ms)+1i*randn(Mrx,Mtx*Ms))*(sqrt(sigma_m_sq));

        g_sk =sqrt(1/2)*(randn(Nrx,K_sbs*Ntx)+1i*randn(Nrx,K_sbs*Ntx))*(sqrt(sigma_sk_sq));
        g_ckj_min_k =sqrt(1/2)*(randn(Nrx,(K_sbs-1)*Ntx)+1i*randn(Nrx,(K_sbs-1)*Ntx)).*(sqrt(sigma_ck_sq));







        Ak_ul = zeros(Ntx,K_sbs);
        Aak_ul = zeros(Ntx,K_sbs);
        bk_ul = zeros(Ntx,K_sbs);
        Ck_ul = zeros(Ntx,K_sbs);
        Dk_ul = zeros(Ntx,K_sbs);
        Ek_ul = zeros(Ntx,K_sbs);


        Ak_dl = zeros(Ms,K_sbs);
        Aak_dl = zeros(Ms,K_sbs);
        bk_dl = zeros(Ms,K_sbs);
        Ck_dl = zeros(Ms,K_sbs);
        Dk_dl = zeros(Ms,K_sbs);
        Ek_dl = zeros(Ms,K_sbs);
        Fk_dl = zeros(Ms,K_sbs);

        for k_ind = 1:K_sbs

            for nn=1:Ntx
                H_ulk=H_ul(:,(k_ind-1)*Ntx+1:(k_ind*Ntx));
                h_ulkn=H_ulk(:,nn);
                h_ulkn_bar=H_ulk(:,[1:nn-1 nn+1:end]);

                hat_H_ulk=hat_H_ul(:,(k_ind-1)*Ntx+1:(k_ind*Ntx));
                hat_h_ulkn=hat_H_ulk(:,nn);

                if strcmp('MRC',LP)

                        wul_kn=hat_h_ulkn;

                else
                        Wzf_ul=hat_H_ul*inv(  hat_H_ul'*  hat_H_ul);
                        wzf_ul_k=Wzf_ul(:,(k_ind-1)*Ntx+1:(k_ind*Ntx));
                        wul_kn=wzf_ul_k(:,nn);
                end



                H_ul_i=H_ul(:,[1:(k_ind-1)*Ntx k_ind*Ntx+1:end]);





           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                Ak_ul(nn,k_ind) = p1*abs((wul_kn'*(h_ulkn*h_ulkn')*(wul_kn)));


                Aak_ul(nn,k_ind) = p1*abs((wul_kn'*(h_ulkn_bar*h_ulkn_bar')*(wul_kn)));

                bk_ul(nn,k_ind)= sqrt(p1)*abs(wul_kn'*h_ulkn);


                Ck_ul(nn,k_ind) = p1*abs((wul_kn'*(H_ul_i*H_ul_i')*(wul_kn)));


                Dk_ul(nn,k_ind) = (p2)*abs(wul_kn'*(G_m*G_m')*(wul_kn));

                Ek_ul(nn,k_ind) = abs(wul_kn'*(wul_kn));
            end
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            hat_H_dl_k=hat_H_dl((k_ind-1)*Nrx+1:(k_ind*Nrx),:);
            H_dl_k=H_dl((k_ind-1)*Nrx+1:(k_ind*Nrx),:);
            H_dl_kmt=H_dl_k(:,(k_ind-1)*Ms+1:(k_ind*Ms));

            for ms=1:Ms

                if strcmp('MRC',LP)
                       
                        Wmrc_dl_k=hat_H_dl_k;
                        wmrc_dl_kms=Wmrc_dl_k(:,(k_ind-1)*Ms+1:(k_ind*Ms));
                        w_dl_km=wmrc_dl_kms(:,ms);

                else
                        Wzf_dl_k=hat_H_dl_k*inv(hat_H_dl_k'* hat_H_dl_k);
                        wzf_dl_kms=Wzf_dl_k(:,(k_ind-1)*Ms+1:(k_ind*Ms));
                        w_dl_km=wzf_dl_kms(:,ms);

                end

                h_dlkm=H_dl_kmt(:,ms);
                h_dlk_mbar=H_dl_kmt(:,[1:ms-1 ms+1:end]);


                H_dl_ki=H_dl_k(:,[1:(k_ind-1)*Ms k_ind*Ms+1:end]);

                g_sk_n=g_sk(:,(k_ind-1)*Ntx+1:(k_ind*Ntx));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                Ak_dl(ms,k_ind) = p2*abs((w_dl_km'*(h_dlkm*h_dlkm')*(w_dl_km)));

                Aak_dl(ms,k_ind) =  p2*abs((w_dl_km'*(h_dlk_mbar*h_dlk_mbar')*(w_dl_km)));

                bk_dl(ms,k_ind)= sqrt(p2)*abs(w_dl_km'*h_dlkm);

                Ck_dl(ms,k_ind) = p2*abs((w_dl_km'*(H_dl_ki*H_dl_ki')*(w_dl_km)));

                Dk_dl(ms,k_ind) = (p1)*abs(w_dl_km'*(g_sk_n *g_sk_n' )*(w_dl_km));

                Ek_dl(ms,k_ind) = abs(w_dl_km'*(w_dl_km));

                Fk_dl(ms,k_ind) = p1*abs((w_dl_km'*(g_ckj_min_k*g_ckj_min_k')*(w_dl_km)));


            end
        end
        Av_ul=Av_ul+Ak_ul;
        Aav_ul=Aav_ul+Aak_ul;
        Bv_ul=Bv_ul+bk_ul;
        Cv_ul=Cv_ul+Ck_ul;
        Dv_ul=Dv_ul+Dk_ul;
        Ev_ul=Ev_ul+Ek_ul;

        Av_dl=Av_dl+Ak_dl;
        Aav_dl=Aav_dl+Aak_dl;
        Bv_dl=Bv_dl+bk_dl;
        Cv_dl=Cv_dl+Ck_dl;
        Dv_dl=Dv_dl+Dk_dl;
        Ev_dl=Ev_dl+Ek_dl;
        Fv_dl=Fv_dl+Fk_dl;

    end


    A_ul=Av_ul/simul;
    Aa_ul=Aav_ul/simul;
    b_ul=Bv_ul/simul;
    C_ul=Cv_ul/simul;
    D_ul = abs(Dv_ul/simul);
    E_ul = abs(Ev_ul/simul);


    A_dl=Av_dl/simul;
    Aa_dl=Aav_dl/simul;
    b_dl=Bv_dl/simul;
    C_dl=Cv_dl/simul;
    D_dl= abs(Dv_dl/simul);
    E_dl = abs(Ev_dl/simul);
    F_dl = abs(Fv_dl/simul);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vars=A_ul-b_ul.^2;
    Denums= vars+Aa_ul+C_ul+D_ul+E_ul;

    Nums=b_ul.^2;
    USRa1vs = log2(1+(Nums./((Denums))));


    var2s=A_dl-b_dl.^2;
    Denum2s= var2s+Aa_dl+C_dl+D_dl+E_dl+F_dl;

    Num2s=b_dl.^2;
    USRa2vs = log2(1+(Num2s./((Denum2s))));



   rate_sum=USRa1vs+USRa2vs;
    sumt=sum(sum(rate_sum));

    USR(ind)=abs(sumt);
    
end
























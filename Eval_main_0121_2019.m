%%
% Sample codes for the VR-APSA algorithm
%
%Publication
%IEICE TRANSACTIONS on Fundamentals of Electronics, Communications and Computer Sciences   Vol.E102-A    No.5    pp.725-728
%Publication Date: 2019/05/01
%Online ISSN: 1745-1337
%DOI: 10.1587/transfun.E102.A.725
%Type of Manuscript: LETTER
%Category: Digital Signal Processing

% Please cite the following work:
% Y.-R. Chien, "Variable regularization affine projection sign algorithm in impulsive
% noisy environment", IEICE Trans. Fundam. Electron. Commun. Comput. Sci., vol. 102,
% no. 5, pp. 725-728, May 2019.
%%
clear;
close all;
clc;
%Definitions:
% number of realizations within the ensemble
ensemble= 10;
% Number of iterations for each ensemble
Tran_len= 8e4*2;

% Four unknown system models can be used by setting SYS_MODE
% SYS_MODE=1 ==> load from h0.mat with length 30
% SYS_MODE=2 ==> two fixed unknown models, each with length 16
% SYS_MODE=3 ==> load from echo2.mat with length 512
% SYS_MODE=4 ==> random generated unknown models with length 128

SYS_MODE=4;

% BG IN is assumed.
% We further assumed that the impulsive noise occurs in some certain
% regions (please refer to the paper)
% IN_ON =1 ==> IN appeared
% IN_ON =0 ==> No IN

IN_ON=1;

% Two types of inputs can be used.
% INPUT_MODE=1 ==> AR(0.9)
% INPUT_MODE=2 ==> ARMA(2,2)
INPUT_MODE=1;

switch SYS_MODE
    case 1
        load('ho.mat');
        H1=wo;
        H1=H1./sum(H1);
        H2=-H1;
        Flt_Len=length(H1);% number of coefficients of the adaptive filter
        mu= 0.2;% 0.3 convergence factor (step)  (0 < mu < 1)
        mu_APSA=mu;

    case 2

        H1= [0.0003,0.0007,-0.0005,-0.0009,0.0019,0.0019,-0.0028,-0.0032,0.0051,0.0071,-0.0077,-0.0165,0.0086,0.0309,-0.0032,-0.0464].';
        H2=[0.0002,0.0008,-0.0006,0.0015,-0.0022,-0.0010,-0.0031,-0.0015,-0.0034,-0.0026,-0.0024,-0.0034,-0.0008,-0.0038,0.0011,-0.0037].';
        Flt_Len=length(H1);
        mu= 0.005;% convergence factor (step)  (0 < mu < 1)
        mu_APSA=mu;

    case 3 %

        load('echo2.mat');
        H1=echo2;
        H1=H1./(norm(H1));
        H2=-H1;
        Flt_Len=length(H1);
        %         H1=randn(Flt_Len,1);H1=H1./sum(H1);
        %         H2=randn(Flt_Len,1);H2=H2./sum(H2);
        mu_APSA=0.2;
        mu= 0.05;% convergence factor (step)  (0 < mu < 1)

    case 4 % Random

        Flt_Len=128;
        H1=randn(Flt_Len,1);H1=H1./sum(H1);
        %  H2=randn(Flt_Len,1);H2=H2./sum(H2);
        H2=-H1;
        mu_APSA=0.1;
        mu= 0.2;% 0.15 convergence factor (step)  (0 < mu < 1)

    otherwise
        disp('SYS_MODE error\n');

end


gamma3=1e-1; % for curve-(iii), delta=0.1
gamma2=1e5;  % for curve-(iv), delta=1e5
rho=9;   %
P=3;% The real projection order is P+1 data reuse factor (L=0 -> NLMS; L=1 -> BNLMS; ...)

%   Initializing & Allocating memory:
W=zeros(Flt_Len,(Tran_len+1),ensemble);
W2=zeros(Flt_Len,(Tran_len+1),ensemble);
W3=zeros(Flt_Len,(Tran_len+1),ensemble);
W4=zeros(Flt_Len,(Tran_len+1),ensemble);
W5=zeros(Flt_Len,(Tran_len+1),ensemble);

%   Computing:
for m=1:ensemble

    x1=randn(Tran_len,1);
    if(INPUT_MODE==1)
        b=1;a=[1,-0.9]; %AR(0.9)
    else
        b=[1,0.5,0.81];a=[1,-0.59,0.4]; %ARMA(2,2)
    end
    xout=filter(b,a,x1); % Creating the input signal (normalized)

    % Variable H
    changed_time= Tran_len/2;
    ddout=zeros(length(xout),1);
    dd_reg=zeros(Flt_Len,1);
    for i=1:length(ddout)
        dd_reg=[xout(i);dd_reg(1:end-1)];
        if i<changed_time
            ddout(i)=H1.'*dd_reg;
        else
            ddout(i)=H2.'*dd_reg;
        end
    end

    % Setting variable SNR
    SNRdB = 30*ones(1,8);
    %  SNRdB = [25 15 20 10 25 14 20 18];
    % SNRdB = [30 28 25 22 19 16 13 10];
    %SNRdB = [30 28 25 22 19 16 13 10]-5;
    %SNRdB=SNRdB(end:-1:1);
    varNoise = var(ddout)./(10.^(SNRdB/10));
    ideal_beta = Flt_Len*(1+sqrt(10.^(SNRdB/10)+1))./(10.^(SNRdB/10));
    noise=zeros(1,length(ddout));
    block_size=1e4*2;

    % Setting Impulsive noise
    GINR = 1e-3;
    pb = 0.01;
    sigma = sqrt(varNoise);
    imp=zeros(1,length(ddout));
    gamma_vector=zeros(length(ddout),1);
    SNR_vector=ones(length(ddout),1);
    for i=1:length(SNRdB)
        SNR_vector((i-1)*block_size+1 : i*block_size ) = SNRdB(i);
        noise((i-1)*block_size+1 : i*block_size ) = sqrt(varNoise(i)).*randn(1,block_size);
        imp((i-1)*block_size+1 : i*block_size ) = BG_noise(pb,sigma(i),GINR,block_size);
        gamma_vector((i-1)*block_size+1 : i*block_size ) = ideal_beta(i)*var(xout);
    end

    if(IN_ON==1)
        imp(1:2e4-1)=0;
        imp(2.5e4:3e4-1)=0;
        imp(3.5e4:4e4)=0;
        imp(4e4:5e4)=0;
        imp(5.5e4:6e4)=0;
        imp(6.5e4:end)=0;
        ddout = ddout+imp'+noise';
    else
        ddout = ddout+noise';
    end


    S   =   struct('step',mu,'filterOrderNo',(Flt_Len-1),'initialCoefficients',...
        W(:,1,m),'memoryLength',P);
    S2   =   struct('step',mu,'filterOrderNo',(Flt_Len-1),'initialCoefficients',...
        W2(:,1,m),'memoryLength',P);
    S3   =   struct('step',mu,'filterOrderNo',(Flt_Len-1),'initialCoefficients',...
        W2(:,1,m),'memoryLength',P);

    % Conventional VR with small rho_2 = 9
    [y,e,W(:,:,m),gamma] = APSA_VARgamma_conv(ddout,transpose(xout),S,rho);

    % Conventional VR with large rho_2 = 9e3
    [y5,e5,W5(:,:,m),gamma5] = APSA_VARgamma_conv(ddout,transpose(xout),S,rho*1000);

    % APSA with delta=0.1
    [y3,e3,W3(:,:,m)]=APSA(ddout.',transpose(xout),S,gamma3);

    % APSA with delta 1e5
    [y2,e2,W2(:,:,m)]=APSA(ddout.',transpose(xout),S,gamma2);

    %Proposed method
    [y4,e4,W4(:,:,m),gamma4,delta_gamma4] = VRAPSA_proposed(ddout,transpose(xout),S,rho);

    disp(sprintf('Ensemble: %d \n',m));

end

misaAPSA=zeros(Tran_len+1,ensemble);
misaAPSA2=zeros(Tran_len+1,ensemble);
misaAPSA3=zeros(Tran_len+1,ensemble);
misaAPSA4=zeros(Tran_len+1,ensemble);
misaAPSA5=zeros(Tran_len+1,ensemble);

parfor k = 1: ensemble
    for j = 1:Tran_len
        if j<changed_time
            misaAPSA(j,k) = norm(H1-W (:,j,k))/norm(H1);
            misaAPSA2(j,k) = norm(H1-W2 (:,j,k))/norm(H1);
            misaAPSA3(j,k) = norm(H1-W3 (:,j,k))/norm(H1);
            misaAPSA4(j,k) = norm(H1-W4 (:,j,k))/norm(H1);
            misaAPSA5(j,k) =norm(H1-W5 (:,j,k))/norm(H1); %APA
        else
            misaAPSA(j,k) = norm(H2-W (:,j,k))/norm(H2);
            misaAPSA2(j,k) = norm(H2-W2 (:,j,k))/norm(H2);
            misaAPSA3(j,k) = norm(H2-W3 (:,j,k))/norm(H2);
            misaAPSA4(j,k) = norm(H2-W4 (:,j,k))/norm(H2);
            misaAPSA5(j,k) =norm(H2-W5 (:,j,k))/norm(H2); %APA
        end
    end
    disp(sprintf('AVG is calculated @ ensemble %d \n',k));
end

avg_misaAPSA=20*log10(sum(misaAPSA,2)/ensemble);
avg_misaAPSA2=20*log10(sum(misaAPSA2,2)/ensemble);
avg_misaAPSA3=20*log10(sum(misaAPSA3,2)/ensemble);
avg_misaAPSA4=20*log10(sum(misaAPSA4,2)/ensemble);
avg_misaAPSA5=20*log10(sum(misaAPSA5,2)/ensemble);

figure
plot(avg_misaAPSA,'b','linewidth',1)
hold on
plot(avg_misaAPSA5,'y','linewidth',1)%APA
plot(avg_misaAPSA3,'c','linewidth',1)
plot(avg_misaAPSA2,'g','linewidth',1)
plot(avg_misaAPSA4,'r','linewidth',1)


legend('Conventional VR [7] with small \rho_2','Conventional VR [7] with large \rho_2',sprintf('APSA with \\delta=%1.1e',gamma3),sprintf('APSA with \\delta=%1.1e',gamma2),'Proposed VR-APSA');
xlabel('Time Index (n)'); %
ylabel('NMSD (dB)'); %
title('NMSD learning curves');


figure;plot(20*log10(gamma));hold on;plot(20*log10(gamma5),'b--');plot(20*log10(gamma4),'r');
legend('Conventional VR [7] with small \rho_2','Conventional VR [7] with large \rho_2','Proposed VR-APSA');
title('Evolution of \delta(n)')
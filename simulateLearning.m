%% simulate estimated model
clear


original = csvread('parameters_original_value.csv');
to_be = csvread('parameters_to_be_estimated.csv');
estimated = csvread('estimated_parameters.csv');

param = original;
param(to_be) = estimated;
warning('off','MATLAB:nearlySingularMatrix')

S = 2;
T = 2700;
Learn=1;
Act=0;
lastshock=0;
%parfor s = 1:S
for s = 1:S   
    [Y_real,EPS,EAS,UPS,UAS,WPS,WAS,SPS,SAS,EPLS,EALS,UPLS,UALS,WPLS,WALS,SPLS,SALS,DPS,DAS,NPS,NAS,DIPS,DIAS,DPLS,DALS,NPLS,NALS,DIPLS,DIALS,DPKS,DAKS,NPKS,NAKS,DIKPS,DIKAS,DPKLS,DAKLS,NPKLS,NAKLS,DIKPLS,DIKALS,DPBS,DABS,NPBS,NABS,DIBPS,DIBAS,DPBLS,DABLS,NPBLS,NABLS,DIBPLS,DIBALS,predictions_ws,predictions_wls,predictions_fs,predictions_fls,predictions_ks,predictions_kls,actual_ws,actual_wls,actual_fs,actual_fls,actual_ks,actual_kls,welfare_c,welfare_dc1,welfare_dc2,gperiods,pub_exp_cr,prob_EE,prob_EU,prob_UU,prob_UE,prob_DD,prob_DN,prob_NN,prob_ND,prob_DD_k,prob_DN_k,prob_NN_k,prob_ND_k,price,EXPcontrol,Invent,Assets,baryk,valI,actualEXP, gdp_deflator, Investment,I, consumption, Prod_k, Prod_c, Un, totalDeb, totalDeb_k,stock_bonds,GB,TA,G,wages_t, desired_consumption,rwage,finalshock,et] = learningModel(s,T, param,Learn,Act,lastshock);
    Gov(:,s)=pub_exp_cr;
    DC(:,s)=desired_consumption;
    Y(:,s)=Y_real;
    C(:,s)=consumption;
    Ig(:,s)=I;
    Gov2(:,s)=actualEXP;
    FinShock(s)=finalshock;
    Welfare_c(:,s)=welfare_c;
    Welfare_dc1(:,s)=welfare_dc1;
    Welfare_dc2(:,s)=welfare_dc2;
    PWS(:,:,s)=predictions_ws;
    PWLS(:,:,s)=predictions_wls;
    PFS(:,:,s)=predictions_fs;
    PFLS(:,:,s)=predictions_fls;
    PKS(:,:,s)=predictions_ks;
    PKLS(:,:,s)=predictions_kls;
    AWS(:,:,s)=actual_ws;
    AWLS(:,:,s)=actual_wls;
    AFS(:,:,s)=actual_fs;
    AFLS(:,:,s)=actual_fls;
    AKS(:,:,s)=actual_ks;
    AKLS(:,:,s)=actual_kls;
    eps(:,s)=EPS;
    eas(:,s)=EAS;
    ups(:,s)=UPS;
    uas(:,s)=UAS;
    wps(:,s)=WPS;
    was(:,s)=WAS;
    sps(:,s)=SPS;
    sas(:,s)=SAS;
    epls(:,s)=EPLS;
    eals(:,s)=EALS;
    upls(:,s)=UPLS;
    uals(:,s)=UALS;
    wpls(:,s)=WPLS;
    wals(:,s)=WALS;
    spls(:,s)=SPLS;
    sals(:,s)=SALS;
    dps(:,s)=DPS;
    das(:,s)=DAS;
    nps(:,s)=NPS;
    nas(:,s)=NAS;
    dips(:,s)=DIPS;
    dias(:,s)=DIAS;
    dpls(:,s)=DPLS;
    dals(:,s)=DALS;
    npls(:,s)=NPLS;
    nals(:,s)=NALS;
    dipls(:,s)=DIPLS;
    dials(:,s)=DIALS;
    dpks(:,s)=DPKS;
    daks(:,s)=DAKS;
    npks(:,s)=NPKS;
    naks(:,s)=NAKS;
    dikps(:,s)=DIKPS;
    dikas(:,s)=DIKAS;
    dpkls(:,s)=DPKLS;
    dakls(:,s)=DAKLS;
    npkls(:,s)=NPKLS;
    nakls(:,s)=NAKLS;
    dikpls(:,s)=DIKPLS;
    dikals(:,s)=DIKALS;
    dpbs(:,s)=DPBS;
    dabs(:,s)=DABS;
    npbs(:,s)=NPBS;
    nabs(:,s)=NABS;
    dibps(:,s)=DIBPS;
    dibas(:,s)=DIBAS;
    dpbls(:,s)=DPBLS;
    dabls(:,s)=DABLS;
    npbls(:,s)=NPBLS;
    nabls(:,s)=NABLS;
    dibpls(:,s)=DIBPLS;
    dibals(:,s)=DIBALS;
end


filename = 'learn_noact.mat';
save(filename,'-v7.3')
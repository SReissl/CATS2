
function  [Y_real,EPS,EAS,UPS,UAS,WPS,WAS,SPS,SAS,EPLS,EALS,UPLS,UALS,WPLS,WALS,SPLS,SALS,DPS,DAS,NPS,NAS,DIPS,DIAS,DPLS,DALS,NPLS,NALS,DIPLS,DIALS,DPKS,DAKS,NPKS,NAKS,DIKPS,DIKAS,DPKLS,DAKLS,NPKLS,NAKLS,DIKPLS,DIKALS,DPBS,DABS,NPBS,NABS,DIBPS,DIBAS,DPBLS,DABLS,NPBLS,NABLS,DIBPLS,DIBALS,predictions_ws,predictions_wls,predictions_fs,predictions_fls,predictions_ks,predictions_kls,predictions_bs,predictions_bls,actual_ws,actual_wls,actual_fs,actual_fls,actual_ks,actual_kls,actual_bs,actual_bls,welfare_c,welfare_dc1,welfare_dc2,gperiods,pub_exp_cr,prob_EE,prob_EU,prob_UU,prob_UE,prob_DD,prob_DN,prob_NN,prob_ND,prob_DD_k,prob_DN_k,prob_NN_k,prob_ND_k,price,EXPcontrol,Invent,Assets,baryk,valI,actualEXP, gdp_deflator, Investment,I, consumption, Prod_k, Prod_c, Un, totalDeb, totalDeb_k,stock_bonds,GB,TA,G,wages_t, DC,rwage,finalshock,et] = learningModel(seed, T, par,Learn,Act,lastshock)
tic 

Bk=1;                               %no. of banks
F=100; %200;%                         %no. of firms/capitalists
W=1000; % 3000;%                        %no. of workers
N = 20;% 50;%                          %no. of capital producing firms
z_c=par(1);                              %no. of aplications in consumption good market
z_k = par(2);                         %no. of aplications in capital good market
z_e = par(3);                            %number of job applications  

rng(seed,'twister')

%SET RANDOM NUMBERS
shock_pk = rand(T,N);
shock_p  = rand(T,F);
prob_k = rand(T,F);
permutations_consumption = NaN(T,W+N+F,z_c); %might happen that two are the same...
permutations_un = NaN(T,W,z_e);
permutations_capital = NaN(T,F,z_k);

for tt = 1:T
    for ii = 1:W+N+F
        permutations_consumption(tt,ii,:) = randperm(F,z_c);
    end
    for ii = 1:W
        permutations_un(tt,ii,:) = randperm(F+N,z_e);
    end
    for ii = 1:F
        permutations_capital(tt,ii,:) = randperm(N,z_k);
    end
    
end

seeds_unemployed = randi(2^30,T,1);
seeds_capital = randi(2^30,T,1);
seed_consumption = randi(2^30,T,1);
seeds_surplusk = randi(2^30,T,N);
seeds_surplus = randi(2^30,T,F);

Gov=zeros(1,T+30);

warning('off','MATLAB:nearlySingularMatrix')


coeffsd=rand(F,1);
md=zeros(1,T);
coeffsdk=rand(N,1);
mdk=zeros(1,T);
YPtest=1;
horizon=400;
tau=2/3;
burnin=horizon+1;
r_f=0.01;%   0.015;% 0.005          %general refinancing rate
r_d=r_f/2;
counter=0;
trans_UE=zeros(1,T);
trans_UE(1)=6000;
trans_EU=zeros(1,T);
trans_EU(1)=500;
trans_UU=zeros(1,T);
trans_UU(1)=4000;
trans_EE=zeros(1,T);
trans_EE(1)=9500;
trans_ND=zeros(1,T);
trans_ND(1)=2500;
trans_DN=zeros(1,T);
trans_DN(1)=5000;
trans_DD=zeros(1,T);
trans_DD(1)=5000;
trans_NN=zeros(1,T);
trans_NN(1)=7500;

trans_ND_k=zeros(1,T);
trans_ND_k(1)=6500;
trans_DN_k=zeros(1,T);
trans_DN_k(1)=2000;
trans_DD_k=zeros(1,T);
trans_DD_k(1)=8000;
trans_NN_k=zeros(1,T);
trans_NN_k(1)=3500;

trans_ND_b=zeros(1,T);
trans_ND_b(1)=6000;
trans_DN_b=zeros(1,T);
trans_DN_b(1)=3500;
trans_DD_b=zeros(1,T);
trans_DD_b(1)=6500;
trans_NN_b=zeros(1,T);
trans_NN_b(1)=4000;

prob_EE=zeros(1,T);
prob_EU=zeros(1,T);
prob_UE=zeros(1,T);
prob_UU=zeros(1,T);
prob_EE_s=zeros(1,T);
prob_EU_s=zeros(1,T);
prob_UE_s=zeros(1,T);
prob_UU_s=zeros(1,T);
prob_EE_ls=zeros(1,T);
prob_EU_ls=zeros(1,T);
prob_UE_ls=zeros(1,T);
prob_UU_ls=zeros(1,T);
prob_EE_ns=zeros(1,T);
prob_EU_ns=zeros(1,T);
prob_UE_ns=zeros(1,T);
prob_UU_ns=zeros(1,T);
prob_DD=zeros(1,T);
prob_DN=zeros(1,T);
prob_ND=zeros(1,T);
prob_NN=zeros(1,T);
prob_DD_s=zeros(1,T);
prob_DN_s=zeros(1,T);
prob_ND_s=zeros(1,T);
prob_NN_s=zeros(1,T);
prob_DD_ls=zeros(1,T);
prob_DN_ls=zeros(1,T);
prob_ND_ls=zeros(1,T);
prob_NN_ls=zeros(1,T);
prob_DD_ns=zeros(1,T);
prob_DN_ns=zeros(1,T);
prob_ND_ns=zeros(1,T);
prob_NN_ns=zeros(1,T);
prob_DD_k=zeros(1,T);
prob_DN_k=zeros(1,T);
prob_ND_k=zeros(1,T);
prob_NN_k=zeros(1,T);
prob_DD_ks=zeros(1,T);
prob_DN_ks=zeros(1,T);
prob_ND_ks=zeros(1,T);
prob_NN_ks=zeros(1,T);
prob_DD_kls=zeros(1,T);
prob_DN_kls=zeros(1,T);
prob_ND_kls=zeros(1,T);
prob_NN_kls=zeros(1,T);
prob_DD_kns=zeros(1,T);
prob_DN_kns=zeros(1,T);
prob_ND_kns=zeros(1,T);
prob_NN_kns=zeros(1,T);

prob_DD_b=zeros(1,T);
prob_DN_b=zeros(1,T);
prob_ND_b=zeros(1,T);
prob_NN_b=zeros(1,T);
prob_DD_bs=zeros(1,T);
prob_DN_bs=zeros(1,T);
prob_ND_bs=zeros(1,T);
prob_NN_bs=zeros(1,T);
prob_DD_bls=zeros(1,T);
prob_DN_bls=zeros(1,T);
prob_ND_bls=zeros(1,T);
prob_NN_bls=zeros(1,T);
prob_DD_bns=zeros(1,T);
prob_DN_bns=zeros(1,T);
prob_ND_bns=zeros(1,T);
prob_NN_bns=zeros(1,T);

rw=zeros(1,T);
rw(1)=0.5;
rwage=zeros(1,T);
rwage(1)=rw(1);
rdivs=zeros(F,T);
rdivs(:,1)=0.1;
rdivs_pos=zeros(1,T);
rdivs_k=zeros(N,T);
rdivs_k(:,1)=0.1;
rdivs_kpos=zeros(1,T);
rdivs_b=zeros(1,T);
rdivs_b(1)=0.1;
rd=zeros(1,F);
rdk=zeros(1,N);
yd=zeros(1,F);
ydk=zeros(1,N);
residualC=zeros(1,T);
bonds_real=zeros(1,T);
deficit_real=zeros(1,T);
Ygap=zeros(1,T);
rw_s=zeros(1,T);
rw_ls=zeros(1,T);
rw_ns=zeros(1,T);
rdb_s=zeros(1,T);
rdb_ls=zeros(1,T);
rdb_ns=zeros(1,T);
rd_s=zeros(1,F);
rd_ls=zeros(1,F);
rd_ns=zeros(1,F);
md_s=zeros(1,T);
md_ls=zeros(1,T);
md_ns=zeros(1,T);
rdk_s=zeros(1,N);
rdk_ls=zeros(1,N);
rdk_ns=zeros(1,N);
mdk_s=zeros(1,T);
mdk_ls=zeros(1,T);
mdk_ns=zeros(1,T);

cdemand1=zeros(W+F+N,T);
cdemand2=zeros(W+F+N,T);
cactual=zeros(W+F+N,T);
welfare_dc1=zeros(W+F+N,1);
welfare_dc2=zeros(W+F+N,1);
welfare_c=zeros(W+F+N,1);

EPS=NaN(101,1);
EAS=NaN(101,1);
UPS=NaN(101,1);
UAS=NaN(101,1);
WPS=NaN(101,1);
WAS=NaN(101,1);
SPS=NaN(101,1);
SAS=NaN(101,1);
   
EPLS=NaN(101,1);
EALS=NaN(101,1);
UPLS=NaN(101,1);
UALS=NaN(101,1);
WPLS=NaN(101,1);
WALS=NaN(101,1);
SPLS=NaN(101,1);
SALS=NaN(101,1);
   
DPS=NaN(101,1);
DAS=NaN(101,1);
NPS=NaN(101,1);
NAS=NaN(101,1);
DIPS=NaN(101,1);
DIAS=NaN(101,1);
   
DPLS=NaN(101,1);
DALS=NaN(101,1);
NPLS=NaN(101,1);
NALS=NaN(101,1);
DIPLS=NaN(101,1);
DIALS=NaN(101,1);
   
DPKS=NaN(101,1);
DAKS=NaN(101,1);
NPKS=NaN(101,1);
NAKS=NaN(101,1);
DIKPS=NaN(101,1);
DIKAS=NaN(101,1);
   
DPKLS=NaN(101,1);
DAKLS=NaN(101,1);
NPKLS=NaN(101,1);
NAKLS=NaN(101,1);
DIKPLS=NaN(101,1);
DIKALS=NaN(101,1);

DPBS=NaN(101,1);
DABS=NaN(101,1);
NPBS=NaN(101,1);
NABS=NaN(101,1);
DIBPS=NaN(101,1);
DIBAS=NaN(101,1);
   
DPBLS=NaN(101,1);
DABLS=NaN(101,1);
NPBLS=NaN(101,1);
NABLS=NaN(101,1);
DIBPLS=NaN(101,1);
DIBALS=NaN(101,1);


%set by the sensitivity
xi = par(4);                          %memory parameter human wealth
chi = par(5);                            %fraction of wealth devoted to consumption
q_adj = par(6);                            %quantity adjustment parameter
p_adj = par(7);                            %price adjustment parameter    
mu =  par(8);                           %bank's gross mark-up
eta =par(9);                        %capital depreciation
Iprob=par(10);                         %probability of investing
phi =  par(11);                        %bank's leverage parameter
theta=par(12);                         %rate of debt reimbursment
delta = par(13);                        %memory parameter in the capital utilization rate
alpha = par(14);                              %labour productivity
k = par(15);                            %capital productivity
div =par(16);                            %share of dividends
div_B=0.3;
barX=par(17);%0.85;                          %desired capital utilization
inventory_depreciation = par(18);%0.3;           %rate at which capital firms' inventories depreciate
b1 = par(19);   %-15;
b2 = par(20);   %13;                      %Parameters for risk evaluation by banks
b_k1 = par(21); %-5;
b_k2 = par(22); %5 ;

interest_rate = par(23);
subsidy = par(24);
tax_rate = par(27);

wage_update_up = par(28);
wage_update_down = par(29);
u_target = par(30);

%phillips curve
wb=1.5;                                % initial wage rate


%else government will use unemployment subsidies as line below.
bond_interest_rate = interest_rate;   %% ricordare di aggiungere anche questi con liquidità
unemployment_subsidy_init = subsidy;

G=zeros(1,T);                       %government expenditures
TA=zeros(1,T);                      %government income
GB=zeros(1,T);                      %governament budget  GB = TA - G - EXP -->???
GB(1)=-6440;
EXP = Gov;       % spesa pubblica EROGABILE, update erogata in fondo
actualEXP=EXP;
pub_exp_cr=EXP;
EXP(1,1)=0;              % inizializzazione, vedi sotto
Gshock=1;
EXPcontrol = zeros(1,T);

pub_exp_c= zeros(1,T);    % fu il valore totale EROGATO alle consumption firms 
tax_rate_d = 0;             %taxes on dividends

exp_c=zeros(F,T);         % valore EROGATO individualmente per updating liquidity
quota_exp_c=zeros(F,T);   % quota singola erogabile impresa per totale
quota_exp_c(:,1)=1/F;     % time 1: occhio se cambi iniziando con t=2
public_dem_c=zeros(F,T);  % domanda pubblica per le consumption

bonds = zeros(1,T);
bonds(1)=6440;             
stock_bonds=zeros(1,T);
stock_bonds(1)=6440;

deficitPil=zeros(1,T);
primary_deficit_pil = zeros(1,T);
primary_GB= zeros(1,T);
dividendsB= zeros(1,T);

average_interest_rate=zeros(1,T);
average_interest_rate_k=zeros(1,T);
unfilledVacancies_k=zeros(1,T);
unfilledVacancies=zeros(1,T);
Occ_status=zeros(1,T,W);
divstatus=zeros(1,F);
divstatus_k=zeros(1,N);
divstatus_b=0;
Div_prob=zeros(1,T,F);
Div_probk=zeros(1,T,N);
Div_probb=zeros(1,T);
bankrupt=zeros(1,T,F+N);


%%Bank's Parameters
b = [b1;b2];
b_k = [b_k1;b_k2];
                         
%%Initial conditions
%capital firm
Leff_k = zeros(1,N);        %employees
Y_k =zeros(1,N);
Y_prev_k=3*ones(1,N);       %why is this not initialised to 0?

Y_kd=Y_prev_k;

P_k=3*ones(1,N);
A_k =10*ones(1,N);
liquidity_k = A_k;
De_k=ones(1,N);                 %expected demand of k-Firm
deb_k=zeros(1,N);
price_k=zeros(1,T+1);
price_k(1:2)=mean(P_k);                %capital price index
Q_k=zeros(1,N);                     %sales 
Ftot_k=zeros(1,N);                  %borrowing by k firms

interest_r_k=zeros(1,N);
initialA_k=zeros(1,T);
%firms
value_investments=zeros(1,F);
investment=zeros(1,F);              %pshysical capital acquired in the period
K=10*ones(1,F);
A=10*ones(1,F)+ K*price_k(1);                     %firm equity
liquidity=A-K*price_k(1);           %firm liquid resources

capital_value = K*price_k(1);

P=3*ones(1,F);                        %prices
Y_prev=5*ones(1,F);                   %past production

Yd=Y_prev;

Q=zeros(1,F);                       %sales   
Leff=zeros(1,F);                    %current employees
De=ones(1,F);                       %expected demand --> why is this not set to 5???
deb=zeros(1,F);                     %firm debts

barK=K;                             %long-run K
barYK=Y_prev/k;                     %
x = barYK./barK;                    %this gives initial capacity utilisation of 150? Why would you initialise it like that???
X1=zeros(1,T);                       %Macro-time series for capacity utilisation                        
X2=zeros(1,T);

interest_r=zeros(1,F);
%households
w=zeros(1,W);                       %earned wages
PA=ones(1,W+F+N)+1;                 %household personal assets (saving stock)
Oc=zeros(1,W);                      %employment status: Oc(i)=j --> worker j employed by firm i; if i=0, j is unemployed

totE=zeros(1,T);                    %bank equity time series
E=ones(1,Bk)*3000;                   %bank's equity
loans=0;                            %total loans
totE(1:2)=E;


%macro time series
consumption=zeros(1,T);             %total cunsumption 
price=zeros(1,T+1);
price(1:2)=P(1);                    %consumer price index time series
Un=zeros(1,T);                      %unemployment
wages_t = ones(1,T);                 %wages
dividends=zeros(1,T);               %total dividends

Ftot=zeros(1,F);                    %borrowings
defaults=zeros(1,T);                %number of bankruptcies
profitsB=zeros(T,Bk);               %bank profit time series
unsatisfiedDemand = zeros(1,T);
totK = zeros(1,T);
Investment = NaN(1,T);
inflationRate=zeros(1,T); %%time series of inflation rate
totalDeb=NaN(1,T);
totalDeb_k=NaN(1,T);
defaults_k=zeros(1,T); 

Y_nominal_k = NaN(1,T);
Y_nominal_c = NaN(1,T);
Y_nominal_tot= NaN(1,T);
Y_real =NaN(1,T);
gdp_deflator = NaN(1,T);
I = NaN(1,T);
Prod_c = NaN(1,T);
Prod_k = NaN(1,T);
baryk = NaN(1,T);
Kdes = NaN(1,T);
Kdem = NaN(1,T);
unsatisfiedDemandK = NaN(1,T);
valI=NaN(1,T);
Dexp = NaN(1,T);
Dek=NaN(1,T);
Dek2=NaN(1,T);
YPE=NaN(W,T);
YP_e=zeros(1,T);
YPU=NaN(W,T);
YP_u=zeros(1,T);
YP_d=zeros(1,T);
YP_dk=zeros(1,T);
YP_nd=zeros(1,T);
YP_ndk=zeros(1,T);
YPE(:,1)=rw(1);
YPU(:,1)=rw(1)*subsidy;
shares = NaN(1,T,F);
Assets = NaN(1,T);
Invent = NaN(1,T);
Kdivs_p=NaN(1,T,N);
Cdivs_p=NaN(1,T,F);
net_money=NaN(1,T);
bankruptcy_rate=NaN(1,T);
P_lower=NaN(1,T);
Kdem2=NaN(1,T);
DC=NaN(1,T);
Growth=NaN(1,T);
credit_mismatch=NaN(1,T);
Occ_status_prev=NaN(T,W);
%%%END OF MODEL INITIALIZATION

% cons_budget = ones(1,F+N+W)*1;
dividends_income_k = zeros(1,N);
dividends_income   = zeros(1,F);
permanent_income = ones(F+W+N,1)./price(1);
%%%HEREAFTER THE MODEL STARTS RUNNING FOR T PERIODS

money(1:2) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)+totE(1);  %What is this supposed to be? Why is bank's equity part of "money"? Why is debt subtracted???
learn=Learn;
learn_act=Act;
endshock=lastshock;
finalshock=randi([2520 2620],1);
qshocks=0;

%%%%%%%%%%%%%%%%%%
 for t=2:T
     disp(t)
     if learn==1 && t>=700 && t<finalshock-15
         if t/5==round(t/5) && t/10~=round(t/10)
             per=t+5;
             pub_exp_cr(per)=normrnd(0.1*alpha*1000,0.5);
             EXP(per)=pub_exp_cr(per);
             if t>=1500
             qshocks=qshocks+1;
             end
         end
     end
     
     if endshock==1 && t==finalshock-5
        per=finalshock;
        pub_exp_cr(per)=normrnd(0.1*alpha*1000,0.5);
        qshocks=qshocks+1;
     end
     
     gperiods=find(pub_exp_cr>0);
     
     currentshock=gperiods(gperiods>=(t-1));
     if isempty(currentshock)
        currentshock=T+20;
     end
     
     if A<=0                                                             %if all the firms are defaulted, then stop simulation 
        break
     end
    
     EXP(t)=EXP(t)*price(t);
     if Gshock==1
     gexp=transpose(quota_exp_c(:,t-1).*(EXP(t))./P(:));
     else
     gexp=0;
     end
     
    stock=Y_prev-Yd;
    Invent(t) = sum(stock)*3;
        for i=1:F
            if stock(i)<=0 && P(i)>=price(t)
                if Gshock==1 && gexp(i)>0
                    De(i) = Y_prev(i) + (-stock(i))*q_adj + gexp(i);
                else
                    De(i) = Y_prev(i) + (-stock(i))*q_adj;
                end
            elseif stock(i)<=0 && P(i)<price(t)
                P(i)=P(i)*(1+shock_p(t,i)*p_adj);
                if Gshock==1 && gexp(i)>0
                    De(i)=Y_prev(i);
                else
                    De(i)=Y_prev(i);
                end
            elseif stock(i)>0 &&P(i)>price(t)
                P(i)=P(i)*(1-shock_p(t,i)*p_adj);
                if Gshock==1 && gexp(i)>0
                    De(i)= Y_prev(i);
                else
                    De(i)= Y_prev(i);
                end
            elseif stock(i)>0 && P(i)<=price(t)
                if Gshock==1 && gexp(i)>0
                    De(i) = Y_prev(i) - stock(i)*q_adj + gexp(i);
                else
                    De(i) = Y_prev(i) - stock(i)*q_adj;
                end
            end    
            if De(i)<alpha                                                  %in order to hire at least one worker
                De(i)=alpha;    
            end
        end
    Dexp(t)= sum(De)*3;
    
    %CAPITAL PRODUCTION DECISION
    inv_dep=inventory_depreciation*Y_k;
    inv_dep=inv_dep.*P_k;
    inventory_k=(1-inventory_depreciation)*Y_k;
    invent_start=inventory_k.*P_k;
   
    stock_k = Y_prev_k-Y_kd; 
    for i=1:N
        if stock_k(i)<=0 && P_k(i)>=price_k(t)
             De_k(i) = (Y_prev_k(i) + (-stock_k(i))*q_adj)-inventory_k(i);
        elseif stock_k(i)<=0 && P_k(i)<price_k(t)
            De_k(i)=Y_prev_k(i);
            P_k(i)=P_k(i)*(1+shock_pk(t,i)*p_adj);
            
        elseif stock_k(i)>0 && P_k(i)>price_k(t)
            De_k(i)=Y_prev_k(i); 
            P_k(i)=P_k(i)*(1-shock_pk(t,i)*p_adj);
           
        elseif stock_k(i)>0 && P_k(i)<=price_k(t)
             De_k(i) = (Y_prev_k(i) - (stock_k(i))*q_adj)-inventory_k(i);
        end    
        
        if De_k(i)<alpha                                                  %in order to hire at least one worker
            De_k(i)=alpha;    
        end
    end    
    
   
    
    %% investments

    prob=prob_k(t,:);
    K_dem=zeros(1,F);
    K_des = barYK/barX;
    Kdes(t)= sum(K_des);
    depreciation = K*eta*1/(Iprob);
    K_dem(prob<Iprob)= K_des(prob<Iprob) - K(prob<Iprob)+depreciation(prob<Iprob);
    K_dem(K_dem<0)=0;
    Kdem2(t) = sum(K_dem);
    
    %labour requirement (consumption good)
    Ld=min(ceil(De./alpha), ceil(K.*k/alpha));  
    wages=wb*Ld;   
    p_k=price_k(t);
    
   
    %labour requirement (capital good)
    Ld_k=ceil(De_k./alpha);
    wages_k = wb*Ld_k;

%% CREDIT MARKET OPENS
             
    
    Ftot(:)=0;
    Ftot_k(:)=0;
   
    
    %compute financial gap
    B=wages+K_dem*p_k-liquidity;                                        %financial gap          
    B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
    B(B<0)=0;
    B_k(B_k<0)=0;
    
    lev=(deb+B)./(A+B+deb);                                             %leverage
    lev_k=(deb_k+B_k)./(A_k+B_k+deb_k);
    
    
    loan_applications=find(B>0);                                        %only firms with positive credit demand will apply for an additional loan
    loan_applications_k=find(B_k>0); 
    %evaluate bankruptcy probability and expected survival
    pr(:,1)= exp(b(1)+b(2)*lev)./(1+exp(b(1)+b(2)*lev)); %glmval(b,lev,'logit');  %zeros(F,1);%                           %banks evaluate bankruptcy probability of each firm, given the estimated
    pr_k(:,1)=exp(b_k(1)+b_k(2)*lev_k)./(1+exp(b_k(1)+b_k(2)*lev_k));%glmval(b_k,lev_k,'logit');%zeros(N,1);%                        %parameters (b,b_k)and computed leverage                                                                        
                                                           
    Xi=(1-(1-theta).^(1+1./pr))./(1-(1-theta));                         %this is Xi in the paper
    Xi_k=(1-(1-theta).^(1+1./pr_k))./(1-(1-theta));
    %proposed rate depends on the estimated bankruptcy probability
    proposed_rate=mu*((1+r_f/theta)./Xi - theta)';
    proposed_rate_k=mu*((1+r_f/theta)./Xi_k - theta)';
    
    

    %for each firm the bank computes the maximum loan and gives loans up to the maximum amount 
    for i=loan_applications
        credit=B(i);
        %the bank gives a maximum credit depending
        %on  maximum expected loss
        
        maxL = (phi*totE(t)-pr(i)*deb(i))/pr(i);
        maxL=max(0,maxL); %maxL never negative        
        credit=min(credit,maxL); %%credit given to firm i           

        
        deb_0=deb(i);
        deb(i)=deb(i)+credit;                                           %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot(i)=credit;                                                 %record flow of new credit for firm i
        %compute new average interest rate
        if deb(i)>0
            interest_r(i)=(deb_0*interest_r(i)+proposed_rate(i)*credit)/deb(i);
        end
    
    end       
    %weighted average interest rate
     average_interest_rate(t)=sum(interest_r.*deb)/sum(deb);

     
    %mutatis mutandis for capital firms
    for i=loan_applications_k
        credit=B_k(i);
        maxL = (phi*totE(t)-pr_k(i)*deb_k(i))/pr_k(i);
        maxL=max(0,maxL);
        credit=min(credit,maxL);
         
        deb_0=deb_k(i);
        deb_k(i)=deb_k(i)+credit;                                       %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot_k(i)=credit;                                               %record flow of new credit for firm i
        if deb_k(i)>0
            interest_r_k(i)=(deb_0*interest_r_k(i)+proposed_rate_k(i)*credit)/deb_k(i);
        end
       
    end 
    average_interest_rate_k(t)=sum(interest_r_k.*deb_k)/sum(deb_k);
    credit_mismatch(t) = sum(B(loan_applications))+sum(B_k(loan_applications_k)) - sum(Ftot)-sum(Ftot_k);

    %%CREDIT MARKET CLOSES    

    %% JOB MARKET OPENS
    Occ_status_prev(t,:)=Oc;
    Occ_status_prev(Occ_status_prev>0)=1;
    
    %determine desired labour and vacancies given available liquidity
    Ld_k=min(Ld_k, (Ftot_k+liquidity_k)/wb);
    Ld_k=floor(Ld_k);
    Ld_k(Ld_k<1)=1;
    vacancies_k=Ld_k-Leff_k;
    surplus_k=find(vacancies_k<0);                                          %firms with too many workers
    
    %%CONSUPTION GOOD

    %%re-define labour demand given available liquidity
    Ld=min(Ld,(Ftot+liquidity)/wb);                                     %demand (stock)     
    Ld=floor(Ld);
    Ld(Ld<1)=1; %%since there is the capital the firms can have positive equity and negative liquidity, in this latter case Ld would be negative, which is impossible
    vacancies=Ld-Leff;                                                  %definitive labour demand (flow)    
    
    %%JOB MARKET OPENS
    surplus=find(vacancies<0);                                          %firms with too many workers
    

    for i=surplus_k
        workforce_k=find(Oc==F+i);
        rng(seeds_surplusk(t,i))%pick all firm i's workers
        f_k=randperm(length(workforce_k));
        f_k=f_k(1:-vacancies_k(i));                                     %take randomly "-vacancies(i)" workers and fire them
        fired_k=workforce_k(f_k);  
        Oc(fired_k)=0;
        w(fired_k)=0;
        Leff_k(i)=Ld_k(i);                                              %update no. of workers
     end


 %firms with excess workforce fire
    for i=surplus   
        workforce=find(Oc==i);  
        rng(seeds_surplus(t,i))%pick all firm i's workers
        f=randperm(length(workforce));
        f=f(1:-vacancies(i));                                           %take randomly "-vacancies(i)" workers and fire them
        fired=workforce(f);  
        Oc(fired)=0;
        w(fired)=0;
        Leff(i)=Ld(i);                                                  %update no. of workers
    end
    
%% UNEMPLOYED WORKERS LOOK FOR A JOB
    
    unemployed=find(Oc==0);
    
    
    rng(seeds_unemployed(t))
    vec = randperm(length(unemployed));
    for un=vec      
        j=unemployed(un);                                               %randomly pick an unemployed worker

        Z_e = permutations_un(t,j,:); 
        flag=1;
        
        while (Oc(j)==0 && flag<=z_e)                              %continue searching until you are unemployed and didn't search at all available firms
            f=Z_e(flag);                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            if f>F %selected firm is a capital firm
                if vacancies_k(f-F)>0                                    %if the selected firm has an open vacancy, take the job
                Oc(j)=f;                                                 %update employed status
                w(j)=wb;                                                 %salary   
                Leff_k(f-F)=Leff_k(f-F)+1;                               %firm's workforce   
                vacancies_k(f-F)=vacancies_k(f-F)-1;                     %firm's vacancies   
                end
            else %selected firm is a consuption firm
                if vacancies(f)>0                                        %if the selected firm has an open vacancy, take the job
                Oc(j)=f;
                w(j)=wb;
                Leff(f)=Leff(f)+1;
                vacancies(f)=vacancies(f)-1;
                end                
            end
            flag=flag+1;                                                %increase counter
        end 
        
    end  
    
    %%JOB MARKET CLOSES
    
    unfilledVacancies_k(t) = sum(vacancies_k(vacancies_k>0)); 
    unfilledVacancies(t) = sum(vacancies(vacancies_k>0));
    
    %% production
    %produce capital
    Dek(t)=sum(De_k);
    Dek2(t)=sum(Leff_k*alpha);
    Y_k=min(De_k,Leff_k*alpha);                              
    Y_prev_k=Y_k;
    Prod_k(t)=sum(Y_prev_k);
  
    Y_k = Y_k+inventory_k;                                   %Y_k is increased by the inventories (capital good is durable)   
    
    %produce consuption
    Yp = min(Leff*alpha, K*k);                               %production frontier given available resources (Leontieff)
    Y=min(De,Yp);                                            %actual production
    
    Ygap(t)=sum(Yp)-sum(Y);
    Y_prev=Y;
    Prod_c(t)=sum(Y_prev);
    X1(t)=sum(Y_prev./sum(Y_prev).*Y_prev./(K*k));
    X2(t)=sum(Y_prev)/(sum(K)*k);
    totK(t)=sum(K);
 
    interests=interest_r.*deb; 
    interests_k=interest_r_k.*deb_k;
    wages=wb*Leff;
    wages_k=wb*Leff_k;

    
    %minimum prices
    
    Pl = 1.01*(wages+interests)./Y_prev;  
    Pl(isinf(Pl)) = 0;
    P(P<Pl) = Pl(P<Pl);

    
    %minimum prices capital
    Pl_k = 1.01*(wages_k+interests_k)./(Y_prev_k);
    Pl_k(isinf(Pl_k)) = 0;
    P_k(P_k<Pl_k) = Pl_k(P_k<Pl_k);
    

  
    %%CAPITAL GOODS MARKET OPENS
    capital_budget=max(0,Ftot+liquidity-Leff*wb);                              %amount of liquidity available to buy capital goods
    %capital_budget=Ftot+liquidity-Leff*wb;                              %amount of liquidity available to buy capital goods
    capdem=min(K_dem,capital_budget./price_k(t));
    Kdem(t) = sum(min(K_dem,capital_budget./price_k(t)));
    capital_demanders=find(K_dem>0);                                    
    investment(:)=0;
    value_investments(:)=0;
    Q_k(:)=0;
    Y_kd=zeros(1,N); 
    
    rng(seeds_capital(t))
    for firm=randperm(length(capital_demanders))
        j=capital_demanders(firm);
        Z = permutations_capital(t,j,:);
        PZ_k=P_k(Z);
        [~,order]=sort(PZ_k);
        flag=1;

        while (capital_budget(j)>0 && K_dem(j)>0 && flag<=z_k) 
            best=Z(order(flag));
            Y_kd(best)=Y_kd(best)+min(capital_budget(j)/P_k(best), K_dem(j));
                if Y_k(best)>0                                                %buy if 'best' firm has still positive stocks 
                    pk=P_k(best);
                    budget=min(capital_budget(j), K_dem(j)*pk);
                    if Y_k(best) > budget/pk           
                        Y_k(best)=Y_k(best)-budget/pk;                   %reduce stocks
                        Q_k(best)=Q_k(best)+budget/pk;                   %update sales
                        K_dem(j)=K_dem(j)-budget/pk;
                        capital_budget(j) = capital_budget(j)-budget;
                        liquidity(j)=liquidity(j)-budget;
                        investment(j)=investment(j)+budget/pk;
                        value_investments(j)= value_investments(j)+budget;
                                                           %j spends all its budget  
                    elseif  Y_k(best) <= budget/pk  
                        K_dem(j)=K_dem(j)-Y_k(best);
                        capital_budget(j)=capital_budget(j) - Y_k(best)*pk;
                        liquidity(j)=liquidity(j)- Y_k(best)*pk;
                        Q_k(best)=Q_k(best)+Y_k(best);    
                        investment(j)=investment(j)+Y_k(best);
                        value_investments(j)= value_investments(j)+Y_k(best)*pk;
                        Y_k(best)=0;    
                    end    
                end
                flag=flag+1;                                                %increase counter
        end  
    end
    
    unsatisfiedDemandK(t) = sum(max(0,capdem-investment));
    valI(t)=sum(value_investments);
    I(t) = sum(investment);

    %%CAPITAL GOOD MARKET CLOSES
    
 
    %%CONSUMPTION GOOD MARKET OPENS
    w(w>0) = wb;
    %taxes on wages
    wn=w*(1-tax_rate);  
    TA(t) = TA(t)+sum(w*tax_rate);
    %workers receive wages
    interest_deposits=sum(r_d*PA(1:W));
    PA(1:W)=PA(1:W)+r_d*PA(1:W);
    if YPtest==0 
        PA(1:W)=PA(1:W)+wn;
    end
    
    %%% SET UNEMPLOYMENT SUBSIDY TO MEET THE PUBLIC
    %%% DEFICIT RATIO REQUIREMENT

    unemployed=find(Oc==0);
    
    unemployment_subsidy = unemployment_subsidy_init;

    
    if YPtest==0
        PA(unemployed)=PA(unemployed)+unemployment_subsidy*wb;
    end
    G(t)=unemployment_subsidy*wb*length(unemployed);
    
    
    workers_income = wn;
    workers_income(unemployed) = unemployment_subsidy*wb;
    
    Occ_status(:,t,:)=Oc;
    Occ_status(Occ_status>0)=1;
    

if YPtest == 1
    for j=1:W
    if Occ_status(1,t,j)==1 && Occ_status_prev(t,j)==0
        trans_UE(t)=trans_UE(t)+1;
    end
    if Occ_status(1,t,j)==0 && Occ_status_prev(t,j)==1
        trans_EU(t)=trans_EU(t)+1;
    end
    if Occ_status(1,t,j)==0 && Occ_status_prev(t,j)==0
        trans_UU(t)=trans_UU(t)+1;
    end
    if Occ_status(1,t,j)==1 && Occ_status_prev(t,j)==1
        trans_EE(t)=trans_EE(t)+1;
    end
    end
    
    if sum(trans_EE(1:t))+sum(trans_EU(1:t))>0
        if t<burnin
        prob_EE(t)=sum(trans_EE(1:t))/(sum(trans_EE(1:t))+sum(trans_EU(1:t)));
        prob_EU(t)=sum(trans_EU(1:t))/(sum(trans_EE(1:t))+sum(trans_EU(1:t)));
        else
        prob_EE(t)=sum(trans_EE((t-(horizon-1)):t))/(sum(trans_EE((t-(horizon-1)):t))+sum(trans_EU((t-(horizon-1)):t)));
        prob_EU(t)=sum(trans_EU((t-(horizon-1)):t))/(sum(trans_EE((t-(horizon-1)):t))+sum(trans_EU((t-(horizon-1)):t)));
        end
    else
    prob_EE(t)=prob_EE(t-1);
    prob_EU(t)=prob_EU(t-1);
    end

    if sum(trans_UE(1:t))+sum(trans_UU(1:t))>0
        if t<burnin
        prob_UU(t)=sum(trans_UU(1:t))/(sum(trans_UE(1:t))+sum(trans_UU(1:t)));
        prob_UE(t)=sum(trans_UE(1:t))/(sum(trans_UE(1:t))+sum(trans_UU(1:t)));
        else
        prob_UU(t)=sum(trans_UU((t-(horizon-1)):t))/(sum(trans_UE((t-(horizon-1)):t))+sum(trans_UU((t-(horizon-1)):t)));
        prob_UE(t)=sum(trans_UE((t-(horizon-1)):t))/(sum(trans_UE((t-(horizon-1)):t))+sum(trans_UU((t-(horizon-1)):t)));   
        end
    else
    prob_UU(t)=prob_UU(t-1);
    prob_UE(t)=prob_UE(t-1);
    end 
    
    
    rwage(t)=(wb*(1-tax_rate))/price(t);
    
    if learn==1 && t>=1000 && t<finalshock+1
        
            shocks=find(pub_exp_cr(1:T)>0);
            lshocks=shocks+1;

            shocks=shocks(shocks>=(t-(horizon-1)) & shocks<=(t));
            lshocks=lshocks(lshocks>=(t-horizon) & lshocks<=(t));
            nshocks=(t-(horizon-1)):t;
            nshocks=setdiff(nshocks,shocks);
            nshocks=setdiff(nshocks,lshocks);


            EE1=trans_EE(shocks);
            EU1=trans_EU(shocks);
            UE1=trans_UE(shocks);
            UU1=trans_UU(shocks);

            EE2=trans_EE(lshocks);
            EU2=trans_EU(lshocks);
            UE2=trans_UE(lshocks);
            UU2=trans_UU(lshocks);

            EE3=trans_EE(nshocks);
            EU3=trans_EU(nshocks);
            UE3=trans_UE(nshocks);
            UU3=trans_UU(nshocks);

            prob_EE_s(t)=sum(EE1)/(sum(EE1)+sum(EU1));
            prob_EU_s(t)=sum(EU1)/(sum(EE1)+sum(EU1));
            prob_UE_s(t)=sum(UE1)/(sum(UE1)+sum(UU1));
            prob_UU_s(t)=sum(UU1)/(sum(UE1)+sum(UU1));

            prob_EE_ls(t)=sum(EE2)/(sum(EE2)+sum(EU2));
            prob_EU_ls(t)=sum(EU2)/(sum(EE2)+sum(EU2));
            prob_UE_ls(t)=sum(UE2)/(sum(UE2)+sum(UU2));
            prob_UU_ls(t)=sum(UU2)/(sum(UE2)+sum(UU2));

            prob_EE_ns(t)=sum(EE3)/(sum(EE3)+sum(EU3));
            prob_EU_ns(t)=sum(EU3)/(sum(EE3)+sum(EU3));
            prob_UE_ns(t)=sum(UE3)/(sum(UE3)+sum(UU3));
            prob_UU_ns(t)=sum(UU3)/(sum(UE3)+sum(UU3));

            rw_s(t)=mean(rwage(shocks));
            rw_ls(t)=mean(rwage(lshocks));
            rw_ns(t)=mean(rwage(nshocks));

    elseif t>=finalshock+1
            prob_EE_s(t)=prob_EE_s(t-1);
            prob_EU_s(t)=prob_EU_s(t-1);
            prob_UE_s(t)=prob_UE_s(t-1);
            prob_UU_s(t)=prob_UU_s(t-1);

            prob_EE_ls(t)=prob_EE_ls(t-1);
            prob_EU_ls(t)=prob_EU_ls(t-1);
            prob_UE_ls(t)=prob_UE_ls(t-1);
            prob_UU_ls(t)=prob_UU_ls(t-1);

            prob_EE_ns(t)=prob_EE_ns(t-1);
            prob_EU_ns(t)=prob_EU_ns(t-1);
            prob_UE_ns(t)=prob_UE_ns(t-1);
            prob_UU_ns(t)=prob_UU_ns(t-1);

            rw_s(t)=rw_s(t-1);
            rw_ls(t)=rw_ls(t-1);
            rw_ns(t)=rw_ns(t-1);
    end
    
    
    if t < burnin
        rw(t)=sum(rwage(1:t))/t;
        rs=rw(t)*unemployment_subsidy;
    else
        rw(t)=sum(rwage(t-(horizon-1):t))/horizon;
        rs=rw(t)*unemployment_subsidy;
    end
    
    if t<burnin+1 && t>3
    depend=rwage(2:t)-rw(t);
    ind=depend(1:t-2);
    depend=depend(2:t-1);
    end
    if t>=burnin+1
    depend=rwage(t-horizon:t)-rw(t);
    ind=depend(1:horizon-1);
    depend=depend(2:horizon);
    end
    if t>10
    coeffw=inv(ind*ind')*ind*depend';
    else
    coeffw=0.95;
    end
    coeffw=min(coeffw,0.999);
    if t<burnin
    pay1=rwage(t)*ones(2000,1);
    pay2=rwage(t)*unemployment_subsidy*ones(2000,1);
    else
    pay1=rw(t)*ones(2000,1);
    pay2=rs*ones(2000,1);
    pay1(1)=rwage(t);
    pay2(1)=rwage(t)*unemployment_subsidy;
    q=1;
    while abs(pay1(q)-rw(t))>0.0005
       q=q+1;
       pay1(q)=rw(t)+coeffw*(pay1(q-1)-rw(t));
       pay2(q)=rs+coeffw*(pay2(q-1)-rs);
    end
    end
    pay1(:)=pay1(:);
    pay2(:)=pay2(:);
    if learn_act && t>=1500
    shockperiods=gperiods-t+1;
    lshockperiods=shockperiods+1;
    shockperiods=shockperiods(shockperiods>0);
    lshockperiods=lshockperiods(lshockperiods>0);
    probs11=prob_EE_ns(t)*ones(2000,1);
    probs12=prob_EU_ns(t)*ones(2000,1);
    probs21=prob_UE_ns(t)*ones(2000,1);
    probs22=prob_UU_ns(t)*ones(2000,1);
    probs11(shockperiods)=prob_EE_s(t);
    probs12(shockperiods)=prob_EU_s(t);
    probs21(shockperiods)=prob_UE_s(t);
    probs22(shockperiods)=prob_UU_s(t);
    probs11(lshockperiods)=prob_EE_ls(t);
    probs12(lshockperiods)=prob_EU_ls(t);
    probs21(lshockperiods)=prob_UE_ls(t);
    probs22(lshockperiods)=prob_UU_ls(t);
    pay1(:)=pay1(:)-rw(t)+rw_ns(t);
    pay2(:)=pay2(:)-rs+unemployment_subsidy*rw_ns(t);
    pay1(shockperiods)=pay1(shockperiods)-rw_ns(t)+rw_s(t);
    pay2(shockperiods)=pay2(shockperiods)-unemployment_subsidy*rw_ns(t)+unemployment_subsidy*rw_s(t);
    pay1(lshockperiods)=pay1(lshockperiods)-rw_ns(t)+rw_ls(t);
    pay2(lshockperiods)=pay2(lshockperiods)-unemployment_subsidy*rw_ns(t)+unemployment_subsidy*rw_ls(t);
    pay1(1)=rwage(t);
    pay2(1)=rwage(t)*unemployment_subsidy;
    permw=projectIncome(probs11,probs12,probs21,probs22,pay1,pay2,r_d);
    else
    permw=projectIncome(prob_EE(t),prob_EU(t),prob_UE(t),prob_UU(t),pay1,pay2,r_d);
    end
    YPE(:,t)=(PA(1:W)/price(t)+permw(1))/(1/r_d);
    YPU(:,t)=(PA(1:W)/price(t)+permw(2))/(1/r_d);
    
    YP_e(t)=permw(1);
    YP_u(t)=permw(2);
    
    if t==currentshock-1 && t>=1500 && learn_act
       pay1w=pay1;
       pay2w=pay2;
    end
    
      
end


    if YPtest==0
        income =  ([workers_income,dividends_income,dividends_income_k]')./price(t);
        permanent_income = permanent_income*xi + (1-xi)*income;
    else
        wi=zeros(1,W);
        for i=1:W
            if Oc(i)>0
                wi(i)=YPE(i,t);
            else
                wi(i)=YPU(i,t);
            end
        end
       income = [wi,yd,ydk]';
       permanent_income = income;
    end
    
    if YPtest==1
        PA(1:W)=PA(1:W)+wn;
        PA(unemployed)=PA(unemployed)+unemployment_subsidy*wb;
    end

    if YPtest==0
        target = 1*permanent_income' + chi*PA./price(t) ; %0.05
    else
        target = 1*permanent_income';
    end
    
    cdemand1(:,t)=target;
    cons_budget = target.*price(t);
    cons_budget = min(PA,cons_budget);
    PA=PA-cons_budget;        
    consumers=find(cons_budget>0);
    cdemand2(:,t)=(cons_budget)./price(t);
    DC(t) = (sum(cons_budget)./price(t))*3;

    Q(:)=0;
    Yd=zeros(1,F);
   %search and matching starts

    C = zeros(1,W+F+N);
    
    
    rng(seed_consumption(t))
    vec = randperm(length(consumers));
    for wor=vec     
        j=consumers(wor);                                               %randomly pick a consumer
        Z = permutations_consumption(t,j,:);
        PZ=P(Z);  
        [~,order]=sort(PZ);                                            %sort prices in ascending order
        flag=1;
        
        while (cons_budget(j)>0 && flag<=z_c)                              %continue buying till budget is positive and there are firms available
            best=Z(order(flag));                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            Yd(best)=Yd(best)+cons_budget(j)/P(best);
            if Y(best)>0                                                %buy if 'best' firm has still positive stocks 
                p=P(best);
                if Y(best) > cons_budget(j)/p           
                    Y(best)=Y(best)-cons_budget(j)/p;                   %reduce stocks
                    Q(best)=Q(best)+cons_budget(j)/p; 
                    C(j) = C(j)+cons_budget(j)/p; %update sales
                    consumption(t)=consumption(t)+cons_budget(j)/p;     
                    cons_budget(j)=0;                                   %j spends all its budget  
                elseif  Y(best) <= cons_budget(j)/p  
                    cons_budget(j)=cons_budget(j)- Y(best)*p;   
                    Q(best)=Q(best)+Y(best);
                    C(j) = C(j)+Y(best); 
                    consumption(t)=consumption(t)+Y(best);
                    Y(best)=0;    
                end    
            end
            flag=flag+1;                                                %increase counter
       end 
        
    end    
    
    cactual(:,t)=C(:);
    
    consumption(t)=consumption(t)*3;
    unsatisfiedDemand(t) = sum(cons_budget);            
    PA=PA+cons_budget;
    
    %%%INTRO SPESA PUBBLICA
    pub_exp_c(t)= 1*EXP(t);
    public_dem_c(:,t)=quota_exp_c(:,t-1).*(pub_exp_c(t)./P');  
    EXPcontrol(t)=sum(public_dem_c(:,t));
    for i=1:F
        %Yd(i)=Yd(i)+public_dem_c(i,t);
        if Y(i)>= public_dem_c(i,t)
           Y(i)=Y(i)-public_dem_c(i,t);   
           Q(i)= Q(i)+public_dem_c(i,t);   
        else 
           Q(i)=Q(i)+Y(i);
           public_dem_c(i,t)=Y(i);              
           Y(i)= 0;      
        end                                     
        %Yd(i)=Yd(i)+public_dem_c(i,t);
    end
    exp_c(:,t)= public_dem_c(:,t).*P';     
    pub_exp_c(t)=sum(exp_c(:,t));
    EXP(t)=pub_exp_c(t);
    actualEXP(t)=pub_exp_c(t)/price(t);
    residualC(t)=sum(Y);
    
%%CONSUMPTION GOOD MARKET CLOSES
    %Capital price index
    if sum(Q_k)>0 
        share_k=(Q_k)/sum(Q_k);
    else
        share_k=ones(1,N)/N;
    end
    
    price_k(t+1)=P_k*share_k';

%%ACCOUNTING

    %capital capacity update
    barYK = delta*barYK + (1-delta)*Y_prev/k;
    baryk(t) = mean(barYK);
    %capital depreciation
    dep = eta*Y_prev/k; %%only the used capital is depreciating
    %capital reduced by depreciation
    K = K - dep;
    %capital increased by bought capital 
    K = K + investment;

    %update capital value in the book
    depreciation_value =  (dep).*capital_value./(K+dep-investment);
    capital_value = capital_value - depreciation_value + value_investments;
    
    %firm revenues
    RIC=P.*Q;   
    %firm update their liquidity, pay wages, interests and installments
    liquidity=liquidity+RIC+Ftot-wages-interests-theta*deb; %%investments already paid
    loans=loans-sum(deb)*theta;
    deb=(1-theta)*deb;
    
    %consumption firm profits
    pi = RIC-wages-interests - depreciation_value;
    %equity law of motion!
    A = A+pi;
    
    %Capital producer accounting
    RIC_k = P_k.*Q_k;
    %update liquidity
    liquidity_k = liquidity_k + RIC_k + Ftot_k-wages_k-interests_k-theta*deb_k;
    %update loans
    loans=loans-sum(deb_k)*theta;
    deb_k=(1-theta)*deb_k;
    %profits
    invent_end=Y_k.*P_k;
    pi_k=RIC_k+(invent_end-invent_start)-wages_k-interests_k-inv_dep;
    
    totalDeb(t)=sum(deb);
    totalDeb_k(t)=sum(deb_k);
    
    %dividends
    divstatus_prev=divstatus;
    divstatus_prev(divstatus_prev>0)=1;
    pospi=find(pi>0);                                                   %pick firms with positive profits
    Div_prob(:,t,:)=pi;
    Div_prob(Div_prob<0)=0;
    Div_prob(Div_prob>0)=1;
    dividends_income(:)=0;
    for i=pospi
        if liquidity(i)>0
        di=min(div*pi(i),liquidity(i));  %%dividends                                              %compute dividends paid by firm i
        divi=di*(1-tax_rate_d); %dividends after taxes
        TA(t)=TA(t)+di*tax_rate_d;
        PA(W+i)=PA(W+i)+divi;                                          %dividends paid to firm i owner
        dividends_income(i) = divi;
        liquidity(i)=liquidity(i)-di;
        A(i)=A(i)-di;
        dividends(t)=dividends(t)+di; %lordi
        pi(i)=pi(i)-divi;
        end
    end
    interest_deposits=interest_deposits+sum(r_d*PA(W+1:W+F));
    PA(W+1:W+F)=PA(W+1:W+F)+r_d*PA(W+1:W+F);
    
    divstatus_prev_k=divstatus_k;
    divstatus_prev_k(divstatus_prev_k>0)=1;
    pospi_k=find(pi_k>0); 
    Div_probk(:,t,:)=pi_k;
    Div_probk(Div_probk<0)=0;
    Div_probk(Div_probk>0)=1;
    dividends_income_k(:)=0;%pick firms with positive profits
    for i=pospi_k
        if liquidity_k(i)>0
        di=min(div*pi_k(i),liquidity_k(i));   
        divi=di*(1-tax_rate_d); %dividends after taxes
        TA(t)=TA(t)+di*tax_rate_d;%compute dividends paid by firm i
        PA(W+F+i)=PA(W+F+i)+divi;                                          %dividends paid to firm i owner
        dividends_income_k(i)=divi;
        liquidity_k(i)=liquidity_k(i)-di;
        dividends(t)=dividends(t)+di;
        pi_k(i)=pi_k(i)-divi;
        end
    end
    interest_deposits=interest_deposits+sum(r_d*PA(W+F+1:W+F+N));
    PA(W+F+1:W+F+N)=PA(W+F+1:W+F+N)+r_d*PA(W+F+1:W+F+N);
    
   
    A_k=liquidity_k+Y_k.*P_k-deb_k;

    
    %replacement of bankrupted consumption firms 
    piB=0;                                                              %reset bank's profits
       
    
    %% time series (before bankruptcies)
    %inflation rate
    if sum(Q)>0 
        share=(Q)/sum(Q);
    else
        share=ones(1,F)/F;
    end
    P_lower(t)= Pl*share';
    RPI=P*share';                                                       %retail price index
    infla_rate=RPI/price(t);
    %price=[price RPI];
    price(t+1)=RPI;
    inflationRate(t) = infla_rate;
    %unemployment rate
    disocct=(W-length(Oc(Oc>0)))/W;
    Un(t)=disocct;
    
    Y_nominal_k(t) = sum(Y_prev_k)*price_k(t);
    Y_nominal_c(t) = sum(Y_prev)*price(t);
    Y_nominal_tot(t)= Y_nominal_k(t)+Y_nominal_c(t);
    Y_real(t) = sum(Y_prev)*price(1) + sum(Y_prev_k)*price_k(1);
    Growth(t) = (Y_real(t)-Y_real(t-1))/Y_real(t-1);
    gdp_deflator(t) = Y_nominal_tot(t)/ Y_real(t);
     
    Investment(t)=I(t)*price_k(1)+sum(Y_k-inventory_k)*price_k(1)+price(1)*sum(Y); %total investment is investment plus inventory variation
 
    negcash_k=find(A_k<=0|liquidity_k<0);
    negcash=find(A<=0|liquidity<0);
    
    NetEq = liquidity-deb;
    Y_prevp=Y_prev(NetEq>0);  
    bankruptcy_rate(t) = (length(negcash_k) + length(negcash))/(F+N);

    %update bankrupted firms!
    for i=negcash                                                       %pick sequentially failed firms
        if liquidity(i)<0 && A(i)>0
           liquidity(i)=liquidity(i)+max(0,PA(W+i));
           PA(W+i)=min(0,PA(W+i));
            if liquidity(i)<0
               piB=piB+liquidity(i);
               liquidity(i)=0;
            end
           A(i)=liquidity(i)+capital_value(i)-deb(i);
        else
        bankrupt(1,t,i)=1;
        defaults(t)=defaults(t)+1;
        zzz=deb(i);                                                     %take residual debts
        if zzz>0
           
            piB=piB+(liquidity(i)-deb(i));                                               %account for bad debts
            loans=loans-zzz;
        end
        A(i)=max(0,PA(W+i))+K(i)*price_k(t+1);                                                   %initialize new firm 
        capital_value(i)=K(i)*price_k(t+1);
        PA(W+i)=min(0,PA(W+i));
        liquidity(i)=A(i)-K(i)*price_k(t+1);
        deb(i)=0;
        P(i)=mean(P);
        targetLev=0.2;
        mxY=((A(i)+targetLev*A(i)/(1-targetLev))*1/wb)*alpha;
        Y_prev(i)=min(trimmean(Y_prevp,10),mxY);
        Yd(i)=Y_prev(i);
        x(i)=Y_prev(i)/k/K(i);
        barK(i)=K(i);
        barYK(i)=Y_prev(i)/k;
        Y(i)=0;
        stock(i)=0;
        interest_r(i)=r_f;
        %%fire workers
        workforce=find(Oc==i);                                          %pick all firm i's workers
        fired=workforce;  
        if YPtest==1
        trans_EU(t)=trans_EU(t)+length(fired);
        end
        Oc(fired)=0;
        w(fired)=0;
        Leff(i)=0;
        end
    end
    
    NetEq_k = A_k;
    if isempty(find(NetEq_k>0, 1))
          warning('all capital firms are bankrupted')
%          display(seed)
       % break
        %%keep the same variable of last year        
    else
        Y_prevp_k=Y_prev_k(NetEq_k>0);  
    end
    
    initialA_k(t)=sum(PA(W+F+negcash_k))/length(negcash_k);
    for i=negcash_k                                                      %pick sequentially failed firms
        if liquidity_k(i)<0 && A_k(i)>0
           liquidity_k(i)=liquidity_k(i)+max(0,PA(W+F+i));
           PA(W+F+i)=min(0,PA(W+F+i));
            if liquidity_k(i)<0
               piB=piB+liquidity_k(i);
               liquidity_k(i)=0;
            end
           A_k(i)=liquidity_k(i)+Y_k(i)*P_k(i)-deb_k(i);
        else
        
        defaults_k(t)=defaults_k(t)+1;
        counter=counter+1;
        bankrupt(1,t,F+i)=1;
        zzz=deb_k(i);                                                     %take residual debts
        if zzz>0
           
            piB=piB+(liquidity_k(i)-deb_k(i));                                               %account for bad debts
            loans=loans-zzz;
        end
        A_k(i)=max(0,PA(W+F+i));                                                   %initialize new firm 
        PA(W+F+i)=min(0,PA(W+F+i));
        liquidity_k(i)=A_k(i);
        deb_k(i)=0;
        P_k(i)=mean(P_k);
        %maximum initial productin is given by the leverage
        targetLev=0.2;
        mxY=((A_k(i)+targetLev*A_k(i)/(1-targetLev))*1/wb)*alpha;
        
        Y_prev_k(i)=min(trimmean(Y_prevp_k,10),mxY);
        Y_kd(i)=Y_prev_k(i);
        Y_k(i)=0;
        stock_k(i)=0;
        interest_r_k(i)=r_f;
        workforce_k=find(Oc==F+i);                 %pick all firm i's workers
        fired_k=workforce_k;  
        if YPtest==1
        trans_EU(t)=trans_EU(t)+length(fired_k);
        end
        Oc(fired_k)=0;
        w(fired_k)=0;
        Leff_k(i)=0; 
        end
    end
    
    
    %% bank accounting
    piB=piB+sum(interests)+sum(interests_k) + bond_interest_rate*bonds(t-1)-interest_deposits;                                             %bank profits  

    
    if piB>0 && totE(t)>0
        Div_probb(t)=1;
        dividendsB(t)=div_B*piB;
        piB=(1-div_B)*piB;
        PA(W+1:W+F+N)=PA(W+1:W+F+N)+dividendsB(t)/(F+N);
    else 
         dividendsB(t) = 0;
         Div_probb(t)=0;
    end    
    E=E+piB;                                                            %update bank capital

    %%add bank's dividends to income of capitalists
    divstatus=dividends_income;	
    divstatus(divstatus>0)=1;	
    divstatus_k=dividends_income_k;	
    divstatus_k(divstatus_k>0)=1;	
    dividends_income_pb=dividends_income;	
    dividends_income_pbk=dividends_income_k;	
    dividends_income_k = dividends_income_k + dividendsB(t)/(F+N);	
    dividends_income = dividends_income + dividendsB(t)/(F+N);	
    divstatus_b_prev=divstatus_b;	
    if dividendsB(t)>0	
       divstatus_b=1;	
    else	
       divstatus_b=0;	
    end	
    Kdivs_p(:,t,:)=dividends_income_k;	
    Cdivs_p(:,t,:)=dividends_income;
    
    profitsB(t)=piB;

    totE(t+1)=E;
  
  
 GB(t)= TA(t) - G(t)- EXP(t)- bond_interest_rate*bonds(t-1); %JAKOB bonds(t-1) mi impalla tutto ma non potevo
 primary_GB(t) = TA(t) - G(t)- EXP(t);

 stock_bonds(t) =sum(-GB(1:t));
 bonds(t) = max(0,stock_bonds(t));
 bonds_real(t) =stock_bonds(t)/price(t);
 
 quota_exp_c(:,t)= (RIC./sum(RIC));  % quota di mercato calcolata come share entrate
 shares(:,t,:) = quota_exp_c(:,t);
 
 deficitPil(t)= -GB(t)/Y_nominal_tot(t); %% ES POLITICA FISCALE
 deficit_real(t)=-GB(t)/price(t);
 primary_deficit_pil(t) = -primary_GB(t)/Y_nominal_tot(t);

money(t) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)+E;
net_money(t) = money(t) - stock_bonds(t);

Assets(t) = sum(PA)/price(t);

if u_target - Un(t)>0
    wb = wb * (1+ wage_update_up * (u_target - Un(t))); 
else
    wb = wb * (1+ wage_update_down * (u_target - Un(t))); 
end

if YPtest==1
       if divstatus_b==1 && divstatus_b_prev==0
           trans_ND_b(t)=trans_ND_b(t)+1;
       end
       if divstatus_b==1 && divstatus_b_prev==1
           trans_DD_b(t)=trans_DD_b(t)+1;
       end
       if divstatus_b==0 && divstatus_b_prev==0
           trans_NN_b(t)=trans_NN_b(t)+1;
       end
       if divstatus_b==0 && divstatus_b_prev==1
           trans_DN_b(t)=trans_DN_b(t)+1;
       end
       
        if t<burnin
            if sum(trans_DD_b(1:t))+sum(trans_DN_b(1:t))>0
            prob_DD_b(t)=sum(trans_DD_b(1:t))/(sum(trans_DD_b(1:t))+sum(trans_DN_b(1:t)));
            prob_DN_b(t)=sum(trans_DN_b(1:t))/(sum(trans_DD_b(1:t))+sum(trans_DN_b(1:t)));
            else
            prob_DD_b(t)=prob_DD_b(t-1);
            prob_DN_b(t)=prob_DN_b(t-1);
            end
        else
            if sum(trans_DD_b((t-(horizon-1)):t))+sum(trans_DN_b((t-(horizon-1)):t))>0
            prob_DD_b(t)=sum(trans_DD_b((t-(horizon-1)):t))/(sum(trans_DD_b((t-(horizon-1)):t))+sum(trans_DN_b((t-(horizon-1)):t)));
            prob_DN_b(t)=sum(trans_DN_b((t-(horizon-1)):t))/(sum(trans_DD_b((t-(horizon-1)):t))+sum(trans_DN_b((t-(horizon-1)):t)));
            else
            prob_DD_b(t)=prob_DD_b(t-1);
            prob_DN_b(t)=prob_DN_b(t-1);
            end
        end

        if t < burnin
            if sum(trans_ND_b(1:t))+sum(trans_NN_b(1:t))>0
            prob_ND_b(t)=sum(trans_ND_b(1:t))/(sum(trans_ND_b(1:t))+sum(trans_NN_b(1:t)));
            prob_NN_b(t)=sum(trans_NN_b(1:t))/(sum(trans_ND_b(1:t))+sum(trans_NN_b(1:t)));
            else
            prob_ND_b(t)=prob_ND_b(t-1);
            prob_NN_b(t)=prob_NN_b(t-1);
            end
        else
            if sum(trans_ND_b((t-(horizon-1)):t))+sum(trans_NN_b((t-(horizon-1)):t))>0      
            prob_ND_b(t)=sum(trans_ND_b((t-(horizon-1)):t))/(sum(trans_ND_b((t-(horizon-1)):t))+sum(trans_NN_b((t-(horizon-1)):t)));
            prob_NN_b(t)=sum(trans_NN_b((t-(horizon-1)):t))/(sum(trans_ND_b((t-(horizon-1)):t))+sum(trans_NN_b((t-(horizon-1)):t)));
            else
            prob_ND_b(t)=prob_ND_b(t-1);
            prob_NN_b(t)=prob_NN_b(t-1);
            end
       end
    
    
    for i=1:F
       if divstatus(i)==1 && divstatus_prev(i)==0
           trans_ND(t)=trans_ND(t)+1;
       end
       if divstatus(i)==1 && divstatus_prev(i)==1
           trans_DD(t)=trans_DD(t)+1;
       end
       if divstatus(i)==0 && divstatus_prev(i)==0
           trans_NN(t)=trans_NN(t)+1;
       end
       if divstatus(i)==0 && divstatus_prev(i)==1
           trans_DN(t)=trans_DN(t)+1;
       end
    end
    
    
    if t<burnin
        if sum(trans_DD(1:t))+sum(trans_DN(1:t))>0
        prob_DD(t)=sum(trans_DD(1:t))/(sum(trans_DD(1:t))+sum(trans_DN(1:t)));
        prob_DN(t)=sum(trans_DN(1:t))/(sum(trans_DD(1:t))+sum(trans_DN(1:t)));
        else
        prob_DD(t)=prob_DD(t-1);
        prob_DN(t)=prob_DN(t-1);
        end
    else
        if sum(trans_DD((t-(horizon-1)):t))+sum(trans_DN((t-(horizon-1)):t))>0
        prob_DD(t)=sum(trans_DD((t-(horizon-1)):t))/(sum(trans_DD((t-(horizon-1)):t))+sum(trans_DN((t-(horizon-1)):t)));
        prob_DN(t)=sum(trans_DN((t-(horizon-1)):t))/(sum(trans_DD((t-(horizon-1)):t))+sum(trans_DN((t-(horizon-1)):t)));
        else
        prob_DD(t)=prob_DD(t-1);
        prob_DN(t)=prob_DN(t-1);
        end
    end

    if t < burnin
        if sum(trans_ND(1:t))+sum(trans_NN(1:t))>0
        prob_ND(t)=sum(trans_ND(1:t))/(sum(trans_ND(1:t))+sum(trans_NN(1:t)));
        prob_NN(t)=sum(trans_NN(1:t))/(sum(trans_ND(1:t))+sum(trans_NN(1:t)));
        else
        prob_ND(t)=prob_ND(t-1);
        prob_NN(t)=prob_NN(t-1);
        end
    else
        if sum(trans_ND((t-(horizon-1)):t))+sum(trans_NN((t-(horizon-1)):t))>0      
        prob_ND(t)=sum(trans_ND((t-(horizon-1)):t))/(sum(trans_ND((t-(horizon-1)):t))+sum(trans_NN((t-(horizon-1)):t)));
        prob_NN(t)=sum(trans_NN((t-(horizon-1)):t))/(sum(trans_ND((t-(horizon-1)):t))+sum(trans_NN((t-(horizon-1)):t)));
        else
        prob_ND(t)=prob_ND(t-1);
        prob_NN(t)=prob_NN(t-1);
        end
    end
    
    rdivs(:,t)=dividends_income_pb./price(t);	
    rdivs_pos(t)=sum(rdivs(:,t))/sum(rdivs(:,t)>0);
    rdivs_b(t)=(dividendsB(t)/(F+N))/price(t);
    
    if learn==1 && t>=1000 && t<finalshock+1
        
            shocks=find(pub_exp_cr(1:T)>0);
            lshocks=shocks+1;

            shocks=shocks(shocks>=(t-(horizon-1)) & shocks<=(t));
            lshocks=lshocks(lshocks>=(t-horizon) & lshocks<=(t));
            nshocks=(t-(horizon-1)):t;
            nshocks=setdiff(nshocks,shocks);
            nshocks=setdiff(nshocks,lshocks);

            DD1=trans_DD(shocks);
            DN1=trans_DN(shocks);
            ND1=trans_ND(shocks);
            NN1=trans_NN(shocks);

            DD2=trans_DD(lshocks);
            DN2=trans_DN(lshocks);
            ND2=trans_ND(lshocks);
            NN2=trans_NN(lshocks);

            DD3=trans_DD(nshocks);
            DN3=trans_DN(nshocks);
            ND3=trans_ND(nshocks);
            NN3=trans_NN(nshocks);
            
            
            DDb1=trans_DD_b(shocks);
            DNb1=trans_DN_b(shocks);
            NDb1=trans_ND_b(shocks);
            NNb1=trans_NN_b(shocks);

            DDb2=trans_DD_b(lshocks);
            DNb2=trans_DN_b(lshocks);
            NDb2=trans_ND_b(lshocks);
            NNb2=trans_NN_b(lshocks);

            DDb3=trans_DD_b(nshocks);
            DNb3=trans_DN_b(nshocks);
            NDb3=trans_ND_b(nshocks);
            NNb3=trans_NN_b(nshocks);
            

            prob_DD_s(t)=sum(DD1)/(sum(DD1)+sum(DN1));
            prob_DN_s(t)=sum(DN1)/(sum(DD1)+sum(DN1));
            prob_ND_s(t)=sum(ND1)/(sum(ND1)+sum(NN1));
            prob_NN_s(t)=sum(NN1)/(sum(ND1)+sum(NN1));

            prob_DD_ls(t)=sum(DD2)/(sum(DD2)+sum(DN2));
            prob_DN_ls(t)=sum(DN2)/(sum(DD2)+sum(DN2));
            prob_ND_ls(t)=sum(ND2)/(sum(ND2)+sum(NN2));
            prob_NN_ls(t)=sum(NN2)/(sum(ND2)+sum(NN2));

            prob_DD_ns(t)=sum(DD3)/(sum(DD3)+sum(DN3));
            prob_DN_ns(t)=sum(DN3)/(sum(DD3)+sum(DN3));
            prob_ND_ns(t)=sum(ND3)/(sum(ND3)+sum(NN3));
            prob_NN_ns(t)=sum(NN3)/(sum(ND3)+sum(NN3));
            
            
            prob_DD_bs(t)=sum(DDb1)/(sum(DDb1)+sum(DNb1));
            prob_DN_bs(t)=sum(DNb1)/(sum(DDb1)+sum(DNb1));
            if sum(NDb1)+sum(NNb1)>0
            prob_ND_bs(t)=sum(NDb1)/(sum(NDb1)+sum(NNb1));
            prob_NN_bs(t)=sum(NNb1)/(sum(NDb1)+sum(NNb1));
            else
            prob_ND_bs(t)=1;
            prob_NN_bs(t)=0;
            end
            
            prob_DD_bls(t)=sum(DDb2)/(sum(DDb2)+sum(DNb2));
            prob_DN_bls(t)=sum(DNb2)/(sum(DDb2)+sum(DNb2));
            if sum(NDb2)+sum(NNb2)>0
            prob_ND_bls(t)=sum(NDb2)/(sum(NDb2)+sum(NNb2));
            prob_NN_bls(t)=sum(NNb2)/(sum(NDb2)+sum(NNb2));
            else
            prob_ND_bls(t)=1;
            prob_NN_bls(t)=0;
            end
            
            
            prob_DD_bns(t)=sum(DDb3)/(sum(DDb3)+sum(DNb3));
            prob_DN_bns(t)=sum(DNb3)/(sum(DDb3)+sum(DNb3));
            if sum(NDb3)+sum(NNb3)>0
            prob_ND_bns(t)=sum(NDb3)/(sum(NDb3)+sum(NNb3));
            prob_NN_bns(t)=sum(NNb3)/(sum(NDb3)+sum(NNb3));
            else
            prob_ND_bns(t)=1;
            prob_NN_bns(t)=0;
            end

            for i=1:F
                rd_s(i)=sum(rdivs(i,shocks))/sum(rdivs(i,shocks)>0);
                rd_ls(i)=sum(rdivs(i,lshocks))/sum(rdivs(i,lshocks)>0);
                rd_ns(i)=sum(rdivs(i,nshocks))/sum(rdivs(i,nshocks)>0);
            end
            
            rdb_s(t)=sum(rdivs_b(shocks))/sum(rdivs_b(shocks)>0);
            rdb_ls(t)=sum(rdivs_b(lshocks))/sum(rdivs_b(lshocks)>0);
            rdb_ns(t)=sum(rdivs_b(nshocks))/sum(rdivs_b(nshocks)>0);
            
            md_s(t)=mean(rd_s);
            md_ls(t)=mean(rd_ls);
            md_ns(t)=mean(rd_ns);

    elseif t>=finalshock+1
            
            prob_DD_s(t)=prob_DD_s(t-1);
            prob_DN_s(t)=prob_DN_s(t-1);
            prob_ND_s(t)=prob_ND_s(t-1);
            prob_NN_s(t)=prob_NN_s(t-1);

            prob_DD_ls(t)=prob_DD_ls(t-1);
            prob_DN_ls(t)=prob_DN_ls(t-1);
            prob_ND_ls(t)=prob_ND_ls(t-1);
            prob_NN_ls(t)=prob_NN_ls(t-1);

            prob_DD_ns(t)=prob_DD_ns(t-1);
            prob_DN_ns(t)=prob_DN_ns(t-1);
            prob_ND_ns(t)=prob_ND_ns(t-1);
            prob_NN_ns(t)=prob_NN_ns(t-1);
            
            prob_DD_bs(t)=prob_DD_bs(t-1);
            prob_DN_bs(t)=prob_DN_bs(t-1);
            prob_ND_bs(t)=prob_ND_bs(t-1);
            prob_NN_bs(t)=prob_NN_bs(t-1);

            prob_DD_bls(t)=prob_DD_bls(t-1);
            prob_DN_bls(t)=prob_DN_bls(t-1);
            prob_ND_bls(t)=prob_ND_bls(t-1);
            prob_NN_bls(t)=prob_NN_bls(t-1);

            prob_DD_bns(t)=prob_DD_bns(t-1);
            prob_DN_bns(t)=prob_DN_bns(t-1);
            prob_ND_bns(t)=prob_ND_bns(t-1);
            prob_NN_bns(t)=prob_NN_bns(t-1);
            
            rdb_s(t)=rdb_s(t-1);
            rdb_ls(t)=rdb_ls(t-1);
            rdb_ns(t)=rdb_ns(t-1);
            
            md_s(t)=md_s(t-1);
            md_ls(t)=md_ls(t-1);
            md_ns(t)=md_ns(t-1);

    end
    
    if t < burnin
        if sum(rdivs_b(2:t)>0)>0
                rd_b=sum(rdivs_b(2:t))/sum(rdivs_b(2:t)>0);
        else
                rd_b=0;
        end
        for i=1:F
            if sum(rdivs(i,2:t)>0)>0
                rd(i)=sum(rdivs(i,2:t))/sum(rdivs(i,2:t)>0);
            else
                rd(i)=0;
            end
        end
    else
        for i=1:F
            rd(i)=sum(rdivs(i,(t-(horizon-1)):t))/sum(rdivs(i,(t-(horizon-1)):t)>0);
        end
            rd_b=sum(rdivs_b((t-(horizon-1)):t))/sum(rdivs_b((t-(horizon-1)):t)>0);
    end
    
        if t<burnin+1 && t>3
        depend=rdivs_b(2:t);
        depend=depend(depend>0);
        depend=depend-rd_b;
        ind=depend(1:length(depend)-2);
        depend=depend(2:length(depend)-1);
        end
        if t>=burnin+1
        depend=rdivs_b(t-horizon:t);
        depend=depend(depend>0);
        depend=depend-rd_b;
        ind=depend(1:length(depend)-1);
        depend=depend(2:length(depend));
        end
        if t>10
        coeffdb=inv(ind*ind')*ind*depend';
        else
        coeffdb=0.75;
        end
        coeffdb=min(coeffdb,0.999);
        
        if t<burnin
            pay1b=rdivs_b(t)*ones(2000,1);
            pay2b=zeros(2000,1);
        else
        pay1b=rd_b*ones(2000,1);
        pay1b(1)=depend(length(depend));
        q=1;
        while abs(pay1b(q)-rd_b)>0.0005
           q=q+1;
           pay1b(q)=rd_b+coeffdb*(pay1b(q-1)-rd_b);
        end
        pay2b=zeros(length(pay1b),1);
        end
        if learn_act && t>=1500
        shockperiods=gperiods-t+1;
        lshockperiods=shockperiods+1;
        shockperiods=shockperiods(shockperiods>0);
        lshockperiods=lshockperiods(lshockperiods>0);
        probs11=prob_DD_bns(t)*ones(2000,1);
        probs12=prob_DN_bns(t)*ones(2000,1);
        probs21=prob_ND_bns(t)*ones(2000,1);
        probs22=prob_NN_bns(t)*ones(2000,1);
        probs11(shockperiods)=prob_DD_bs(t);
        probs12(shockperiods)=prob_DN_bs(t);
        probs21(shockperiods)=prob_ND_bs(t);
        probs22(shockperiods)=prob_NN_bs(t);
        probs11(lshockperiods)=prob_DD_bls(t);
        probs12(lshockperiods)=prob_DN_bls(t);
        probs21(lshockperiods)=prob_ND_bls(t);
        probs22(lshockperiods)=prob_NN_bls(t);
        pay1b(:)=pay1b(:)-rd_b+rdb_ns(t);
        pay1b(shockperiods)=pay1b(shockperiods)-rdb_ns(t)+rdb_s(t);
        pay1b(lshockperiods)=pay1b(lshockperiods)-rdb_ns(t)+rdb_ls(t);
        pay1b(1,:)=rdivs_b(t);
        permdb=projectIncome(probs11,probs12,probs21,probs22,pay1b,pay2b,r_d);
        else
        permdb=projectIncome(prob_DD_b(t),prob_DN_b(t),prob_ND_b(t),prob_NN_b(t),pay1b,pay2b,r_d);
        end
    
        md(t)=mean(rd);
        if t<burnin+1 && t>3
        depend=rdivs_pos(2:t)-md(t);
        ind=depend(1:t-2);
        depend=depend(2:t-1);
        end
        if t>=burnin+1
        depend=rdivs_pos(t-horizon:t)-md(t);
        ind=depend(1:horizon-1);
        depend=depend(2:horizon);
        end
        if t>10
        for i=1:F
        coeffsd(i)=inv(ind*ind')*ind*depend';
        coeffsd(i)=min(coeffsd(i),0.999);
        end
        else
        coeffsd(:)=0.75;
        end
        if t<burnin
            pay1=ones(2000,F);
            for i=1:F
            pay1(:,i)=rdivs(i,t);
            end
            pay2=zeros(2000,F);
        else
        pay1=md(t)*ones(2000,F);
        currentpay=rdivs(:,t);
        currentpay(currentpay==0)=rdivs_pos(t);
        pay1(1,:)=currentpay;
        q=1;
        while abs(mean(pay1(q,:))-md(t))>0.0005
           q=q+1;
           pay1(q,:)=md(t)+coeffsd'.*(pay1(q-1,:)-md(t));
        end
        pay2=zeros(length(pay1),F);
        end
        pay1(:,:)=pay1(:,:);
        pay2(:,:)=pay2(:,:);
        if learn_act && t>=1500
        shockperiods=gperiods-t+1;
        lshockperiods=shockperiods+1;
        shockperiods=shockperiods(shockperiods>0);
        lshockperiods=lshockperiods(lshockperiods>0);
        probs11=prob_DD_ns(t)*ones(2000,1);
        probs12=prob_DN_ns(t)*ones(2000,1);
        probs21=prob_ND_ns(t)*ones(2000,1);
        probs22=prob_NN_ns(t)*ones(2000,1);
        probs11(shockperiods)=prob_DD_s(t);
        probs12(shockperiods)=prob_DN_s(t);
        probs21(shockperiods)=prob_ND_s(t);
        probs22(shockperiods)=prob_NN_s(t);
        probs11(lshockperiods)=prob_DD_ls(t);
        probs12(lshockperiods)=prob_DN_ls(t);
        probs21(lshockperiods)=prob_ND_ls(t);
        probs22(lshockperiods)=prob_NN_ls(t);
        pay1(:,:)=pay1(:,:)-md(t)+md_ns(t);
        pay1(shockperiods,:)=pay1(shockperiods,:)-md_ns(t)+md_s(t);
        pay1(lshockperiods,:)=pay1(lshockperiods,:)-md_ns(t)+md_ls(t);
        pay1(1,:)=rdivs(:,t);
        permds=projectIncomeMult(probs11,probs12,probs21,probs22,pay1,pay2,r_d);
        else
        permds=projectIncomeMult(prob_DD(t),prob_DN(t),prob_ND(t),prob_NN(t),pay1,pay2,r_d);
        end
        for i=1:F
            permd=permds(:,i);
            V_div=permd;
            if dividends_income(i)>0
                if dividendsB(t)>0
                yd(i)=(PA(W+i)/price(t)+V_div(1)+permdb(1))/(1/r_d);
                else
                yd(i)=(PA(W+i)/price(t)+V_div(1)+permdb(2))/(1/r_d);
                end
            else
                if dividendsB(t)>0
                yd(i)=(PA(W+i)/price(t)+V_div(2)+permdb(1))/(1/r_d);
                else
                yd(i)=(PA(W+i)/price(t)+V_div(2)+permdb(2))/(1/r_d);
                end
            end
            YP_d(t)=YP_d(t)+V_div(1);
            YP_nd(t)=YP_nd(t)+V_div(2);
        end
        YP_d(t)=YP_d(t)/F;
        YP_nd(t)=YP_nd(t)/F;
        
       if t==currentshock-1 && t>=1500 && learn_act
           divi=sum(divstatus==1);
           ndivi=sum(divstatus==0);
           diviN=round(divi*prob_DN_s(t));
           ndiviD=round(ndivi*prob_ND_s(t));
           divips=divi-diviN+ndiviD;
           ndivips=ndivi-ndiviD+diviN;

           divipsN=round(divips*prob_DN_ls(t));
           ndivipsD=round(ndivips*prob_ND_ls(t));

           divipls=divips-divipsN+ndivipsD;
           ndivipls=ndivips-ndivipsD+divipsN;

           div_s=mean(pay1(2,:));
           div_ls=mean(pay1(3,:));
           counter=counter+1;
           DPS(qshocks)=divips;
           NPS(qshocks)=ndivips;
           DPLS(qshocks)=divipls;
           NPLS(qshocks)=ndivipls;
           DIPS(qshocks)=div_s;
           DIPLS(qshocks)=div_ls;
           
           
           divib=sum(divstatus_b==1);
           ndivib=sum(divstatus_b==0);
           divibN=round(divib*prob_DN_bs(t));
           ndivibD=round(ndivib*prob_ND_bs(t));
           divibps=divib-divibN+ndivibD;
           ndivibps=ndivib-ndivibD+divibN;

           divibpsN=round(divibps*prob_DN_bls(t));
           ndivibpsD=round(ndivibps*prob_ND_bls(t));

           divibpls=divibps-divibpsN+ndivibpsD;
           ndivibpls=ndivibps-ndivibpsD+divibpsN;
           
           div_bs=pay1b(2);
           div_bls=pay1b(3);
           DPBS(qshocks)=divibps;
           NPBS(qshocks)=ndivibps;
           DPBLS(qshocks)=divibpls;
           NPBLS(qshocks)=ndivibpls;
           DIBPS(qshocks)=div_bs;
           DIBPLS(qshocks)=div_bls;
       end
        
    for l=1:N
       if divstatus_k(l)==1 && divstatus_prev_k(l)==0
           trans_ND_k(t)=trans_ND_k(t)+1;
       end
       if divstatus_k(l)==1 && divstatus_prev_k(l)==1
           trans_DD_k(t)=trans_DD_k(t)+1;
       end
       if divstatus_k(l)==0 && divstatus_prev_k(l)==0
           trans_NN_k(t)=trans_NN_k(t)+1;
       end
       if divstatus_k(l)==0 && divstatus_prev_k(l)==1
           trans_DN_k(t)=trans_DN_k(t)+1;
       end
    end
    

    if t<=horizon
        if sum(trans_DD_k(1:t))+sum(trans_DN_k(1:t))>0
        prob_DD_k(t)=sum(trans_DD_k(1:t))/(sum(trans_DD_k(1:t))+sum(trans_DN_k(1:t)));
        prob_DN_k(t)=sum(trans_DN_k(1:t))/(sum(trans_DD_k(1:t))+sum(trans_DN_k(1:t)));
        else
        prob_DD_k(t)=prob_DD_k(t-1);
        prob_DN_k(t)=prob_DN_k(t-1);
        end
    else
        if sum(trans_DD_k((t-(horizon-1)):t))+sum(trans_DN_k((t-(horizon-1)):t))>0
        prob_DD_k(t)=sum(trans_DD_k((t-(horizon-1)):t))/(sum(trans_DD_k((t-(horizon-1)):t))+sum(trans_DN_k((t-(horizon-1)):t)));
        prob_DN_k(t)=sum(trans_DN_k((t-(horizon-1)):t))/(sum(trans_DD_k((t-(horizon-1)):t))+sum(trans_DN_k((t-(horizon-1)):t)));
        else
        prob_DD_k(t)=prob_DD_k(t-1);
        prob_DN_k(t)=prob_DN_k(t-1);
        end
    end

    if t<burnin
        if sum(trans_ND_k(1:t))+sum(trans_NN_k(1:t))>0
        prob_ND_k(t)=sum(trans_ND_k(1:t))/(sum(trans_ND_k(1:t))+sum(trans_NN_k(1:t)));
        prob_NN_k(t)=sum(trans_NN_k(1:t))/(sum(trans_ND_k(1:t))+sum(trans_NN_k(1:t)));
        else
        prob_ND_k(t)=prob_ND_k(t-1);
        prob_NN_k(t)=prob_NN_k(t-1);
        end
    else
        if sum(trans_ND_k((t-(horizon-1)):t))+sum(trans_NN_k((t-(horizon-1)):t))>0
        prob_ND_k(t)=sum(trans_ND_k((t-(horizon-1)):t))/(sum(trans_ND_k((t-(horizon-1)):t))+sum(trans_NN_k((t-(horizon-1)):t)));
        prob_NN_k(t)=sum(trans_NN_k((t-(horizon-1)):t))/(sum(trans_ND_k((t-(horizon-1)):t))+sum(trans_NN_k((t-(horizon-1)):t)));
        else
        prob_ND_k(t)=prob_ND_k(t-1);
        prob_NN_k(t)=prob_NN_k(t-1);
        end
    end
    
    
    rdivs_k(:,t)=dividends_income_pbk./price(t);
    rdivs_kpos(t)=sum(rdivs_k(:,t))/sum(rdivs_k(:,t)>0);
    
    if learn==1 && t>=1000 && t<finalshock+1
            shocks=find(pub_exp_cr(1:T)>0);
            lshocks=shocks+1;

            shocks=shocks(shocks>=(t-(horizon-1)) & shocks<=(t));
            lshocks=lshocks(lshocks>=(t-horizon) & lshocks<=(t));
            nshocks=(t-(horizon-1)):t;
            nshocks=setdiff(nshocks,shocks);
            nshocks=setdiff(nshocks,lshocks);

            DDk1=trans_DD_k(shocks);
            DNk1=trans_DN_k(shocks);
            NDk1=trans_ND_k(shocks);
            NNk1=trans_NN_k(shocks);

            DDk2=trans_DD_k(lshocks);
            DNk2=trans_DN_k(lshocks);
            NDk2=trans_ND_k(lshocks);
            NNk2=trans_NN_k(lshocks);

            DDk3=trans_DD_k(nshocks);
            DNk3=trans_DN_k(nshocks);
            NDk3=trans_ND_k(nshocks);
            NNk3=trans_NN_k(nshocks);

            prob_DD_ks(t)=sum(DDk1)/(sum(DDk1)+sum(DNk1));
            prob_DN_ks(t)=sum(DNk1)/(sum(DDk1)+sum(DNk1));
            prob_ND_ks(t)=sum(NDk1)/(sum(NDk1)+sum(NNk1));
            prob_NN_ks(t)=sum(NNk1)/(sum(NDk1)+sum(NNk1));

            prob_DD_kls(t)=sum(DDk2)/(sum(DDk2)+sum(DNk2));
            prob_DN_kls(t)=sum(DNk2)/(sum(DDk2)+sum(DNk2));
            prob_ND_kls(t)=sum(NDk2)/(sum(NDk2)+sum(NNk2));
            prob_NN_kls(t)=sum(NNk2)/(sum(NDk2)+sum(NNk2));

            prob_DD_kns(t)=sum(DDk3)/(sum(DDk3)+sum(DNk3));
            prob_DN_kns(t)=sum(DNk3)/(sum(DDk3)+sum(DNk3));
            prob_ND_kns(t)=sum(NDk3)/(sum(NDk3)+sum(NNk3));
            prob_NN_kns(t)=sum(NNk3)/(sum(NDk3)+sum(NNk3));

            for i=1:N
                rdk_s(i)=sum(rdivs_k(i,shocks))/sum(rdivs_k(i,shocks)>0);
                rdk_ls(i)=sum(rdivs_k(i,lshocks))/sum(rdivs_k(i,lshocks)>0);
                rdk_ns(i)=sum(rdivs_k(i,nshocks))/sum(rdivs_k(i,nshocks)>0);
            end

            mdk_s(t)=mean(rdk_s);
            mdk_ls(t)=mean(rdk_ls);
            mdk_ns(t)=mean(rdk_ns);

    elseif t>=finalshock+1
            
            prob_DD_ks(t)=prob_DD_ks(t-1);
            prob_DN_ks(t)=prob_DN_ks(t-1);
            prob_ND_ks(t)=prob_ND_ks(t-1);
            prob_NN_ks(t)=prob_NN_ks(t-1);

            prob_DD_kls(t)=prob_DD_kls(t-1);
            prob_DN_kls(t)=prob_DN_kls(t-1);
            prob_ND_kls(t)=prob_ND_kls(t-1);
            prob_NN_kls(t)=prob_NN_kls(t-1);

            prob_DD_kns(t)=prob_DD_kns(t-1);
            prob_DN_kns(t)=prob_DN_kns(t-1);
            prob_ND_kns(t)=prob_ND_kns(t-1);
            prob_NN_kns(t)=prob_NN_kns(t-1);

            mdk_s(t)=mdk_s(t-1);
            mdk_ls(t)=mdk_ls(t-1);
            mdk_ns(t)=mdk_ns(t-1);
    end
    
      
    if t < burnin
        for l=1:N
            if sum(rdivs_k(l,1:t)>0)>0
                rdk(l)=sum(rdivs_k(l,1:t))/sum(rdivs_k(l,1:t)>0);
            else
                rdk(l)=0;
            end
        end
    else
       for l=1:N
            rdk(l)=sum(rdivs_k(l,(t-(horizon-1)):t))/sum(rdivs_k(l,(t-(horizon-1)):t)>0);
       end 
        
    end
    
        mdk(t)=mean(rdk);
        if t<burnin+1 && t>3
        depend=rdivs_kpos(2:t)-mdk(t);
        ind=depend(1:t-2);
        depend=depend(2:t-1);
        end
        if t>=burnin+1
        depend=rdivs_kpos(t-horizon:t)-mdk(t);
        ind=depend(1:horizon-1);
        depend=depend(2:horizon);
        end
        if t>10
        for i=1:N
        coeffsdk(i)=inv(ind*ind')*ind*depend';
        coeffsdk(i)=min(coeffsdk(i),0.999);
        end
        else
        coeffsdk(:)=0.75;
        end
        if t<burnin
            pay1=ones(2000,N);
            for i=1:N
            pay1(:,i)=rdivs_k(i,t);
            end
            pay2=zeros(2000,N);
        else
        pay1=mdk(t)*ones(2000,N);
        currentpay=rdivs_k(:,t);
        currentpay(currentpay==0)=rdivs_kpos(t);
        pay1(1,:)=currentpay;
        q=1;
        while abs(mean(pay1(q,:))-mdk(t))>0.0005
           q=q+1;
           pay1(q,:)=mdk(t)+coeffsdk'.*(pay1(q-1,:)-mdk(t));
        end
        pay2=zeros(length(pay1),N);
        end
        pay1(:,:)=pay1(:,:);
        pay2(:,:)=pay2(:,:);
        if learn_act && t>=1500
        shockperiods=gperiods-t+1;
        lshockperiods=shockperiods+1;
        shockperiods=shockperiods(shockperiods>0);
        lshockperiods=lshockperiods(lshockperiods>0);
        probs11=prob_DD_kns(t)*ones(2000,1);
        probs12=prob_DN_kns(t)*ones(2000,1);
        probs21=prob_ND_kns(t)*ones(2000,1);
        probs22=prob_NN_kns(t)*ones(2000,1);
        probs11(shockperiods)=prob_DD_ks(t);
        probs12(shockperiods)=prob_DN_ks(t);
        probs21(shockperiods)=prob_ND_ks(t);
        probs22(shockperiods)=prob_NN_ks(t);
        probs11(lshockperiods)=prob_DD_kls(t);
        probs12(lshockperiods)=prob_DN_kls(t);
        probs21(lshockperiods)=prob_ND_kls(t);
        probs22(lshockperiods)=prob_NN_kls(t);
        pay1(:,:)=pay1(:,:)-mdk(t)+mdk_ns(t);
        pay1(shockperiods,:)=pay1(shockperiods,:)-mdk_ns(t)+mdk_s(t);
        pay1(lshockperiods,:)=pay1(lshockperiods,:)-mdk_ns(t)+mdk_ls(t);
        pay1(1,:)=rdivs_k(:,t);
        permdks=projectIncomeMult(probs11,probs12,probs21,probs22,pay1,pay2,r_d);
        else
        permdks=projectIncomeMult(prob_DD_k(t),prob_DN_k(t),prob_ND_k(t),prob_NN_k(t),pay1,pay2,r_d);
        end
        for l=1:N
            permdk=permdks(:,l);
            V_divk=permdk;
            if dividends_income_k(l)>0
                if dividendsB(t)>0
                ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdb(1))/(1/r_d);
                else
                ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdb(2))/(1/r_d);
                end
            else
                if dividendsB(t)>0
                ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdb(1))/(1/r_d);
                else
                ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdb(2))/(1/r_d);
                end
            end
            YP_dk(t)=YP_dk(t)+V_divk(1);
            YP_ndk(t)=YP_ndk(t)+V_divk(2);
        end
        YP_dk(t)=YP_dk(t)/N;
        YP_ndk(t)=YP_ndk(t)/N;
        
        if t==currentshock-1 && t>=1500 && learn_act
           divik=sum(divstatus_k==1);
           ndivik=sum(divstatus_k==0);
           divikN=round(divik*prob_DN_ks(t));
           ndivikD=round(ndivik*prob_ND_ks(t));
           divikps=divik-divikN+ndivikD;
           ndivikps=ndivik-ndivikD+divikN;

           divikpsN=round(divikps*prob_DN_kls(t));
           ndivikpsD=round(ndivikps*prob_ND_kls(t));

           divikpls=divikps-divikpsN+ndivikpsD;
           ndivikpls=ndivikps-ndivikpsD+divikpsN;

           divk_s=mean(pay1(2,:));
           divk_ls=mean(pay1(3,:));
           
           DPKS(qshocks)=divikps;
           NPKS(qshocks)=ndivikps;
           DPKLS(qshocks)=divikpls;
           NPKLS(qshocks)=ndivikpls;
           DIKPS(qshocks)=divk_s;
           DIKPLS(qshocks)=divk_ls;
           
        end
       
     if t==currentshock-1 && t>=1500 && learn_act
       empl=sum(Oc>0);
       unempl=sum(Oc==0);
       emplU=round(empl*prob_EU_s(t));
       unemplE=round(unempl*prob_UE_s(t));
       emplps=empl-emplU+unemplE;
       unemplps=unempl-unemplE+emplU;
       
       emplpsU=round(emplps*prob_EU_ls(t));
       unemplpsE=round(unemplps*prob_UE_ls(t));
       
       emplpls=emplps-emplpsU+unemplpsE;
       unemplpls=unemplps-unemplpsE+emplpsU;
       
       wage_s=pay1w(2);
       subsidy_s=pay2w(2);
       wage_ls=pay1w(3);
       subsidy_ls=pay2w(3);
       
       EPS(qshocks)=emplps;
       UPS(qshocks)=unemplps;
       EPLS(qshocks)=emplpls;
       UPLS(qshocks)=unemplpls;
       WPS(qshocks)=wage_s;
       SPS(qshocks)=subsidy_s;
       WPLS(qshocks)=wage_ls;
       SPLS(qshocks)=subsidy_ls;
       
     end
    
     if t==currentshock && t>1500 && learn_act
       actualemps=sum(Oc>0);
       actualunemps=sum(Oc==0);

       actualwages=rwage(t);
       actualsubsidys=rwage(t)*unemployment_subsidy;
       
       EAS(qshocks)=actualemps;
       UAS(qshocks)=actualunemps;
       WAS(qshocks)=actualwages;
       SAS(qshocks)=actualsubsidys;

    end

    if t==currentshock+1 && t>1501 && learn_act
       actualempls=sum(Oc>0);
       actualunempls=sum(Oc==0);

       actualwagels=rwage(t);
       actualsubsidyls=rwage(t)*unemployment_subsidy;
       
       EALS(qshocks)=actualempls;
       UALS(qshocks)=actualunempls;
       WALS(qshocks)=actualwagels;
       SALS(qshocks)=actualsubsidyls;

    end
end

if t==currentshock && t>1500 && learn_act
   actualdivs=sum(divstatus==1);
   actualndivs=sum(divstatus==0);
   
   divpaids=rdivs(:,t);
   divpaids=mean(divpaids(divpaids>0));
   
   actualdivks=sum(divstatus_k==1);
   actualndivks=sum(divstatus_k==0);
   
   divpaidks=rdivs_k(:,t);
   divpaidks=mean(divpaidks(divpaidks>0));
   
   if divstatus_b==1
   actualdivbs=1;
   actualndivbs=0;
   else
   actualdivbs=1;
   actualndivbs=0;
   end
   
   divpaidbs=rdivs_b(t);
   
   DAS(qshocks)=actualdivs;
   NAS(qshocks)=actualndivs;
   DIAS(qshocks)=divpaids;
   
   DAKS(qshocks)=actualdivks;
   NAKS(qshocks)=actualndivks;
   DIKAS(qshocks)=divpaidks;
   
   DABS(qshocks)=actualdivbs;
   NABS(qshocks)=actualndivbs;
   DIBAS(qshocks)=divpaidbs;
   
end

if t==currentshock+1 && t>1501 && learn_act
   actualdivls=sum(divstatus==1);
   actualndivls=sum(divstatus==0);
   
   divpaidls=rdivs(:,t);
   divpaidls=mean(divpaidls(divpaidls>0));
   
   actualdivkls=sum(divstatus_k==1);
   actualndivkls=sum(divstatus_k==0);
   
   divpaidkls=rdivs_k(:,t);
   divpaidkls=mean(divpaidkls(divpaidkls>0));
    
   if divstatus_b==1
   actualdivbls=1;
   actualndivbls=0;
   else
   actualdivbls=1;
   actualndivbls=0;
   end
   
   divpaidbls=rdivs_b(t);
   
   DALS(qshocks)=actualdivls;
   NALS(qshocks)=actualndivls;
   DIALS(qshocks)=divpaidls;
   
   DAKLS(qshocks)=actualdivkls;
   NAKLS(qshocks)=actualndivkls;
   DIKALS(qshocks)=divpaidkls;
   
   DABLS(qshocks)=actualdivbls;
   NABLS(qshocks)=actualndivbls;
   DIBALS(qshocks)=divpaidbls;
   
end

wages_t(t) = wb;

 end
 
 for i=1:(W+F+N)
   cactual(i,:)=cactual(i,:);
   cdemand1(i,:)=cdemand1(i,:);
   cdemand2(i,:)=cdemand2(i,:);
   welfare_c(i)=mean((cactual(i,1500:2500).^(1-tau)-1)/(1-tau));
   welfare_dc1(i)=mean((cdemand1(i,1500:2500).^(1-tau)-1)/(1-tau));
   welfare_dc2(i)=mean((cdemand2(i,1500:2500).^(1-tau)-1)/(1-tau));
 end

predictions_ws=zeros(4,3);
actual_ws=zeros(4,3);
predictions_wls=zeros(4,3);
actual_wls=zeros(4,3);

predictions_fs=zeros(3,3);
actual_fs=zeros(3,3);
predictions_fls=zeros(3,3);
actual_fls=zeros(3,3);

predictions_ks=zeros(3,3);
actual_ks=zeros(3,3);
predictions_kls=zeros(3,3);
actual_kls=zeros(3,3);

predictions_bs=zeros(3,3);
actual_bs=zeros(3,3);
predictions_bls=zeros(3,3);
actual_bls=zeros(3,3);

if endshock==1
   EPS=rmmissing(EPS);
   predictions_ws(1,1)=mean(EPS);
   SEM = std(EPS)/sqrt(length(EPS));               
   ts = tinv([0.025  0.975],length(EPS)-1);      
   CI = mean(EPS) + ts*SEM;
   predictions_ws(1,2)=CI(1);
   predictions_ws(1,3)=CI(2);
   
   
   UPS=rmmissing(UPS);
   predictions_ws(2,1)=mean(UPS);
   SEM = std(UPS)/sqrt(length(UPS));               
   ts = tinv([0.025  0.975],length(UPS)-1);      
   CI = mean(UPS) + ts*SEM;
   predictions_ws(2,2)=CI(1);
   predictions_ws(2,3)=CI(2);
   
   
   WPS=rmmissing(WPS);
   predictions_ws(3,1)=mean(WPS);
   SEM = std(WPS)/sqrt(length(WPS));               
   ts = tinv([0.025  0.975],length(WPS)-1);      
   CI = mean(WPS) + ts*SEM;
   predictions_ws(3,2)=CI(1);
   predictions_ws(3,3)=CI(2);
   
   
   SPS=rmmissing(SPS);
   predictions_ws(4,1)=mean(SPS);
   SEM = std(SPS)/sqrt(length(SPS));               
   ts = tinv([0.025  0.975],length(SPS)-1);      
   CI = mean(SPS) + ts*SEM;
   predictions_ws(4,2)=CI(1);
   predictions_ws(4,3)=CI(2);
   
   
   EAS=rmmissing(EAS);
   actual_ws(1,1)=mean(EAS);
   SEM = std(EAS)/sqrt(length(EAS));               
   ts = tinv([0.025  0.975],length(EAS)-1);      
   CI = mean(EAS) + ts*SEM;
   actual_ws(1,2)=CI(1);
   actual_ws(1,3)=CI(2);
   
   
   UAS=rmmissing(UAS);
   actual_ws(2,1)=mean(UAS);
   SEM = std(UAS)/sqrt(length(UAS));               
   ts = tinv([0.025  0.975],length(UAS)-1);      
   CI = mean(UAS) + ts*SEM;
   actual_ws(2,2)=CI(1);
   actual_ws(2,3)=CI(2);
   
   
   WAS=rmmissing(WAS);
   actual_ws(3,1)=mean(WAS);
   SEM = std(WAS)/sqrt(length(WAS));               
   ts = tinv([0.025  0.975],length(WAS)-1);      
   CI = mean(WAS) + ts*SEM;
   actual_ws(3,2)=CI(1);
   actual_ws(3,3)=CI(2);
   
   
   SAS=rmmissing(SAS);
   actual_ws(4,1)=mean(SAS);
   SEM = std(SAS)/sqrt(length(SAS));               
   ts = tinv([0.025  0.975],length(SAS)-1);      
   CI = mean(SAS) + ts*SEM;
   actual_ws(4,2)=CI(1);
   actual_ws(4,3)=CI(2);
   
   
   EPLS=rmmissing(EPLS);
   predictions_wls(1,1)=mean(EPLS);
   SEM = std(EPLS)/sqrt(length(EPLS));               
   ts = tinv([0.025  0.975],length(EPLS)-1);      
   CI = mean(EPLS) + ts*SEM;
   predictions_wls(1,2)=CI(1);
   predictions_wls(1,3)=CI(2);
   
   
   UPLS=rmmissing(UPLS);
   predictions_wls(2,1)=mean(UPLS);
   SEM = std(UPLS)/sqrt(length(UPLS));               
   ts = tinv([0.025  0.975],length(UPLS)-1);      
   CI = mean(UPLS) + ts*SEM;
   predictions_wls(2,2)=CI(1);
   predictions_wls(2,3)=CI(2);
   
   
   WPLS=rmmissing(WPLS);
   predictions_wls(3,1)=mean(WPLS);
   SEM = std(WPLS)/sqrt(length(WPLS));               
   ts = tinv([0.025  0.975],length(WPLS)-1);      
   CI = mean(WPLS) + ts*SEM;
   predictions_wls(3,2)=CI(1);
   predictions_wls(3,3)=CI(2);
   
   
   SPLS=rmmissing(SPLS);
   predictions_wls(4,1)=mean(SPLS);
   SEM = std(SPLS)/sqrt(length(SPLS));               
   ts = tinv([0.025  0.975],length(SPLS)-1);      
   CI = mean(SPLS) + ts*SEM;
   predictions_wls(4,2)=CI(1);
   predictions_wls(4,3)=CI(2);
   
   
   EALS=rmmissing(EALS);
   actual_wls(1,1)=mean(EALS);
   SEM = std(EALS)/sqrt(length(EALS));               
   ts = tinv([0.025  0.975],length(EALS)-1);      
   CI = mean(EALS) + ts*SEM;
   actual_wls(1,2)=CI(1);
   actual_wls(1,3)=CI(2);
   
   
   UALS=rmmissing(UALS);
   actual_wls(2,1)=mean(UALS);
   SEM = std(UALS)/sqrt(length(UALS));               
   ts = tinv([0.025  0.975],length(UALS)-1);      
   CI = mean(UALS) + ts*SEM;
   actual_wls(2,2)=CI(1);
   actual_wls(2,3)=CI(2);
   
   
   WALS=rmmissing(WALS);
   actual_wls(3,1)=mean(WALS);
   SEM = std(WALS)/sqrt(length(WALS));               
   ts = tinv([0.025  0.975],length(WALS)-1);      
   CI = mean(WALS) + ts*SEM;
   actual_wls(3,2)=CI(1);
   actual_wls(3,3)=CI(2);
   
   
   SALS=rmmissing(SALS);
   actual_wls(4,1)=mean(SALS);
   SEM = std(SALS)/sqrt(length(SALS));               
   ts = tinv([0.025  0.975],length(SALS)-1);      
   CI = mean(SALS) + ts*SEM;
   actual_wls(4,2)=CI(1);
   actual_wls(4,3)=CI(2);
   
   
   DPS=rmmissing(DPS);
   predictions_fs(1,1)=mean(DPS);
   SEM = std(DPS)/sqrt(length(DPS));               
   ts = tinv([0.025  0.975],length(DPS)-1);      
   CI = mean(DPS) + ts*SEM;
   predictions_fs(1,2)=CI(1);
   predictions_fs(1,3)=CI(2);
   
   NPS=rmmissing(NPS);
   predictions_fs(2,1)=mean(NPS);
   SEM = std(NPS)/sqrt(length(NPS));               
   ts = tinv([0.025  0.975],length(NPS)-1);      
   CI = mean(NPS) + ts*SEM;
   predictions_fs(2,2)=CI(1);
   predictions_fs(2,3)=CI(2);

   
   DIPS=rmmissing(DIPS);
   predictions_fs(3,1)=mean(DIPS);
   SEM = std(DIPS)/sqrt(length(DIPS));               
   ts = tinv([0.025  0.975],length(DIPS)-1);      
   CI = mean(DIPS) + ts*SEM;
   predictions_fs(3,2)=CI(1);
   predictions_fs(3,3)=CI(2);
   
   
   DAS=rmmissing(DAS);
   actual_fs(1,1)=mean(DAS);
   SEM = std(DAS)/sqrt(length(DAS));               
   ts = tinv([0.025  0.975],length(DAS)-1);      
   CI = mean(DAS) + ts*SEM;
   actual_fs(1,2)=CI(1);
   actual_fs(1,3)=CI(2);
   
   NAS=rmmissing(NAS);
   actual_fs(2,1)=mean(NAS);
   SEM = std(NAS)/sqrt(length(NAS));               
   ts = tinv([0.025  0.975],length(NAS)-1);      
   CI = mean(NAS) + ts*SEM;
   actual_fs(2,2)=CI(1);
   actual_fs(2,3)=CI(2);

   
   DIAS=rmmissing(DIAS);
   actual_fs(3,1)=mean(DIAS);
   SEM = std(DIAS)/sqrt(length(DIAS));               
   ts = tinv([0.025  0.975],length(DIAS)-1);      
   CI = mean(DIAS) + ts*SEM;
   actual_fs(3,2)=CI(1);
   actual_fs(3,3)=CI(2);
   
   
   DPLS=rmmissing(DPLS);
   predictions_fls(1,1)=mean(DPLS);
   SEM = std(DPLS)/sqrt(length(DPLS));               
   ts = tinv([0.025  0.975],length(DPLS)-1);      
   CI = mean(DPLS) + ts*SEM;
   predictions_fls(1,2)=CI(1);
   predictions_fls(1,3)=CI(2);
   
   NPLS=rmmissing(NPLS);
   predictions_fls(2,1)=mean(NPLS);
   SEM = std(NPLS)/sqrt(length(NPLS));               
   ts = tinv([0.025  0.975],length(NPLS)-1);      
   CI = mean(NPLS) + ts*SEM;
   predictions_fls(2,2)=CI(1);
   predictions_fls(2,3)=CI(2);

   
   DIPLS=rmmissing(DIPLS);
   predictions_fls(3,1)=mean(DIPLS);
   SEM = std(DIPLS)/sqrt(length(DIPLS));               
   ts = tinv([0.025  0.975],length(DIPLS)-1);      
   CI = mean(DIPLS) + ts*SEM;
   predictions_fls(3,2)=CI(1);
   predictions_fls(3,3)=CI(2);
   
   
   DALS=rmmissing(DALS);
   actual_fls(1,1)=mean(DALS);
   SEM = std(DALS)/sqrt(length(DALS));               
   ts = tinv([0.025  0.975],length(DALS)-1);      
   CI = mean(DALS) + ts*SEM;
   actual_fls(1,2)=CI(1);
   actual_fls(1,3)=CI(2);
   
   NALS=rmmissing(NALS);
   actual_fls(2,1)=mean(NALS);
   SEM = std(NALS)/sqrt(length(NALS));               
   ts = tinv([0.025  0.975],length(NALS)-1);      
   CI = mean(NALS) + ts*SEM;
   actual_fls(2,2)=CI(1);
   actual_fls(2,3)=CI(2);

   
   DIALS=rmmissing(DIALS);
   actual_fls(3,1)=mean(DIALS);
   SEM = std(DIALS)/sqrt(length(DIALS));               
   ts = tinv([0.025  0.975],length(DIALS)-1);      
   CI = mean(DIALS) + ts*SEM;
   actual_fls(3,2)=CI(1);
   actual_fls(3,3)=CI(2);
    
   
   
   
   DPKS=rmmissing(DPKS);
   predictions_ks(1,1)=mean(DPKS);
   SEM = std(DPKS)/sqrt(length(DPKS));               
   ts = tinv([0.025  0.975],length(DPKS)-1);      
   CI = mean(DPKS) + ts*SEM;
   predictions_ks(1,2)=CI(1);
   predictions_ks(1,3)=CI(2);
   
   NPKS=rmmissing(NPKS);
   predictions_ks(2,1)=mean(NPKS);
   SEM = std(NPKS)/sqrt(length(NPKS));               
   ts = tinv([0.025  0.975],length(NPKS)-1);      
   CI = mean(NPKS) + ts*SEM;
   predictions_ks(2,2)=CI(1);
   predictions_ks(2,3)=CI(2);

   
   DIKPS=rmmissing(DIKPS);
   predictions_ks(3,1)=mean(DIKPS);
   SEM = std(DIKPS)/sqrt(length(DIKPS));               
   ts = tinv([0.025  0.975],length(DIKPS)-1);      
   CI = mean(DIKPS) + ts*SEM;
   predictions_ks(3,2)=CI(1);
   predictions_ks(3,3)=CI(2);
   
   
   DAKS=rmmissing(DAKS);
   actual_ks(1,1)=mean(DAKS);
   SEM = std(DAKS)/sqrt(length(DAKS));               
   ts = tinv([0.025  0.975],length(DAKS)-1);      
   CI = mean(DAKS) + ts*SEM;
   actual_ks(1,2)=CI(1);
   actual_ks(1,3)=CI(2);
   
   NAKS=rmmissing(NAKS);
   actual_ks(2,1)=mean(NAKS);
   SEM = std(NAKS)/sqrt(length(NAKS));               
   ts = tinv([0.025  0.975],length(NAKS)-1);      
   CI = mean(NAKS) + ts*SEM;
   actual_ks(2,2)=CI(1);
   actual_ks(2,3)=CI(2);

   
   DIKAS=rmmissing(DIKAS);
   actual_ks(3,1)=mean(DIKAS);
   SEM = std(DIKAS)/sqrt(length(DIKAS));               
   ts = tinv([0.025  0.975],length(DIKAS)-1);      
   CI = mean(DIKAS) + ts*SEM;
   actual_ks(3,2)=CI(1);
   actual_ks(3,3)=CI(2);
   
   
   DPKLS=rmmissing(DPKLS);
   predictions_kls(1,1)=mean(DPKLS);
   SEM = std(DPKLS)/sqrt(length(DPKLS));               
   ts = tinv([0.025  0.975],length(DPKLS)-1);      
   CI = mean(DPKLS) + ts*SEM;
   predictions_kls(1,2)=CI(1);
   predictions_kls(1,3)=CI(2);
   
   NPKLS=rmmissing(NPKLS);
   predictions_kls(2,1)=mean(NPKLS);
   SEM = std(NPKLS)/sqrt(length(NPKLS));               
   ts = tinv([0.025  0.975],length(NPKLS)-1);      
   CI = mean(NPKLS) + ts*SEM;
   predictions_kls(2,2)=CI(1);
   predictions_kls(2,3)=CI(2);

   
   DIKPLS=rmmissing(DIKPLS);
   predictions_kls(3,1)=mean(DIKPLS);
   SEM = std(DIKPLS)/sqrt(length(DIKPLS));               
   ts = tinv([0.025  0.975],length(DIKPLS)-1);      
   CI = mean(DIKPLS) + ts*SEM;
   predictions_kls(3,2)=CI(1);
   predictions_kls(3,3)=CI(2);
   
   
   DAKLS=rmmissing(DAKLS);
   actual_kls(1,1)=mean(DAKLS);
   SEM = std(DAKLS)/sqrt(length(DAKLS));               
   ts = tinv([0.025  0.975],length(DAKLS)-1);      
   CI = mean(DAKLS) + ts*SEM;
   actual_kls(1,2)=CI(1);
   actual_kls(1,3)=CI(2);
   
   NAKLS=rmmissing(NAKLS);
   actual_kls(2,1)=mean(NAKLS);
   SEM = std(NAKLS)/sqrt(length(NAKLS));               
   ts = tinv([0.025  0.975],length(NAKLS)-1);      
   CI = mean(NAKLS) + ts*SEM;
   actual_kls(2,2)=CI(1);
   actual_kls(2,3)=CI(2);

   
   DIKALS=rmmissing(DIKALS);
   actual_kls(3,1)=mean(DIKALS);
   SEM = std(DIKALS)/sqrt(length(DIKALS));               
   ts = tinv([0.025  0.975],length(DIKALS)-1);      
   CI = mean(DIKALS) + ts*SEM;
   actual_kls(3,2)=CI(1);
   actual_kls(3,3)=CI(2);
   
   
   DPBS=rmmissing(DPBS);
   predictions_bs(1,1)=mean(DPBS);
   SEM = std(DPBS)/sqrt(length(DPBS));               
   ts = tinv([0.025  0.975],length(DPBS)-1);      
   CI = mean(DPBS) + ts*SEM;
   predictions_bs(1,2)=CI(1);
   predictions_bs(1,3)=CI(2);
   
   NPBS=rmmissing(NPBS);
   predictions_bs(2,1)=mean(NPBS);
   SEM = std(NPBS)/sqrt(length(NPBS));               
   ts = tinv([0.025  0.975],length(NPBS)-1);      
   CI = mean(NPBS) + ts*SEM;
   predictions_bs(2,2)=CI(1);
   predictions_bs(2,3)=CI(2);

   
   DIBPS=rmmissing(DIBPS);
   predictions_bs(3,1)=mean(DIBPS);
   SEM = std(DIBPS)/sqrt(length(DIBPS));               
   ts = tinv([0.025  0.975],length(DIBPS)-1);      
   CI = mean(DIBPS) + ts*SEM;
   predictions_bs(3,2)=CI(1);
   predictions_bs(3,3)=CI(2);
   
   
   DABS=rmmissing(DABS);
   actual_bs(1,1)=mean(DABS);
   SEM = std(DABS)/sqrt(length(DABS));               
   ts = tinv([0.025  0.975],length(DABS)-1);      
   CI = mean(DABS) + ts*SEM;
   actual_bs(1,2)=CI(1);
   actual_bs(1,3)=CI(2);
   
   NABS=rmmissing(NABS);
   actual_bs(2,1)=mean(NABS);
   SEM = std(NABS)/sqrt(length(NABS));               
   ts = tinv([0.025  0.975],length(NABS)-1);      
   CI = mean(NABS) + ts*SEM;
   actual_bs(2,2)=CI(1);
   actual_bs(2,3)=CI(2);

   
   DIBAS=rmmissing(DIBAS);
   actual_bs(3,1)=mean(DIBAS);
   SEM = std(DIBAS)/sqrt(length(DIBAS));               
   ts = tinv([0.025  0.975],length(DIBAS)-1);      
   CI = mean(DIBAS) + ts*SEM;
   actual_bs(3,2)=CI(1);
   actual_bs(3,3)=CI(2);
   
   
   DPBLS=rmmissing(DPBLS);
   predictions_bls(1,1)=mean(DPBLS);
   SEM = std(DPBLS)/sqrt(length(DPBLS));               
   ts = tinv([0.025  0.975],length(DPBLS)-1);      
   CI = mean(DPBLS) + ts*SEM;
   predictions_bls(1,2)=CI(1);
   predictions_bls(1,3)=CI(2);
   
   NPBLS=rmmissing(NPBLS);
   predictions_bls(2,1)=mean(NPBLS);
   SEM = std(NPBLS)/sqrt(length(NPBLS));               
   ts = tinv([0.025  0.975],length(NPBLS)-1);      
   CI = mean(NPBLS) + ts*SEM;
   predictions_bls(2,2)=CI(1);
   predictions_bls(2,3)=CI(2);

   
   DIBPLS=rmmissing(DIBPLS);
   predictions_bls(3,1)=mean(DIBPLS);
   SEM = std(DIBPLS)/sqrt(length(DIBPLS));               
   ts = tinv([0.025  0.975],length(DIBPLS)-1);      
   CI = mean(DIBPLS) + ts*SEM;
   predictions_bls(3,2)=CI(1);
   predictions_bls(3,3)=CI(2);
   
   
   DABLS=rmmissing(DABLS);
   actual_bls(1,1)=mean(DABLS);
   SEM = std(DABLS)/sqrt(length(DABLS));               
   ts = tinv([0.025  0.975],length(DABLS)-1);      
   CI = mean(DABLS) + ts*SEM;
   actual_bls(1,2)=CI(1);
   actual_bls(1,3)=CI(2);
   
   NABLS=rmmissing(NABLS);
   actual_bls(2,1)=mean(NABLS);
   SEM = std(NABLS)/sqrt(length(NABLS));               
   ts = tinv([0.025  0.975],length(NABLS)-1);      
   CI = mean(NABLS) + ts*SEM;
   actual_bls(2,2)=CI(1);
   actual_bls(2,3)=CI(2);

   
   DIBALS=rmmissing(DIBALS);
   actual_bls(3,1)=mean(DIBALS);
   SEM = std(DIBALS)/sqrt(length(DIBALS));               
   ts = tinv([0.025  0.975],length(DIBALS)-1);      
   CI = mean(DIBALS) + ts*SEM;
   actual_bls(3,2)=CI(1);
   actual_bls(3,3)=CI(2);
   
   
end

et = toc;

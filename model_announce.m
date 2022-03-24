
function  [Y_real,noprod,cactual,cdemand1,cdemand2,welfare_cs,welfare_dc1s,welfare_dc2s,welfare_c,welfare_dc1,welfare_dc2,Ygap,rDIV,rDIVK,dincome,dincomek,DIV,DIVK,OCC,wincome,deviation,opt_share,neut_share,pes_share,types_agg,gperiod,Growth,Dek,Dek2,Kdem2,bankruptcy_rate,residualC,YP_e,YP_u,prob_EE,prob_EU,prob_UU,prob_UE,prob_EE_p,prob_EU_p,prob_UU_p,prob_UE_p,prob_DD,prob_DN,prob_NN,prob_ND,prob_DD_k,prob_DN_k,prob_NN_k,prob_ND_k,price,Div_prob,Div_probb,Div_probk,profitsB,Cdivs_p,Kdivs_p,shares,dividendsB,EXPcontrol,Invent,Assets,baryk,Kdes,Kdem,Dexp,Kdivs,Cdivs,valI,unsatisfiedDemandK,actualEXP,X1,X2,totK, gdp_deflator, Investment,I, consumption, Prod_k, Prod_c, Un, totalDeb, totalDeb_k,stock_bonds,GB,TA,G,wages_t, DC,totE,average_interest_rate,impulse_rw_opt,YPE_true,YPE_anchor,YPU_true,YPU_anchor,rwage,Occ_status,Type,et] = model_announce(seed, T, par,Hor)
tic 

Bk=1;                               %no. of banks
F=100; %200;%                         %no. of firms/capitalists
W=1000; % 3000;%                        %no. of workers
N = 20;% 50;%                          %no. of capital producing firms
z_c=par(1);                              %no. of aplications in consumption good market
z_k = par(2);                         %no. of aplications in capital good market
z_e = par(3);                            %number of job applications  

totdiv=zeros(1,F);
totdivk=zeros(1,N);

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

Kshare=0.5;
Nshare=0;
types=zeros(1,W+F+N);
randtypes=rand([1,W+F+N]);
for i=1:(W+F+N)
    if randtypes(i)<=Kshare
        types(i)=1;
    else
        types(i)=-1;
    end
end
randtypes2=rand([1,W+F+N]);
for i=1:(W+F+N)
    if types(i)==-1 && randtypes2(i) < Nshare
        types(i)=0;
    end
end

randswitch_w=rand(T,W);
randswitch_f=rand(T,F);
randswitch_k=rand(T,N);

OCC=zeros(1,T);
wincome=zeros(1,T);
dincome=zeros(1,T);
dincomek=zeros(1,T);
DIV=zeros(1,T);
DIVK=zeros(1,T);
rDIV=zeros(1,T);
rDIVK=zeros(1,T);
YPtest=1;
horizon=Hor;
tau=2/3;
noprod=0;
burnin=horizon+1;
r_f=0.01;%   0.015;% 0.005          %general refinancing rate
r_d=r_f/2;                          %deposit rate
counter=0;
Type=zeros(T,W+F+N);
adapt=0.95;
opt_share=zeros(1,T);
pes_share=zeros(1,T);
neut_share=zeros(1,T);
opt_share(1)=sum(types==1)/(W+F+N);
neut_share(1)=sum(types==0)/(W+F+N);
pes_share(1)=sum(types==-1)/(W+F+N);
Type(1,:)=types;
types_agg=zeros(1,T);
trans_UE=zeros(1,T);
trans_UE(1)=6000;
trans_EU=zeros(1,T);
trans_EU(1)=500;
trans_UU=zeros(1,T);
trans_UU(1)=4000;
trans_EE=zeros(1,T);
trans_EE(1)=9500;
trans_UE_p=zeros(1,T);
trans_EU_p=zeros(1,T);
trans_UU_p=zeros(1,T);
trans_EE_p=zeros(1,T);
trans_ND_p=zeros(1,T);
trans_DN_p=zeros(1,T);
trans_DD_p=zeros(1,T);
trans_NN_p=zeros(1,T);

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

trans_ND_bp=zeros(1,T);
trans_DN_bp=zeros(1,T);
trans_DD_bp=zeros(1,T);
trans_NN_bp=zeros(1,T);
trans_ND_kp=zeros(1,T);
trans_DN_kp=zeros(1,T);
trans_DD_kp=zeros(1,T);
trans_NN_kp=zeros(1,T);
prob_EE=zeros(1,T);
prob_EU=zeros(1,T);
prob_UE=zeros(1,T);
prob_UU=zeros(1,T);
prob_EE_opt=zeros(1,T);
prob_EU_opt=zeros(1,T);
prob_UE_opt=zeros(1,T);
prob_UU_opt=zeros(1,T);
prob_EE_pes=zeros(1,T);
prob_EU_pes=zeros(1,T);
prob_UE_pes=zeros(1,T);
prob_UU_pes=zeros(1,T);
prob_EE_p=zeros(1,T);
prob_EU_p=zeros(1,T);
prob_UE_p=zeros(1,T);
prob_UU_p=zeros(1,T);
prob_DD=zeros(1,T);
prob_DN=zeros(1,T);
prob_ND=zeros(1,T);
prob_NN=zeros(1,T);
prob_DD_p=zeros(1,T);
prob_DN_p=zeros(1,T);
prob_ND_p=zeros(1,T);
prob_NN_p=zeros(1,T);
prob_DD_opt=zeros(1,T);
prob_DN_opt=zeros(1,T);
prob_ND_opt=zeros(1,T);
prob_NN_opt=zeros(1,T);
prob_DD_pes=zeros(1,T);
prob_DN_pes=zeros(1,T);
prob_ND_pes=zeros(1,T);
prob_NN_pes=zeros(1,T);
prob_DD_k=zeros(1,T);
prob_DN_k=zeros(1,T);
prob_ND_k=zeros(1,T);
prob_NN_k=zeros(1,T);
prob_DD_kp=zeros(1,T);
prob_DN_kp=zeros(1,T);
prob_ND_kp=zeros(1,T);
prob_NN_kp=zeros(1,T);
prob_DD_k_opt=zeros(1,T);
prob_DN_k_opt=zeros(1,T);
prob_ND_k_opt=zeros(1,T);
prob_NN_k_opt=zeros(1,T);
prob_DD_k_pes=zeros(1,T);
prob_DN_k_pes=zeros(1,T);
prob_ND_k_pes=zeros(1,T);
prob_NN_k_pes=zeros(1,T);
prob_DD_b=zeros(1,T);
prob_DN_b=zeros(1,T);
prob_ND_b=zeros(1,T);
prob_NN_b=zeros(1,T);
prob_DD_bp=zeros(1,T);
prob_DN_bp=zeros(1,T);
prob_ND_bp=zeros(1,T);
prob_NN_bp=zeros(1,T);
prob_DD_b_opt=zeros(1,T);
prob_DN_b_opt=zeros(1,T);
prob_ND_b_opt=zeros(1,T);
prob_NN_b_opt=zeros(1,T);
prob_DD_b_pes=zeros(1,T);
prob_DN_b_pes=zeros(1,T);
prob_ND_b_pes=zeros(1,T);
prob_NN_b_pes=zeros(1,T);
rd_e=zeros(1,F);
rdk_e=zeros(1,N);
deviation=zeros(W+F+N,T);
YPE_anchor=zeros(1,T);
YPU_anchor=zeros(1,T);
YPE_true=zeros(1,T);
YPU_true=zeros(1,T);
YPE_belief=zeros(W,T);
yd_belief=zeros(F,T);
ydk_belief=zeros(N,T);
coeffw=rand(1,1);
coeffsd=rand(F,1);
md=zeros(1,T);
coeffsdk=rand(N,1);
mdk=zeros(1,T);
step=1;
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
rd_true=zeros(1,F);
rdk=zeros(1,N);
rdk_true=zeros(1,N);
yd=zeros(1,F);
ydk=zeros(1,N);
residualC=zeros(1,T);
bonds_real=zeros(1,T);
deficit_real=zeros(1,T);
Ygap=zeros(1,T);

%set by the sensitivity
xi = par(4);                          %memory parameter human wealth
chi = par(5);                            %fraction of wealth devoted to consumption
q_adj = par(6);                            %quantity adjustment parameter
p_adj = par(7);                            %price adjustment parameter    
mu =  par(8);                           %bank's gross mark-up
eta = par(9);                        %capital depreciation
Iprob=par(10);                         %probability of investing
phi =  par(11);                        %bank's leverage parameter
theta=par(12);                         %rate of debt reimbursment
delta = par(13);                        %memory parameter in the capital utilization rate
alpha = par(14);                             %labour productivity
k = par(15);                            %capital productivity
div =par(16);                            %share of dividends
div_B=0.3;
barX=par(17);                         %desired capital utilization
inventory_depreciation = par(18);           %rate at which capital firms' inventories depreciate
b1 = par(19);   
b2 = par(20);                        %Parameters for risk evaluation by banks
b_k1 = par(21);
b_k2 = par(22);

interest_rate = par(23);
subsidy = par(24);
tax_rate = par(27);

wage_update_up = par(28);
wage_update_down = par(29);
u_target = par(30);
%Government
Gov = 0;

%phillips curve
wb=1.5;                                % initial wage rate


%else government will use unemployment subsidies as line below.
bond_interest_rate = interest_rate;   %% ricordare di aggiungere anche questi con liquidità
unemployment_subsidy_init = subsidy;

G=zeros(1,T);                       %government expenditures
TA=zeros(1,T);                      %government income
GB=zeros(1,T);                      %governament budget  GB = TA - G - EXP -->???
GB(1)=-6440;
EXP = Gov*ones (1,T);       % spesa pubblica EROGABILE, update erogata in fondo
actualEXP=EXP;
EXP(1,1)=0;              % inizializzazione, vedi sotto
Gshock=1;
EXPcontrol = zeros(1,T);
%maxperiods=csvread('maxY.csv');
%minperiods=csvread('minY.csv');
shockperiod = randi([1000,1100]);
announce = shockperiod-10;
%shockperiod = maxperiods(seed);
opinions=1;
switching=0;
belief_shock=0.1;
duration=20;
impulse_rw_opt=ones(1,T+1);
impulse_rw_pes=ones(1,T+1);
impulse_EU_opt=ones(1,T);
impulse_EU_pes=ones(1,T);
impulse_UU_opt=ones(1,T);
impulse_UU_pes=ones(1,T);
impulse_rw_opt(shockperiod)=1+belief_shock;
impulse_rw_pes(shockperiod)=1-belief_shock;
impulse_EU_opt(shockperiod)=1-belief_shock;
impulse_EU_pes(shockperiod)=1+belief_shock;
impulse_UU_opt(shockperiod)=1-belief_shock;
impulse_UU_pes(shockperiod)=1+belief_shock;

impulse_rdb_opt=ones(1,T+1);
impulse_rdb_pes=ones(1,T+1);
impulse_DN_b_opt=ones(1,T);
impulse_DN_b_pes=ones(1,T);
impulse_NN_b_opt=ones(1,T);
impulse_NN_b_pes=ones(1,T);
impulse_DN_b_opt(shockperiod)=1-belief_shock;
impulse_DN_b_pes(shockperiod)=1+belief_shock;
impulse_NN_b_opt(shockperiod)=1-belief_shock;
impulse_NN_b_pes(shockperiod)=1+belief_shock;
impulse_rdb_opt(shockperiod)=1+belief_shock;
impulse_rdb_pes(shockperiod)=1-belief_shock;

impulse_rd=ones(F,T+1);
impulse_DN_opt=ones(1,T);
impulse_DN_pes=ones(1,T);
impulse_NN_opt=ones(1,T);
impulse_NN_pes=ones(1,T);
impulse_DN_opt(shockperiod)=1-belief_shock;
impulse_DN_pes(shockperiod)=1+belief_shock;
impulse_NN_opt(shockperiod)=1-belief_shock;
impulse_NN_pes(shockperiod)=1+belief_shock;

for i=1:F
    if types(W+i)==1
        impulse_rd(i,shockperiod)=1+belief_shock;
    end
    if types(W+i)==-1
        impulse_rd(i,shockperiod)=1-belief_shock;
    end
end

impulse_rdk=ones(N,T+1);
impulse_DN_k_opt=ones(1,T);
impulse_DN_k_pes=ones(1,T);
impulse_NN_k_opt=ones(1,T);
impulse_NN_k_pes=ones(1,T);
impulse_DN_k_opt(shockperiod)=1-belief_shock;
impulse_DN_k_pes(shockperiod)=1+belief_shock;
impulse_NN_k_opt(shockperiod)=1-belief_shock;
impulse_NN_k_pes(shockperiod)=1+belief_shock;

for i=1:N
    if types(W+F+i)==1
        impulse_rdk(i,shockperiod)=1+belief_shock;
    end
    if types(W+F+i)==-1
        impulse_rdk(i,shockperiod)=1-belief_shock;
    end
end

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
YPE(:,1)=rw(1);
YPU(:,1)=rw(1)*subsidy;
Kdivs = NaN(1,T,N);
Cdivs = NaN(1,T,F);
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
cdemand1=zeros(W+F+N,T);
cdemand2=zeros(W+F+N,T);
cactual=zeros(W+F+N,T);
welfare_dc1=zeros(W+F+N,1);
welfare_dc2=zeros(W+F+N,1);
welfare_c=zeros(W+F+N,1);
welfare_dc1s=zeros(W+F+N,1);
welfare_dc2s=zeros(W+F+N,1);
welfare_cs=zeros(W+F+N,1);
Growth=NaN(1,T);
credit_mismatch=NaN(1,T);
Occ_status_prev=NaN(T,W);
divstatus=zeros(1,F);
divstatus_k=zeros(1,N);
divstatus_b=0;
%%%END OF MODEL INITIALIZATION

% cons_budget = ones(1,F+N+W)*1;
dividends_income_k = zeros(1,N);
dividends_income   = zeros(1,F);
permanent_income = ones(F+W+N,1)./price(1);
%%%HEREAFTER THE MODEL STARTS RUNNING FOR T PERIODS

money(1:2) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)+totE(1);  %What is this supposed to be? Why is bank's equity part of "money"? Why is debt subtracted???
%trigger=0;
gperiod=0;
%%%%%%%%%%%%%%%%%%
 for t=2:T
     disp(t)
     if t==shockperiod && Gshock==1	
         %EXP(1,t) = (0.1*Y_real(t-1)/price(1));	
         EXP(1,t)=0.1*alpha*1000;	
         gperiod=shockperiod;	
     end
     
     if A<=0                                                             %if all the firms are defaulted, then stop simulation 
       % break
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
                if Gshock==1 && t==shockperiod
                    De(i) = Y_prev(i) + (-stock(i))*q_adj + gexp(i);
                else
                    De(i) = Y_prev(i) + (-stock(i))*q_adj;
                end
            elseif stock(i)<=0 && P(i)<price(t)
                P(i)=P(i)*(1+shock_p(t,i)*p_adj);
                if Gshock==1 && t==shockperiod
                    De(i)=Y_prev(i);
                else
                    De(i)=Y_prev(i);
                end
            elseif stock(i)>0 &&P(i)>price(t)
                P(i)=P(i)*(1-shock_p(t,i)*p_adj);
                if Gshock==1 && t==shockperiod
                    De(i)= Y_prev(i);
                else
                    De(i)= Y_prev(i);
                end
            elseif stock(i)>0 && P(i)<=price(t)
                if Gshock==1 && t==shockperiod
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
    
    if t<shockperiod && opinions==1
        prob_UU_p(t)=prob_UU(t);
        prob_UE_p(t)=prob_UE(t);
        prob_EE_p(t)=prob_EE(t);
        prob_EU_p(t)=prob_EU(t);
    end
        
    if t>=shockperiod && opinions==1
        for j=1:W
        if Occ_status(1,t,j)==1 && Occ_status_prev(t,j)==0
            trans_UE_p(t)=trans_UE_p(t)+1;
        end
        if Occ_status(1,t,j)==0 && Occ_status_prev(t,j)==1
            trans_EU_p(t)=trans_EU_p(t)+1;
        end
        if Occ_status(1,t,j)==0 && Occ_status_prev(t,j)==0
            trans_UU_p(t)=trans_UU_p(t)+1;
        end
        if Occ_status(1,t,j)==1 && Occ_status_prev(t,j)==1
            trans_EE_p(t)=trans_EE_p(t)+1;
        end
        end

        if sum(trans_EE_p(shockperiod:t))+sum(trans_EU_p(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_EE_p(t)=sum(trans_EE_p(shockperiod:t))/(sum(trans_EE_p(shockperiod:t))+sum(trans_EU_p(shockperiod:t)));
            prob_EU_p(t)=sum(trans_EU_p(shockperiod:t))/(sum(trans_EE_p(shockperiod:t))+sum(trans_EU_p(shockperiod:t)));
            else
            prob_EE_p(t)=sum(trans_EE_p((t-(horizon-1)):t))/(sum(trans_EE_p((t-(horizon-1)):t))+sum(trans_EU_p((t-(horizon-1)):t)));
            prob_EU_p(t)=sum(trans_EU_p((t-(horizon-1)):t))/(sum(trans_EE_p((t-(horizon-1)):t))+sum(trans_EU_p((t-(horizon-1)):t)));
            end
        else
            prob_EE_p(t)=0;
            prob_EU_p(t)=0;
        end

        if sum(trans_UE_p(shockperiod:t))+sum(trans_UU_p(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_UU_p(t)=sum(trans_UU_p(shockperiod:t))/(sum(trans_UE_p(shockperiod:t))+sum(trans_UU_p(shockperiod:t)));
            prob_UE_p(t)=sum(trans_UE_p(shockperiod:t))/(sum(trans_UE_p(shockperiod:t))+sum(trans_UU_p(shockperiod:t)));
            else
            prob_UU_p(t)=sum(trans_UU_p((t-(horizon-1)):t))/(sum(trans_UE_p((t-(horizon-1)):t))+sum(trans_UU_p((t-(horizon-1)):t)));
            prob_UE_p(t)=sum(trans_UE_p((t-(horizon-1)):t))/(sum(trans_UE_p((t-(horizon-1)):t))+sum(trans_UU_p((t-(horizon-1)):t)));  
            end
        else
            prob_UU_p(t)=0;
            prob_UE_p(t)=0;
        end
    end
    
    if t>shockperiod && opinions==1
       tMatE_true=[prob_EE_p(t) prob_EU_p(t); prob_UE_p(t) prob_UU_p(t)];
       if isinf(cond(tMatE_true))
            impulse_EU_opt(t)=impulse_EU_opt(t-1);
            impulse_UU_opt(t)=impulse_UU_opt(t-1);
            impulse_EU_pes(t)=impulse_EU_pes(t-1);
            impulse_UU_pes(t)=impulse_UU_pes(t-1);
       else
            impulse_EU_opt(t)=adapt*impulse_EU_opt(t-1)+(1-adapt)*(prob_EU_p(t)/prob_EU_opt(t-1));
            impulse_UU_opt(t)=adapt*impulse_UU_opt(t-1)+(1-adapt)*(prob_UU_p(t)/prob_UU_opt(t-1));
            impulse_EU_pes(t)=adapt*impulse_EU_pes(t-1)+(1-adapt)*(prob_EU_p(t)/prob_EU_pes(t-1));
            impulse_UU_pes(t)=adapt*impulse_UU_pes(t-1)+(1-adapt)*(prob_UU_p(t)/prob_UU_pes(t-1));
       end
       if impulse_EU_opt(t)>1 || impulse_EU_opt(t-1)==1
          impulse_EU_opt(t)=1;
       end
       if impulse_UU_opt(t)>1 || impulse_UU_opt(t-1)==1
          impulse_UU_opt(t)=1;
       end
       if impulse_EU_pes(t)<1 || impulse_EU_pes(t-1)==1
          impulse_EU_pes(t)=1;
       end
       if impulse_UU_pes(t)<1 || impulse_UU_pes(t-1)==1
          impulse_UU_pes(t)=1;
       end
    end
    
    prob_EU_opt(t)=max(0,min(1,prob_EU(t)*impulse_EU_opt(t)));
    prob_EU_pes(t)=max(0,min(1,prob_EU(t)*impulse_EU_pes(t)));
    prob_UU_opt(t)=max(0,min(1,prob_UU(t)*impulse_UU_opt(t)));
    prob_UU_pes(t)=max(0,min(1,prob_UU(t)*impulse_UU_pes(t)));
    prob_EE_opt(t)=1-prob_EU_opt(t);
    prob_EE_pes(t)=1-prob_EU_pes(t);
    prob_UE_opt(t)=1-prob_UU_opt(t);
    prob_UE_pes(t)=1-prob_UU_pes(t);
    
    rwage(t)=(wb*(1-tax_rate))/price(t);
    
    if t < burnin
        rw(t)=sum(rwage(2:t))/t;
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
    cw(t)=coeffw;
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
    permw=projectIncome(prob_EE(t),prob_EU(t),prob_UE(t),prob_UU(t),pay1,pay2,r_d);
    
    YPE(:,t)=(PA(1:W)/price(t)+permw(1))/(1/r_d);
    YPU(:,t)=(PA(1:W)/price(t)+permw(2))/(1/r_d);
    YPE_belief(:,t)=(permw(1))/(1/r_d);
    
    
    if opinions==1
        if t<shockperiod && t>=announce
            pay1_opt=pay1;
            pay2_opt=pay2;
            pay1_opt(shockperiod-t:shockperiod+duration-t)=pay1_opt(shockperiod-t:shockperiod+duration-t)*impulse_rw_opt(shockperiod);
            pay2_opt(shockperiod-t:shockperiod+duration-t)=pay2_opt(shockperiod-t:shockperiod+duration-t)*impulse_rw_opt(shockperiod);
            probs11=zeros(1,length(pay1_opt));
            probs12=zeros(1,length(pay1_opt));
            probs21=zeros(1,length(pay1_opt));
            probs22=zeros(1,length(pay1_opt));
            probs11(:)=prob_EE(t);
            probs12(:)=prob_EU(t);
            probs21(:)=prob_UE(t);
            probs22(:)=prob_UU(t);
            probs12(shockperiod-t:shockperiod+duration-t)=prob_EU(t)*impulse_EU_opt(shockperiod);
            probs11(shockperiod-t:shockperiod+duration-t)=1-probs12(shockperiod-t:shockperiod+duration-t);
            probs22(shockperiod-t:shockperiod+duration-t)=prob_UU(t)*impulse_UU_opt(shockperiod);
            probs21(shockperiod-t:shockperiod+duration-t)=1-probs22(shockperiod-t:shockperiod+duration-t);
            permwopt=projectIncome(probs11,probs12,probs21,probs22,pay1_opt,pay2_opt,r_d);
            pay1_pes=pay1;
            pay2_pes=pay2;
            pay1_pes(shockperiod-t:shockperiod+duration-t)=pay1_pes(shockperiod-t:shockperiod+duration-t)*impulse_rw_pes(shockperiod);
            pay2_pes(shockperiod-t:shockperiod+duration-t)=pay2_pes(shockperiod-t:shockperiod+duration-t)*impulse_rw_pes(shockperiod);
            probs11=zeros(1,length(pay1_pes));
            probs12=zeros(1,length(pay1_pes));
            probs21=zeros(1,length(pay1_pes));
            probs22=zeros(1,length(pay1_pes));
            probs11(:)=prob_EE(t);
            probs12(:)=prob_EU(t);
            probs21(:)=prob_UE(t);
            probs22(:)=prob_UU(t);
            probs12(shockperiod-t:shockperiod+duration-t)=prob_EU(t)*impulse_EU_pes(shockperiod);
            probs11(shockperiod-t:shockperiod+duration-t)=1-probs12(shockperiod-t:shockperiod+duration-t);
            probs22(shockperiod-t:shockperiod+duration-t)=prob_UU(t)*impulse_UU_pes(shockperiod);
            probs21(shockperiod-t:shockperiod+duration-t)=1-probs22(shockperiod-t:shockperiod+duration-t);
            permwpes=projectIncome(probs11,probs12,probs21,probs22,pay1_pes,pay2_pes,r_d);
            
            for i=1:W
            if types(i)==1
            YPE(i,t)=(PA(i)/price(t)+permwopt(1))/(1/r_d);
            YPU(i,t)=(PA(i)/price(t)+permwopt(2))/(1/r_d);
            YPE_belief(i,t)=(permwopt(1))/(1/r_d);
            end
            if types(i)==-1
            YPE(i,t)=(PA(i)/price(t)+permwpes(1))/(1/r_d);
            YPU(i,t)=(PA(i)/price(t)+permwpes(2))/(1/r_d);
            YPE_belief(i,t)=(permwpes(1))/(1/r_d);
            end
            end
        end
        
        if t>=shockperiod && t< shockperiod+duration
            pay1_opt=pay1;
            pay2_opt=pay2;
            pay1_opt(1:((shockperiod+duration)-t))=pay1_opt(1:((shockperiod+duration)-t))*impulse_rw_opt(t);
            pay2_opt(1:((shockperiod+duration)-t))=pay2_opt(1:((shockperiod+duration)-t))*impulse_rw_opt(t);
            probs11=zeros(1,length(pay1_opt));
            probs12=zeros(1,length(pay1_opt));
            probs21=zeros(1,length(pay1_opt));
            probs22=zeros(1,length(pay1_opt));
            probs11(:)=prob_EE(t);
            probs12(:)=prob_EU(t);
            probs21(:)=prob_UE(t);
            probs22(:)=prob_UU(t);
            probs11(1:((shockperiod+duration)-t))=prob_EE_opt(t);
            probs12(1:((shockperiod+duration)-t))=prob_EU_opt(t);
            probs21(1:((shockperiod+duration)-t))=prob_UE_opt(t);
            probs22(1:((shockperiod+duration)-t))=prob_UU_opt(t);
            permwopt=projectIncome(probs11,probs12,probs21,probs22,pay1_opt,pay2_opt,r_d);
            pay1_pes=pay1;
            pay2_pes=pay2;
            pay1_pes(1:((shockperiod+duration)-t))=pay1_pes(1:((shockperiod+duration)-t))*impulse_rw_pes(t);
            pay2_pes(1:((shockperiod+duration)-t))=pay2_pes(1:((shockperiod+duration)-t))*impulse_rw_pes(t);
            probs11=zeros(1,length(pay1_pes));
            probs12=zeros(1,length(pay1_pes));
            probs21=zeros(1,length(pay1_pes));
            probs22=zeros(1,length(pay1_pes));
            probs11(:)=prob_EE(t);
            probs12(:)=prob_EU(t);
            probs21(:)=prob_UE(t);
            probs22(:)=prob_UU(t);
            probs11(1:((shockperiod+duration)-t))=prob_EE_pes(t);
            probs12(1:((shockperiod+duration)-t))=prob_EU_pes(t);
            probs21(1:((shockperiod+duration)-t))=prob_UE_pes(t);
            probs22(1:((shockperiod+duration)-t))=prob_UU_pes(t);
            permwpes=projectIncome(probs11,probs12,probs21,probs22,pay1_pes,pay2_pes,r_d);
            
            for i=1:W
            if types(i)==1
            YPE(i,t)=(PA(i)/price(t)+permwopt(1))/(1/r_d);
            YPU(i,t)=(PA(i)/price(t)+permwopt(2))/(1/r_d);
            YPE_belief(i,t)=(permwopt(1))/(1/r_d);
            end
            if types(i)==-1
            YPE(i,t)=(PA(i)/price(t)+permwpes(1))/(1/r_d);
            YPU(i,t)=(PA(i)/price(t)+permwpes(2))/(1/r_d);
            YPE_belief(i,t)=(permwpes(1))/(1/r_d);
            end
            end
        end
    end
    
    YP_e(t)=sum(YPE(:,t));
    YP_u(t)=sum(YPU(:,t));
    
    if t>=shockperiod && opinions==1 && switching==1
       if t-shockperiod < burnin
        rw_true=sum(rwage(shockperiod:t))/(t-(shockperiod-1));
       else
        rw_true=sum(rwage(t-(horizon-1):t))/horizon;
       end
       pay1_true=rw_true*ones(2000,1);
       pay2_true=rw_true*unemployment_subsidy*ones(2000,1);
       pay1_true(1)=rwage(t);
       pay2_true(1)=rwage(t)*unemployment_subsidy;
       q=1;
        while abs(pay1_true(q)-rw_true)>0.0005
           q=q+1;
           pay1_true(q)=rw_true+coeffw*(pay1_true(q-1)-rw_true);
           pay2_true(q)=rw_true*unemployment_subsidy+coeffw*(pay2_true(q-1)-rw_true*unemployment_subsidy);
        end
       tMatE_true=[prob_EE_p(t) prob_EU_p(t); prob_UE_p(t) prob_UU_p(t)];
       if isinf(cond(tMatE_true))
       YPE_true(t)=rw_true;
       deviation(1:W,t)=0;
       else
       permwt=projectIncome(prob_EE_p(t),prob_EU_p(t),prob_UE_p(t),prob_UU_p(t),pay1_true,pay2_true,r_d);
       YPE_true(t)=permwt(1)/(1/r_d);
       for i=1:W
           deviation(i,t)=(YPE_true(t)-YPE_belief(i,t))/YPE_belief(i,t);
            if types(i)==1
                switch_index=0.95*deviation(i,t)+0.05*(opt_share(t-1)-(pes_share(t-1)+neut_share(t-1)));
                switch_index=1/(1+exp(20*switch_index+4));                
                if switch_index>=randswitch_w(t,i)
                   types(i)=0;
                end
            else
                if types(i)==0
                    if deviation(i,t)>=0
                        switch_index=0.95*deviation(i,t)+0.05*(opt_share(t-1)-(pes_share(t-1)+neut_share(t-1)));
                        switch_index=1/(1+exp(20*(-switch_index)+4));                        
                        if switch_index>=randswitch_w(t,i)
                            types(i)=1;
                        end
                    else
                        switch_index=0.95*deviation(i,t)+0.05*(neut_share(t-1)+opt_share(t-1)-pes_share(t-1));
                        switch_index=1/(1+exp(20*switch_index+4));
                        if switch_index>=randswitch_w(t,i)
                            types(i)=-1;
                        end
                    end
                else
                switch_index=0.95*deviation(i,t)+0.05*((opt_share(t-1)+neut_share(t-1))-pes_share(t-1));    
                switch_index=1/(1+exp(20*(-switch_index)+4));
                if switch_index>=randswitch_w(t,i)
                    types(i)=0;
                end
                end
            end
       end
       end
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
    cons_budget = max(0,cons_budget);
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
            if Y(best)==0
                noprod=noprod+1;
            end
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
    actualEXP(t)=sum(public_dem_c(:,t));
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
        pi(i)=pi(i)-di;
        totdiv(i)=totdiv(i)+di;
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
        pi_k(i)=pi_k(i)-di;
        totdivk(i)=totdivk(i)+di;
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

%disp(sum(deb)+sum(deb_k)+stock_bonds(t)-sum(PA)-sum(liquidity)-sum(liquidity_k)-E)
%disp(sum(PA)+sum(A)+sum(A_k)+E-stock_bonds(t)-sum(capital_value)-sum(Y_k.*P_k))

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
    
    
    if t<shockperiod && opinions==1
        prob_DD_p(t)=prob_DD(t);
        prob_DN_p(t)=prob_DN(t);
        prob_NN_p(t)=prob_NN(t);
        prob_ND_p(t)=prob_ND(t);
        prob_DD_bp(t)=prob_DD_b(t);
        prob_DN_bp(t)=prob_DN_b(t);
        prob_NN_bp(t)=prob_NN_b(t);
        prob_ND_bp(t)=prob_ND_b(t);
    end
    
    if t>=shockperiod && opinions==1
       if divstatus_b==1 && divstatus_b_prev==0
           trans_ND_bp(t)=trans_ND_bp(t)+1;
       end
       if divstatus_b==1 && divstatus_b_prev==1
           trans_DD_bp(t)=trans_DD_bp(t)+1;
       end
       if divstatus_b==0 && divstatus_b_prev==0
           trans_NN_bp(t)=trans_NN_bp(t)+1;
       end
       if divstatus_b==0 && divstatus_b_prev==1
           trans_DN_bp(t)=trans_DN_bp(t)+1;
       end
       
       if sum(trans_DD_bp(shockperiod:t))+sum(trans_DN_bp(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_DD_bp(t)=sum(trans_DD_bp(shockperiod:t))/(sum(trans_DD_bp(shockperiod:t))+sum(trans_DN_bp(shockperiod:t)));
            prob_DN_bp(t)=sum(trans_DN_bp(shockperiod:t))/(sum(trans_DD_bp(shockperiod:t))+sum(trans_DN_bp(shockperiod:t)));
            else
            prob_DD_bp(t)=sum(trans_DD_bp((t-(horizon-1)):t))/(sum(trans_DD_bp((t-(horizon-1)):t))+sum(trans_DN_bp((t-(horizon-1)):t)));
            prob_DN_bp(t)=sum(trans_DN_bp((t-(horizon-1)):t))/(sum(trans_DD_bp((t-(horizon-1)):t))+sum(trans_DN_bp((t-(horizon-1)):t)));
            end
        else
            prob_DD_bp(t)=0;
            prob_DN_bp(t)=0;
        end

        if sum(trans_ND_bp(shockperiod:t))+sum(trans_NN_bp(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_NN_bp(t)=sum(trans_NN_bp(shockperiod:t))/(sum(trans_ND_bp(shockperiod:t))+sum(trans_NN_bp(shockperiod:t)));
            prob_ND_bp(t)=sum(trans_ND_bp(shockperiod:t))/(sum(trans_ND_bp(shockperiod:t))+sum(trans_NN_bp(shockperiod:t)));
            else
            prob_NN_bp(t)=sum(trans_NN_bp((t-(horizon-1)):t))/(sum(trans_ND_bp((t-(horizon-1)):t))+sum(trans_NN_bp((t-(horizon-1)):t)));
            prob_ND_bp(t)=sum(trans_ND_bp((t-(horizon-1)):t))/(sum(trans_ND_bp((t-(horizon-1)):t))+sum(trans_NN_bp((t-(horizon-1)):t)));  
            end
        else
            prob_NN_bp(t)=0;
            prob_ND_bp(t)=0;
        end
       
        
        for i=1:F
           if divstatus(i)==1 && divstatus_prev(i)==0
               trans_ND_p(t)=trans_ND_p(t)+1;
           end
           if divstatus(i)==1 && divstatus_prev(i)==1
               trans_DD_p(t)=trans_DD_p(t)+1;
           end
           if divstatus(i)==0 && divstatus_prev(i)==0
               trans_NN_p(t)=trans_NN_p(t)+1;
           end
           if divstatus(i)==0 && divstatus_prev(i)==1
               trans_DN_p(t)=trans_DN_p(t)+1;
           end
        end

        if sum(trans_DD_p(shockperiod:t))+sum(trans_DN_p(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_DD_p(t)=sum(trans_DD_p(shockperiod:t))/(sum(trans_DD_p(shockperiod:t))+sum(trans_DN_p(shockperiod:t)));
            prob_DN_p(t)=sum(trans_DN_p(shockperiod:t))/(sum(trans_DD_p(shockperiod:t))+sum(trans_DN_p(shockperiod:t)));
            else
            prob_DD_p(t)=sum(trans_DD_p((t-(horizon-1)):t))/(sum(trans_DD_p((t-(horizon-1)):t))+sum(trans_DN_p((t-(horizon-1)):t)));
            prob_DN_p(t)=sum(trans_DN_p((t-(horizon-1)):t))/(sum(trans_DD_p((t-(horizon-1)):t))+sum(trans_DN_p((t-(horizon-1)):t)));
            end
        else
            prob_DD_p(t)=0;
            prob_DN_p(t)=0;
        end

        if sum(trans_ND_p(shockperiod:t))+sum(trans_NN_p(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_NN_p(t)=sum(trans_NN_p(shockperiod:t))/(sum(trans_ND_p(shockperiod:t))+sum(trans_NN_p(shockperiod:t)));
            prob_ND_p(t)=sum(trans_ND_p(shockperiod:t))/(sum(trans_ND_p(shockperiod:t))+sum(trans_NN_p(shockperiod:t)));
            else
            prob_NN_p(t)=sum(trans_NN_p((t-(horizon-1)):t))/(sum(trans_ND_p((t-(horizon-1)):t))+sum(trans_NN_p((t-(horizon-1)):t)));
            prob_ND_p(t)=sum(trans_ND_p((t-(horizon-1)):t))/(sum(trans_ND_p((t-(horizon-1)):t))+sum(trans_NN_p((t-(horizon-1)):t)));  
            end
        else
            prob_NN_p(t)=0;
            prob_ND_p(t)=0;
        end
    end
    
    tMatB_true=[prob_DD_bp(t) prob_DN_bp(t); prob_ND_bp(t) prob_NN_bp(t)];
    if isinf(prob_DD_bp(t))||isinf(prob_DN_bp(t))||isinf(prob_ND_bp(t))||isinf(prob_NN_bp(t))
        tMatB_true=[0 0; 0 0]; 
    end
    if isnan(prob_DD_bp(t))||isnan(prob_DN_bp(t))||isnan(prob_ND_bp(t))||isnan(prob_NN_bp(t))
        tMatB_true=[0 0; 0 0]; 
    end
    
    
    tMatD_true=[prob_DD_p(t) prob_DN_p(t); prob_ND_p(t) prob_NN_p(t)];
    if isinf(prob_DD_p(t))||isinf(prob_DN_p(t))||isinf(prob_ND_p(t))||isinf(prob_NN_p(t))
        tMatD_true=[0 0; 0 0]; 
    end
    if isnan(prob_DD_p(t))||isnan(prob_DN_p(t))||isnan(prob_ND_p(t))||isnan(prob_NN_p(t))
        tMatD_true=[0 0; 0 0]; 
    end
    
    rdivs(:,t)=dividends_income_pb./price(t);
    rdivs_pos(t)=sum(rdivs(:,t))/sum(rdivs(:,t)>0);
    rdivs_b(t)=(dividendsB(t)/(F+N))/price(t);
    
    if t>shockperiod && opinions==1
        if isinf(cond(tMatB_true))
            impulse_DN_b_opt(t)=impulse_DN_b_opt(t-1);
            impulse_NN_b_opt(t)=impulse_NN_b_opt(t-1);
            impulse_DN_b_pes(t)=impulse_DN_b_pes(t-1);
            impulse_NN_b_pes(t)=impulse_NN_b_pes(t-1);    
        else
            impulse_DN_b_opt(t)=max(0,adapt*impulse_DN_b_opt(t-1)+(1-adapt)*(prob_DN_bp(t)/prob_DN_b_opt(t-1)));
            impulse_NN_b_opt(t)=max(0,adapt*impulse_NN_b_opt(t-1)+(1-adapt)*(prob_NN_bp(t)/prob_NN_b_opt(t-1)));
            impulse_DN_b_pes(t)=max(0,adapt*impulse_DN_b_pes(t-1)+(1-adapt)*(prob_DN_bp(t)/prob_DN_b_pes(t-1)));
            impulse_NN_b_pes(t)=max(0,adapt*impulse_NN_b_pes(t-1)+(1-adapt)*(prob_NN_bp(t)/prob_NN_b_pes(t-1)));
            if impulse_DN_b_opt(t)>1 || impulse_DN_b_opt(t-1)==1
               impulse_DN_b_opt(t)=1;
            end
            if impulse_NN_b_opt(t)>1 || impulse_NN_b_opt(t-1)==1
               impulse_NN_b_opt(t)=1;
            end
            if impulse_DN_b_pes(t)<1 || impulse_DN_b_pes(t-1)==1
               impulse_DN_b_pes(t)=1;
            end
            if impulse_NN_b_pes(t)<1 || impulse_NN_b_pes(t-1)==1
               impulse_NN_b_pes(t)=1;
            end
        end
        
        
        if isinf(cond(tMatD_true))
            impulse_DN_opt(t)=impulse_DN_opt(t-1);
            impulse_NN_opt(t)=impulse_NN_opt(t-1);
            impulse_DN_pes(t)=impulse_DN_pes(t-1);
            impulse_NN_pes(t)=impulse_NN_pes(t-1);    
        else
            impulse_DN_opt(t)=max(0,adapt*impulse_DN_opt(t-1)+(1-adapt)*(prob_DN_p(t)/prob_DN_opt(t-1)));
            impulse_NN_opt(t)=max(0,adapt*impulse_NN_opt(t-1)+(1-adapt)*(prob_NN_p(t)/prob_NN_opt(t-1)));
            impulse_DN_pes(t)=max(0,adapt*impulse_DN_pes(t-1)+(1-adapt)*(prob_DN_p(t)/prob_DN_pes(t-1)));
            impulse_NN_pes(t)=max(0,adapt*impulse_NN_pes(t-1)+(1-adapt)*(prob_NN_p(t)/prob_NN_pes(t-1)));
            if impulse_DN_opt(t)>1 || impulse_DN_opt(t-1)==1
               impulse_DN_opt(t)=1;
            end
            if impulse_NN_opt(t)>1 || impulse_NN_opt(t-1)==1
               impulse_NN_opt(t)=1;
            end
            if impulse_DN_pes(t)<1 || impulse_DN_pes(t-1)==1
               impulse_DN_pes(t)=1;
            end
            if impulse_NN_pes(t)<1 || impulse_NN_pes(t-1)==1
               impulse_NN_pes(t)=1;
            end
        end
        
        if rdivs_b(t-1)>0
            impulse_rdb_opt(t)=adapt*impulse_rdb_opt(t-1)+(1-adapt)*min(2,(rdivs_b(t)/(rdivs_b(t-1)*impulse_rdb_opt(t-1))));
            impulse_rdb_pes(t)=adapt*impulse_rdb_pes(t-1)+(1-adapt)*min(2,(rdivs_b(t)/(rdivs_b(t-1)*impulse_rdb_pes(t-1))));
        else
            impulse_rdb_opt(t)=impulse_rdb_opt(t-1);
            impulse_rdb_pes(t)=impulse_rdb_pes(t-1);
        end
        if impulse_rdb_opt(t-1)>=1 && impulse_rdb_opt(t)<=1 || impulse_rdb_opt(t-1)<=1 && impulse_rdb_opt(t)>=1
            impulse_rdb_opt(t)=1;
        end
        if impulse_rdb_pes(t-1)>=1 && impulse_rdb_pes(t)<=1 || impulse_rdb_pes(t-1)<=1 && impulse_rdb_pes(t)>=1
            impulse_rdb_pes(t)=1;
        end
        
        
        for i=1:F
        if types(W+i)==1 || types(W+i)==-1
            if rdivs(i,t-1)>0
            impulse_rd(i,t)=adapt*impulse_rd(i,t-1)+(1-adapt)*min(2,(rdivs(i,t)/(rdivs(i,t-1)*impulse_rd(i,t-1))));
            else
            impulse_rd(i,t)=impulse_rd(i,t-1);
            end
        end
        if impulse_rd(i,t-1)>=1 && impulse_rd(i,t)<=1 || impulse_rd(i,t-1)<=1 && impulse_rd(i,t)>=1
            impulse_rd(i,t)=1;
        end
        end
    end
    
    prob_DN_b_opt(t)=max(0,min(1,prob_DN_b(t)*impulse_DN_b_opt(t)));
    prob_DD_b_opt(t)=1-prob_DN_b_opt(t);
    prob_DN_b_pes(t)=max(0,min(1,prob_DN_b(t)*impulse_DN_b_pes(t)));
    prob_DD_b_pes(t)=1-prob_DN_b_pes(t);
    prob_NN_b_opt(t)=max(0,min(1,prob_NN_b(t)*impulse_NN_b_opt(t)));
    prob_ND_b_opt(t)=1-prob_NN_b_opt(t);
    prob_NN_b_pes(t)=max(0,min(1,prob_NN_b(t)*impulse_NN_b_pes(t)));
    prob_ND_b_pes(t)=1-prob_NN_b_pes(t);
    
    prob_DN_opt(t)=max(0,min(1,prob_DN(t)*impulse_DN_opt(t)));
    prob_DD_opt(t)=1-prob_DN_opt(t);
    prob_DN_pes(t)=max(0,min(1,prob_DN(t)*impulse_DN_pes(t)));
    prob_DD_pes(t)=1-prob_DN_pes(t);
    prob_NN_opt(t)=max(0,min(1,prob_NN(t)*impulse_NN_opt(t)));
    prob_ND_opt(t)=1-prob_NN_opt(t);
    prob_NN_pes(t)=max(0,min(1,prob_NN(t)*impulse_NN_pes(t)));
    prob_ND_pes(t)=1-prob_NN_pes(t);
    
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
            rd_e(i)=rd(i)*impulse_rd(i,t);
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
        cdb(t)=coeffdb;
        
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
        
        permdb=projectIncome(prob_DD_b(t),prob_DN_b(t),prob_ND_b(t),prob_NN_b(t),pay1b,pay2b,r_d);
        
        if opinions==1
            if t>=shockperiod && t< shockperiod+duration
            pay1b_opt=pay1b;
            pay2b_opt=pay2b;
            pay1b_opt(1:((shockperiod+duration)-t))=pay1b_opt(1:((shockperiod+duration)-t))*impulse_rdb_opt(t);
            pay2b_opt(1:((shockperiod+duration)-t))=pay2b_opt(1:((shockperiod+duration)-t))*impulse_rdb_opt(t);
            probs11=zeros(1,length(pay1b_opt));
            probs12=zeros(1,length(pay1b_opt));
            probs21=zeros(1,length(pay1b_opt));
            probs22=zeros(1,length(pay1b_opt));
            probs11(:)=prob_DD_b(t);
            probs12(:)=prob_DN_b(t);
            probs21(:)=prob_ND_b(t);
            probs22(:)=prob_NN_b(t);
            probs11(1:((shockperiod+duration)-t))=prob_DD_b_opt(t);
            probs12(1:((shockperiod+duration)-t))=prob_DN_b_opt(t);
            probs21(1:((shockperiod+duration)-t))=prob_ND_b_opt(t);
            probs22(1:((shockperiod+duration)-t))=prob_NN_b_opt(t);
            permdbopt=projectIncome(probs11,probs12,probs21,probs22,pay1b_opt,pay2b_opt,r_d);
            pay1b_pes=pay1b;
            pay2b_pes=pay2b;
            pay1b_pes(1:((shockperiod+duration)-t))=pay1b_pes(1:((shockperiod+duration)-t))*impulse_rdb_pes(t);
            pay2b_pes(1:((shockperiod+duration)-t))=pay2b_pes(1:((shockperiod+duration)-t))*impulse_rdb_pes(t);
            probs11=zeros(1,length(pay1b_pes));
            probs12=zeros(1,length(pay1b_pes));
            probs21=zeros(1,length(pay1b_pes));
            probs22=zeros(1,length(pay1b_pes));
            probs11(:)=prob_DD_b(t);
            probs12(:)=prob_DN_b(t);
            probs21(:)=prob_ND_b(t);
            probs22(:)=prob_NN_b(t);
            probs11(1:((shockperiod+duration)-t))=prob_DD_b_pes(t);
            probs12(1:((shockperiod+duration)-t))=prob_DN_b_pes(t);
            probs21(1:((shockperiod+duration)-t))=prob_ND_b_pes(t);
            probs22(1:((shockperiod+duration)-t))=prob_NN_b_pes(t);
            permdbpes=projectIncome(probs11,probs12,probs21,probs22,pay1b_pes,pay2b_pes,r_d);
            end
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
        cd(:,t)=coeffsd;
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
        permds=projectIncomeMult(prob_DD(t),prob_DN(t),prob_ND(t),prob_NN(t),pay1,pay2,r_d);
        for i=1:F
        permd=permds(:,i);
        V_div=permd;
            if dividends_income(i)>0
                if dividendsB(t)>0
                V_div=V_div(1)+permdb(1);
                else
                V_div=V_div(1)+permdb(2);
                end
                yd(i)=(PA(W+i)/price(t)+V_div)/(1/r_d);
                yd_belief(i,t)=(V_div)/(1/r_d);
            else
                if dividendsB(t)>0
                V_div=V_div(2)+permdb(1);
                else
                V_div=V_div(2)+permdb(2);
                end
                yd(i)=(PA(W+i)/price(t)+V_div)/(1/r_d);
                yd_belief(i,t)=(V_div)/(1/r_d);
            end
        end
    
    if opinions==1 && t<shockperiod && t>=announce
        pay1_e=pay1;
        pay2_e=pay2;
        for i=1:F
        pay1_e(shockperiod-t:shockperiod+duration-t,i)=pay1_e(shockperiod-t:shockperiod+duration-t,i)*impulse_rd(i,shockperiod);
        end
        probs11opt=zeros(1,length(pay1_e(:,1)));
        probs12opt=zeros(1,length(pay1_e(:,1)));
        probs21opt=zeros(1,length(pay1_e(:,1)));
        probs22opt=zeros(1,length(pay1_e(:,1)));
        probs11opt(:)=prob_DD(t);
        probs12opt(:)=prob_DN(t);
        probs21opt(:)=prob_ND(t);
        probs22opt(:)=prob_NN(t);
        probs12opt(shockperiod-t:shockperiod+duration-t)=prob_DN(t)*impulse_DN_opt(shockperiod);
        probs11opt(shockperiod-t:shockperiod+duration-t)=1-probs12opt(shockperiod-t:shockperiod+duration-t);
        probs22opt(shockperiod-t:shockperiod+duration-t)=prob_NN(t)*impulse_NN_opt(shockperiod);
        probs21opt(shockperiod-t:shockperiod+duration-t)=1-probs22opt(shockperiod-t:shockperiod+duration-t);
        probs11pes=probs11opt;
        probs12pes=probs12opt;
        probs21pes=probs21opt;
        probs22pes=probs22opt;
        probs12pes(shockperiod-t:shockperiod+duration-t)=prob_DN(t)*impulse_DN_pes(shockperiod);
        probs11pes(shockperiod-t:shockperiod+duration-t)=1-probs12pes(shockperiod-t:shockperiod+duration-t);
        probs22pes(shockperiod-t:shockperiod+duration-t)=prob_NN(t)*impulse_NN_pes(shockperiod);
        probs21pes(shockperiod-t:shockperiod+duration-t)=1-probs22pes(shockperiod-t:shockperiod+duration-t);
        permds_opt=projectIncomeMult(probs11opt,probs12opt,probs21opt,probs22opt,pay1_e,pay2_e,r_d);
        permds_pes=projectIncomeMult(probs11pes,probs12pes,probs21pes,probs22pes,pay1_e,pay2_e,r_d);
        for i=1:F
        if types(W+i)==1
            V_div=permds_opt(:,i);
            if dividendsB(t)>0
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbopt(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbopt(1))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbopt(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbopt(1))/(1/r_d);
                end
            else
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbopt(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbopt(2))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbopt(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbopt(2))/(1/r_d);
                end
            end
        end
        if types(W+i)==-1
            V_div=permds_pes(:,i);
            if dividendsB(t)>0
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbpes(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbpes(1))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbpes(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbpes(1))/(1/r_d);
                end
            else
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbpes(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbpes(2))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbpes(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbpes(2))/(1/r_d);
                end
            end
        end
        end
    end
        
        
    if opinions==1 && t>=shockperiod && t<shockperiod+duration
        pay1_e=pay1;
        pay2_e=pay2;
        for i=1:F
        pay1_e(1:((shockperiod+duration)-t),i)=pay1_e(1:((shockperiod+duration)-t),i)*impulse_rd(i,t);
        end
        probs11opt=zeros(1,length(pay1_e(:,1)));
        probs12opt=zeros(1,length(pay1_e(:,1)));
        probs21opt=zeros(1,length(pay1_e(:,1)));
        probs22opt=zeros(1,length(pay1_e(:,1)));
        probs11opt(:)=prob_DD(t);
        probs12opt(:)=prob_DN(t);
        probs21opt(:)=prob_ND(t);
        probs22opt(:)=prob_NN(t);
        probs11opt(1:((shockperiod+duration)-t))=prob_DD_opt(t);
        probs12opt(1:((shockperiod+duration)-t))=prob_DN_opt(t);
        probs21opt(1:((shockperiod+duration)-t))=prob_ND_opt(t);
        probs22opt(1:((shockperiod+duration)-t))=prob_NN_opt(t);
        probs11pes=probs11opt;
        probs12pes=probs12opt;
        probs21pes=probs21opt;
        probs22pes=probs22opt;
        probs11pes(1:((shockperiod+duration)-t))=prob_DD_pes(t);
        probs12pes(1:((shockperiod+duration)-t))=prob_DN_pes(t);
        probs21pes(1:((shockperiod+duration)-t))=prob_ND_pes(t);
        probs22pes(1:((shockperiod+duration)-t))=prob_NN_pes(t);
        permds_opt=projectIncomeMult(probs11opt,probs12opt,probs21opt,probs22opt,pay1_e,pay2_e,r_d);
        permds_pes=projectIncomeMult(probs11pes,probs12pes,probs21pes,probs22pes,pay1_e,pay2_e,r_d);
        for i=1:F
        if types(W+i)==1
            V_div=permds_opt(:,i);
            if dividendsB(t)>0
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbopt(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbopt(1))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbopt(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbopt(1))/(1/r_d);
                end
            else
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbopt(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbopt(2))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbopt(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbopt(2))/(1/r_d);
                end
            end
        end
        if types(W+i)==-1
            V_div=permds_pes(:,i);
            if dividendsB(t)>0
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbpes(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbpes(1))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbpes(1))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbpes(1))/(1/r_d);
                end
            else
                if dividends_income(i)>0
                    yd(i)=(PA(W+i)/price(t)+V_div(1)+permdbpes(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(1)+permdbpes(2))/(1/r_d);
                else
                    yd(i)=(PA(W+i)/price(t)+V_div(2)+permdbpes(2))/(1/r_d);
                    yd_belief(i,t)=(V_div(2)+permdbpes(2))/(1/r_d);
                end
            end
        end
        end
    end
    
    if opinions==1
       if t>=shockperiod
           if (t-shockperiod)<burnin
                if sum(rdivs_b(shockperiod:t)>0)>0
                    rdb_true=sum(rdivs_b(shockperiod:t))/sum(rdivs_b(shockperiod:t)>0);
                else
                    rdb_true=0;
                end
                for i=1:F
                if sum(rdivs(i,shockperiod:t)>0)>0
                    rd_true(i)=sum(rdivs(i,shockperiod:t))/sum(rdivs(i,shockperiod:t)>0);
                else
                    rd_true(i)=0;
                end
                end
           else
                rdb_true=sum(rdivs_b(t-(horizon-1):t))/sum(rdivs_b(t-(horizon-1):t)>0);
                for i=1:F
                rd_true(i)=sum(rdivs(i,t-(horizon-1):t))/sum(rdivs(i,t-(horizon-1):t)>0);
                end
           end
       end
    end

    if t>=shockperiod && switching==1
        if ~isinf(cond(tMatB_true))
           pay1b_true=rdb_true*ones(2000,1);
           pay2b_true=zeros(2000,1);
           currentpay=rdivs_b(shockperiod:t);
           currentpay=currentpay(currentpay>0);
           if ~isempty(currentpay)
           pay1b_true(1)=currentpay(length(currentpay));
           else
           pay1b_true(1)=rdivs_b(t);
           end
           q=1;
           while abs(pay1b_true(q)-rdb_true)>0.0005
               q=q+1;
               pay1b_true(q)=rdb_true+coeffdb*(pay1b_true(q-1)-rdb_true);
           end
           permdbt=projectIncome(prob_DD_bp(t),prob_DN_bp(t),prob_ND_bp(t),prob_NN_bp(t),pay1b_true,pay2b_true,r_d);
        end
        
        
        if ~isinf(cond(tMatD_true))
        md_true=mean(rd_true);
        pay1_true=md_true*ones(2000,F);
        currentpay=rdivs(:,t);
        currentpay(currentpay==0)=rdivs_pos(t);
        pay1_true(1,:)=currentpay;
        q=1;
        while abs(mean(pay1_true(q,:))-md_true)>0.0005
           q=q+1;
           pay1_true(q,:)=md_true+coeffsd'.*(pay1_true(q-1,:)-md_true);
        end
        pay2_true=zeros(length(pay1),F);
        permdts=projectIncomeMult(prob_DD_p(t),prob_DN_p(t),prob_ND_p(t),prob_NN_p(t),pay1_true,pay2_true,r_d); 
        end
        
        if ~isinf(cond(tMatB_true)) && ~isinf(cond(tMatD_true))
        for i=1:F
        V_div_true=permdts(:,i);
        if dividends_income(i)>0
            if dividendsB(t)>0    
                yd_true=(V_div_true(1)+permdbt(1))/(1/r_d);
            else
                yd_true=(V_div_true(1)+permdbt(2))/(1/r_d);
            end
        else
            if dividendsB(t)>0    
                yd_true=(V_div_true(2)+permdbt(1))/(1/r_d);
            else
                yd_true=(V_div_true(2)+permdbt(2))/(1/r_d);
            end  
        end
        deviation(W+i,t)=(yd_true-yd_belief(i,t))/yd_belief(i,t);       
            if types(W+i)==1
                switch_index=0.99*deviation(W+i,t)+0.01*(opt_share(t-1)-(pes_share(t-1)+neut_share(t-1)));
                switch_index=1/(1+exp(15*switch_index+4));
                if switch_index>=randswitch_f(t,i)
                   types(W+i)=0;
                end
            else
                if types(W+i)==0
                    if deviation(W+i,t)>=0
                        switch_index=0.99*deviation(W+i,t)+0.01*(opt_share(t-1)-(pes_share(t-1)+neut_share(t-1)));
                        switch_index=1/(1+exp(15*(-switch_index)+4));
                        if switch_index>=randswitch_f(t,i)
                            types(W+i)=1;
                        end
                    else
                        switch_index=0.99*deviation(W+i,t)+0.01*(neut_share(t-1)+opt_share(t-1)-pes_share(t-1));
                        switch_index=1/(1+exp(15*switch_index+4));
                        if switch_index>=randswitch_f(t,i)
                            types(W+i)=-1;
                        end
                    end
                else
                switch_index=0.99*deviation(W+i,t)+0.01*((opt_share(t-1)+neut_share(t-1))-pes_share(t-1));
                switch_index=1/(1+exp(15*(-switch_index)+4));
                if switch_index>=randswitch_f(t,i)
                    types(W+i)=0;
                end
                end
            end
        end
        end       
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

    if t<burnin
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
    
    if t<shockperiod && opinions==1
        prob_DD_kp(t)=prob_DD_k(t);
        prob_DN_kp(t)=prob_DN_k(t);
        prob_NN_kp(t)=prob_NN_k(t);
        prob_ND_kp(t)=prob_ND_k(t);
    end
    
    if t>=shockperiod && opinions==1
        for l=1:N
        if divstatus_k(l)==1 && divstatus_prev_k(l)==0
           trans_ND_kp(t)=trans_ND_kp(t)+1;
       end
       if divstatus_k(l)==1 && divstatus_prev_k(l)==1
           trans_DD_kp(t)=trans_DD_kp(t)+1;
       end
       if divstatus_k(l)==0 && divstatus_prev_k(l)==0
           trans_NN_kp(t)=trans_NN_kp(t)+1;
       end
       if divstatus_k(l)==0 && divstatus_prev_k(l)==1
           trans_DN_kp(t)=trans_DN_kp(t)+1;
       end
        end

        if sum(trans_DD_kp(shockperiod:t))+sum(trans_DN_kp(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_DD_kp(t)=sum(trans_DD_kp(shockperiod:t))/(sum(trans_DD_kp(shockperiod:t))+sum(trans_DN_kp(shockperiod:t)));
            prob_DN_kp(t)=sum(trans_DN_kp(shockperiod:t))/(sum(trans_DD_kp(shockperiod:t))+sum(trans_DN_kp(shockperiod:t)));
            else
            prob_DD_kp(t)=sum(trans_DD_kp((t-(horizon-1)):t))/(sum(trans_DD_kp((t-(horizon-1)):t))+sum(trans_DN_kp((t-(horizon-1)):t)));
            prob_DN_kp(t)=sum(trans_DN_kp((t-(horizon-1)):t))/(sum(trans_DD_kp((t-(horizon-1)):t))+sum(trans_DN_kp((t-(horizon-1)):t)));
            end
        else
            prob_DN_kp(t)=0;
            prob_DN_kp(t)=0;
        end

        if sum(trans_ND_kp(shockperiod:t))+sum(trans_NN_kp(shockperiod:t))>0
            if (t-shockperiod)<burnin
            prob_NN_kp(t)=sum(trans_NN_kp(shockperiod:t))/(sum(trans_ND_kp(shockperiod:t))+sum(trans_NN_kp(shockperiod:t)));
            prob_ND_kp(t)=sum(trans_ND_kp(shockperiod:t))/(sum(trans_ND_kp(shockperiod:t))+sum(trans_NN_kp(shockperiod:t)));
            else
            prob_NN_kp(t)=sum(trans_NN_kp((t-(horizon-1)):t))/(sum(trans_ND_kp((t-(horizon-1)):t))+sum(trans_NN_kp((t-(horizon-1)):t)));
            prob_ND_kp(t)=sum(trans_ND_kp((t-(horizon-1)):t))/(sum(trans_ND_kp((t-(horizon-1)):t))+sum(trans_NN_kp((t-(horizon-1)):t)));  
            end
        else
            prob_NN_kp(t)=0;
            prob_ND_kp(t)=0;
        end
    end
    
    tMatDk_true=[prob_DD_kp(t) prob_DN_kp(t); prob_ND_kp(t) prob_NN_kp(t)];
    if isinf(prob_DD_kp(t))||isinf(prob_DN_kp(t))||isinf(prob_ND_kp(t))||isinf(prob_NN_kp(t))
    tMatDk_true=[0 0; 0 0]; 
    end
    if isnan(prob_DD_kp(t))||isnan(prob_DN_kp(t))||isnan(prob_ND_kp(t))||isnan(prob_NN_kp(t))
    tMatDk_true=[0 0; 0 0]; 
    end
    
    rdivs_k(:,t)=dividends_income_pbk./price(t);
    rdivs_kpos(t)=sum(rdivs_k(:,t))/sum(rdivs_k(:,t)>0);
    
    if t>shockperiod && opinions==1
        if isinf(cond(tMatDk_true))
            impulse_DN_k_opt(t)=impulse_DN_k_opt(t-1);
            impulse_NN_k_opt(t)=impulse_NN_k_opt(t-1);
            impulse_DN_k_pes(t)=impulse_DN_k_pes(t-1);
            impulse_NN_k_pes(t)=impulse_NN_k_pes(t-1);    
        else
            impulse_DN_k_opt(t)=max(0,adapt*impulse_DN_k_opt(t-1)+(1-adapt)*(prob_DN_kp(t)/prob_DN_k_opt(t-1)));
            impulse_NN_k_opt(t)=max(0,adapt*impulse_NN_k_opt(t-1)+(1-adapt)*(prob_NN_kp(t)/prob_NN_k_opt(t-1)));
            impulse_DN_k_pes(t)=max(0,adapt*impulse_DN_k_pes(t-1)+(1-adapt)*(prob_DN_kp(t)/prob_DN_k_pes(t-1)));
            impulse_NN_k_pes(t)=max(0,adapt*impulse_NN_k_pes(t-1)+(1-adapt)*(prob_NN_kp(t)/prob_NN_k_pes(t-1)));
            if impulse_DN_k_opt(t)>1 || impulse_DN_k_opt(t-1)==1
               impulse_DN_k_opt(t)=1;
            end
            if impulse_NN_k_opt(t)>1 || impulse_NN_k_opt(t-1)==1
               impulse_NN_k_opt(t)=1;
            end
            if impulse_DN_k_pes(t)<1 || impulse_DN_k_pes(t-1)==1
               impulse_DN_k_pes(t)=1;
            end
            if impulse_NN_k_pes(t)<1 || impulse_NN_k_pes(t-1)==1
               impulse_NN_k_pes(t)=1;
            end
        end
        for l=1:N
        if types(W+F+l)==1 || types(W+F+l)==-1    
            if rdivs_k(l,t-1)>0
            impulse_rdk(l,t)=adapt*impulse_rdk(l,t-1)+(1-adapt)*min(2,(rdivs(l,t)/(rdivs_k(l,t-1)*impulse_rdk(l,t-1))));
            else
            impulse_rdk(l,t)=impulse_rdk(l,t-1);
            end
        else
            impulse_rdk(l,t)=1;
        end
        if impulse_rdk(l,t-1)<=1 && impulse_rdk(l,t)>=1 || impulse_rdk(l,t-1)>=1 && impulse_rdk(l,t)<=1
            impulse_rdk(l,t)=1;
        end
        end
    end
    
    prob_DN_k_opt(t)=max(0,min(1,prob_DN_k(t)*impulse_DN_k_opt(t)));
    prob_DD_k_opt(t)=1-prob_DN_k_opt(t);
    prob_NN_k_opt(t)=max(0,min(1,prob_NN_k(t)*impulse_NN_k_opt(t)));
    prob_ND_k_opt(t)=1-prob_NN_k_opt(t);
    prob_DN_k_pes(t)=max(0,min(1,prob_DN_k(t)*impulse_DN_k_pes(t)));
    prob_DD_k_pes(t)=1-prob_DN_k_pes(t);
    prob_NN_k_pes(t)=max(0,min(1,prob_NN_k(t)*impulse_NN_k_pes(t)));
    prob_ND_k_pes(t)=1-prob_NN_k_pes(t);

    if t < burnin
        for l=1:N
            if sum(rdivs_k(l,2:t)>0)>0
                rdk(l)=sum(rdivs_k(l,2:t))/sum(rdivs_k(l,2:t)>0);
            else
                rdk(l)=0;
            end
        end
    else
       for l=1:N
            rdk(l)=sum(rdivs_k(l,(t-(horizon-1)):t))/sum(rdivs_k(l,(t-(horizon-1)):t)>0);
            rdk_e(l)=rdk(l)*impulse_rdk(l,t);
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
        cdk(:,t)=coeffsdk;
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
        permdks=projectIncomeMult(prob_DD_k(t),prob_DN_k(t),prob_ND_k(t),prob_NN_k(t),pay1,pay2,r_d);
        for i=1:N
        permdk=permdks(:,i);
        V_divk=permdk;
            if dividends_income_k(i)>0
                if dividendsB(t)>0
                V_divk=V_divk(1)+permdb(1);
                else
                V_divk=V_divk(1)+permdb(2);
                end
                ydk(i)=(PA(W+F+i)/price(t)+V_divk)/(1/r_d);
                ydk_belief(i,t)=(V_divk)/(1/r_d);
            else
                if dividendsB(t)>0
                V_divk=V_divk(2)+permdb(1);
                else
                V_divk=V_divk(2)+permdb(2);
                end
                ydk(i)=(PA(W+F+i)/price(t)+V_divk)/(1/r_d);
                ydk_belief(i,t)=(V_divk)/(1/r_d);
            end
        end
    
    if opinions==1 && t<shockperiod && t>=announce
        pay1_e=pay1;
        pay2_e=pay2;
        for i=1:N
        pay1_e(shockperiod-t:shockperiod+duration-t,i)=pay1_e(shockperiod-t:shockperiod+duration-t,i)*impulse_rdk(i,shockperiod);
        end
        probs11opt=zeros(1,length(pay2_e));
        probs12opt=zeros(1,length(pay2_e));
        probs21opt=zeros(1,length(pay2_e));
        probs22opt=zeros(1,length(pay2_e));
        probs11opt(:)=prob_DD_k(t);
        probs12opt(:)=prob_DN_k(t);
        probs21opt(:)=prob_ND_k(t);
        probs22opt(:)=prob_NN_k(t);
        probs12opt(shockperiod-t:shockperiod+duration-t)=prob_DN_k(t)*impulse_DN_k_opt(shockperiod);
        probs11opt(shockperiod-t:shockperiod+duration-t)=1-probs12opt(shockperiod-t:shockperiod+duration-t);
        probs22opt(shockperiod-t:shockperiod+duration-t)=prob_NN_k_opt(t)*impulse_NN_k_opt(shockperiod);
        probs21opt(shockperiod-t:shockperiod+duration-t)=1-probs22opt(shockperiod-t:shockperiod+duration-t);
        probs11pes=probs11opt;
        probs12pes=probs12opt;
        probs21pes=probs21opt;
        probs22pes=probs22opt;
        probs12pes(shockperiod-t:shockperiod+duration-t)=prob_DN_k(t)*impulse_DN_k_pes(shockperiod);
        probs11pes(shockperiod-t:shockperiod+duration-t)=1-probs12pes(shockperiod-t:shockperiod+duration-t);
        probs22pes(shockperiod-t:shockperiod+duration-t)=prob_NN_k(t)*impulse_NN_k_pes(shockperiod);
        probs21pes(shockperiod-t:shockperiod+duration-t)=1-probs22pes(shockperiod-t:shockperiod+duration-t);
        permdks_opt=projectIncomeMult(probs11opt,probs12opt,probs21opt,probs22opt,pay1_e,pay2_e,r_d);
        permdks_pes=projectIncomeMult(probs11pes,probs12pes,probs21pes,probs22pes,pay1_e,pay2_e,r_d);
        for l=1:N
        if types(W+F+l)==1
            V_divk=permdks_opt(:,l);
            if dividends_income_k(l)>0
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbopt(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbopt(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbopt(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbopt(2))/(1/r_d);
                end
            else
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbopt(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbopt(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbopt(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbopt(2))/(1/r_d);
                end
            end
        end
        if types(W+F+l)==-1
            V_divk=permdks_pes(:,l);
            if dividends_income_k(l)>0
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbpes(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbpes(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbpes(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbpes(2))/(1/r_d);
                end
            else
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbpes(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbpes(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbpes(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbpes(2))/(1/r_d);
                end
            end
        end
        end
    end
        
    if opinions==1 && t>=shockperiod && t<shockperiod+duration
        pay1_e=pay1;
        pay2_e=pay2;
        for i=1:N
        pay1_e(1:((shockperiod+duration)-t),i)=pay1_e(1:((shockperiod+duration)-t),i)*impulse_rdk(i,t);
        end
        probs11opt=zeros(1,length(pay2_e));
        probs12opt=zeros(1,length(pay2_e));
        probs21opt=zeros(1,length(pay2_e));
        probs22opt=zeros(1,length(pay2_e));
        probs11opt(:)=prob_DD_k(t);
        probs12opt(:)=prob_DN_k(t);
        probs21opt(:)=prob_ND_k(t);
        probs22opt(:)=prob_NN_k(t);
        probs11opt(1:((shockperiod+duration)-t))=prob_DD_k_opt(t);
        probs12opt(1:((shockperiod+duration)-t))=prob_DN_k_opt(t);
        probs21opt(1:((shockperiod+duration)-t))=prob_ND_k_opt(t);
        probs22opt(1:((shockperiod+duration)-t))=prob_NN_k_opt(t);
        probs11pes=probs11opt;
        probs12pes=probs12opt;
        probs21pes=probs21opt;
        probs22pes=probs22opt;
        probs11pes(1:((shockperiod+duration)-t))=prob_DD_k_pes(t);
        probs12pes(1:((shockperiod+duration)-t))=prob_DN_k_pes(t);
        probs21pes(1:((shockperiod+duration)-t))=prob_ND_k_pes(t);
        probs22pes(1:((shockperiod+duration)-t))=prob_NN_k_pes(t);
        permdks_opt=projectIncomeMult(probs11opt,probs12opt,probs21opt,probs22opt,pay1_e,pay2_e,r_d);
        permdks_pes=projectIncomeMult(probs11pes,probs12pes,probs21pes,probs22pes,pay1_e,pay2_e,r_d);
        for l=1:N
        if types(W+F+l)==1
            V_divk=permdks_opt(:,l);
            if dividends_income_k(l)>0
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbopt(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbopt(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbopt(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbopt(2))/(1/r_d);
                end
            else
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbopt(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbopt(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbopt(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbopt(2))/(1/r_d);
                end
            end
        end
        if types(W+F+l)==-1
            V_divk=permdks_pes(:,l);
            if dividends_income_k(l)>0
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbpes(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbpes(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(1)+permdbpes(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(1)+permdbpes(2))/(1/r_d);
                end
            else
                if dividendsB(t)>0
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbpes(1))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbpes(1))/(1/r_d);
                else
                    ydk(l)=(PA(W+F+l)/price(t)+V_divk(2)+permdbpes(2))/(1/r_d);
                    ydk_belief(l,t)=(V_divk(2)+permdbpes(2))/(1/r_d);
                end
            end
        end
        end
    end
    
    
    if opinions==1
       if t>=shockperiod
           if (t-shockperiod)<burnin
            for i=1:N
            if sum(rdivs_k(i,shockperiod:t)>0)>0
                rdk_true(i)=sum(rdivs_k(i,shockperiod:t))/sum(rdivs_k(i,shockperiod:t)>0);
            else
                rdk_true(i)=0;
            end
            end
           else
            for i=1:N
            rdk_true(i)=sum(rdivs_k(i,t-(horizon-1):t))/sum(rdivs_k(i,t-(horizon-1):t)>0);
            end
           end
       end
    end
 

    if t>=shockperiod && switching==1
        if ~isinf(cond(tMatDk_true))
        mdk_true=mean(rdk_true);
        pay1_true=mdk_true*ones(2000,N);
        currentpay=rdivs_k(:,t);
        currentpay(currentpay==0)=rdivs_kpos(t);
        pay1_true(1,:)=currentpay;
        q=1;
        while abs(mean(pay1_true(q,:))-mdk_true)>0.0005
           q=q+1;
           pay1_true(q,:)=mdk_true+coeffsdk'.*(pay1_true(q-1,:)-mdk_true);
        end
        pay2_true=zeros(length(pay1),N);
        permdkt=projectIncomeMult(prob_DD_kp(t),prob_DN_kp(t),prob_ND_kp(t),prob_NN_kp(t),pay1_true,pay2_true,r_d); 
        end
        
        if ~isinf(cond(tMatB_true)) && ~isinf(cond(tMatDk_true))
        for l=1:N
        V_divk_true=permdkt(:,l);
        if dividends_income_k(l)>0
            if dividendsB(t)>0
                ydk_true=(V_divk_true(1)+permdbt(1))/(1/r_d);
            else
                ydk_true=(V_divk_true(1)+permdbt(2))/(1/r_d);
            end
        else
            if dividendsB(t)>0
                ydk_true=(V_divk_true(2)+permdbt(1))/(1/r_d);
            else
                ydk_true=(V_divk_true(2)+permdbt(2))/(1/r_d);
            end
        end
        deviation(W+F+l,t)=(ydk_true-ydk_belief(l,t))/ydk_belief(l,t);       
            if types(W+F+l)==1
                switch_index=0.99*deviation(W+F+l,t)+0.01*(opt_share(t-1)-(pes_share(t-1)+neut_share(t-1)));
                switch_index=1/(1+exp(15*switch_index+4));
                if switch_index>=randswitch_k(t,l)
                   types(W+F+l)=0;
                end
            else
                if types(W+F+l)==0
                    if deviation(W+F+l,t)>=0
                        switch_index=0.99*deviation(W+F+l,t)+0.01*(opt_share(t-1)-(pes_share(t-1)+neut_share(t-1)));
                        switch_index=1/(1+exp(15*(-switch_index)+4));
                        if switch_index>=randswitch_k(t,l)
                            types(W+F+l)=1;
                        end
                    else
                        switch_index=0.99*deviation(W+F+l,t)+0.01*(neut_share(t-1)+opt_share(t-1)-pes_share(t-1));
                        switch_index=1/(1+exp(15*switch_index+4));
                        if switch_index>=randswitch_k(t,l)
                            types(W+F+l)=-1;
                        end
                    end
                else
                switch_index=0.99*deviation(W+F+l,t)+0.01*((opt_share(t-1)+neut_share(t-1))-pes_share(t-1));
                switch_index=1/(1+exp(15*(-switch_index)+4));
                if switch_index>=randswitch_k(t,l)
                    types(W+F+l)=0;
                end
                end
            end
        end
        end       
    end

end

if t>=shockperiod && opinions==1
    impulse_rw_opt(t+1)=adapt*impulse_rw_opt(t)+(1-adapt)*(((wb*(1-tax_rate))/price(t+1))/(rwage(t)*impulse_rw_opt(t)));
    impulse_rw_pes(t+1)=adapt*impulse_rw_pes(t)+(1-adapt)*(((wb*(1-tax_rate))/price(t+1))/(rwage(t)*impulse_rw_pes(t)));
    if impulse_rw_opt(t+1)<1 || impulse_rw_opt(t)==1 
        impulse_rw_opt(t+1)=1;
    end
    if impulse_rw_pes(t+1)>1 || impulse_rw_pes(t)==1
        impulse_rw_pes(t+1)=1;
    end
end

wages_t(t) = wb;
types_agg(t)=mean(types);
Type(t,:)=types;
opt_share(t)=sum(types==1)/(W+F+N);
neut_share(t)=sum(types==0)/(W+F+N);
pes_share(t)=sum(types==-1)/(W+F+N);
wincome(t)=sum(workers_income)/price(t);
dincome(t)=sum(dividends_income)/price(t);
dincomek(t)=sum(dividends_income_k)/price(t);
OCC(t)=sum(Occ_status(1,t,:))/W;
%rwage(t)=(wb*(1-tax_rate))/price(t);
DIV(t)=sum(divstatus)/F;
DIVK(t)=sum(divstatus_k)/N;
rDIV(t)=sum(pi)/price(t);
rDIVK(t)=sum(pi_k)/price(t);
step=step+1;
 end

for i=1:(W+F+N)
   cactual(i,:)=cactual(i,:);
   cdemand1(i,:)=cdemand1(i,:);
   cdemand2(i,:)=cdemand2(i,:);
   welfare_cs(i)=mean((cactual(i,1000:1040).^(1-tau)-1)/(1-tau));
   welfare_dc1s(i)=mean((cdemand1(i,1000:1040).^(1-tau)-1)/(1-tau));
   welfare_dc2s(i)=mean((cdemand2(i,1000:1040).^(1-tau)-1)/(1-tau));
   welfare_c(i)=mean((cactual(i,500:T).^(1-tau)-1)/(1-tau));
   welfare_dc1(i)=mean((cdemand1(i,500:T).^(1-tau)-1)/(1-tau));
   welfare_dc2(i)=mean((cdemand2(i,500:T).^(1-tau)-1)/(1-tau));
end
 

et = toc;

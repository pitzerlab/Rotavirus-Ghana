global B  wm wi1 wi2 u um beta b1 b2 phi1 phi2 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; % immig deaths b2 phi2 

load('navrongo_demographic_parameters.mat');            %demographic parameters
load('navrongo_model_parameters.mat');                  %model estimated parameters
% parameters definition
%ghana_age_dist - Ghana age distribution
%ghana_birth_rate - crude birth rate
%ghana_cdr - crude death rate
%ghana_imr  -  immigration
%vcov  -  vaccine coverage 
%Differential equations - -rasisV1w.m


% - -  for results  - - output names- -  - 
 cases_U5Y=zeros(50,120);
 cases_U1Y=zeros(50,120);
 cases_1Y=zeros(50,120);
 cases_2_4Y=zeros(50,120);
 age_dist=zeros(50,3);


%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age);
au2=24; 
au5=27; 
avgage=[1:24 30 42 54];
agegroups={'0-11m','12-23m','24-59m'};
agep=[(1/12)*ghana_age_dist(1)*ones(1,12) (1/48)*ghana_age_dist(2)*ones(1,12) (1/4)*ghana_age_dist(2)*ones(1,3) ghana_age_dist(3:17).*ones(1,15)];


agep=agep/sum(agep); %Proportion of pop'n in each age group (square age distribution)
%u=[1/4.3*ones(1,12) 1/52*ones(1,4) 1/(52*5)*ones(1,14) 1/(52*25)];
u=[1/1.0*ones(1,24) 1/12*ones(1,3) 1/(12*5)*ones(1,14) 1/(12*25)];


%% - -  simulation length - - - - - -
t0=round(12*47);                                    %burn-in period
tvacc=63;                                            %previccination period
vacT=225;                                              %length of vaccination simulation from april 2012  
tmax=t0+tvacc+vacT;                                    %length of simulation

%%-- set up date for the simulations (1960 to 2030)
datepop=[(1960:1:2031)' 12*ones(72,1) ones(72,1) zeros(72,3)];
datesim=(datenum([1960 1 1 0 0 0]):30.45:datenum([2030 12 31 0 0 0]))';  %30.495


 %------------ demographic parameters ------------------ ---------------------
N1=30000;               %population
N=N1*agep;
B=interp1(datenum([1960 1 1 0 0 0; datepop]),[ghana_birth_rate(1,1); ghana_birth_rate(:,1)]/1000,datesim); 
B=[B zeros(tmax,al-1)];
cdr=interp1(datenum(datepop),ghana_cdr/1000,datesim);
cdr=log(1+cdr)/12;
immig=.3*ones(length(cdr),1); 
immig=log(1+immig)/12;
um=.0001*ones(1,al);
um=log(1+um-immig-.0001)/12;  

% %%------end  demographic parameters-----------------------------------------




%%% INFECTION PARAMETERS %%%

dur=1; %duration of infectiousness
d1=1/dur; %rate of recovery from primary infection (per week)
d2=2*d1; %rate of recovery from subsequent infection (per week)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
ri3=0.1; %relative infectiousness of asymptomatic infection
wi1=1/2.999; %rate of waning immunity following primary infection
wi2=1/2.999; %rate of waning immunity following 2nd infection
wA=0; %rate of waning immunity against symptomatic infection in adults



%%% ESTIMATED PARAMETERS %%%

pars=[31.529;0.315;0.9992;1.050;6.63e-08;5.125;0.012;2.54e-05;0.203;5.26]; %mean

for icount=1:50
ptrans=R0_par(icount);                  %transmission parameters ~R0
wm=1/wm_par(icount);                     %rate of waning maternal immunity 
b1=pars(3);                              %amplitude of annual seasonal forcing 
phi1=pars(4);                             %annual seasonal offset
b2=pars(5);                               %amplitude of biannual seasonal forcing   
phi2=pars(6);                             %biannual seasonal offset  
ar1=1; 
ar2=1; 

c2=100*[ar1*ones(al,12) ar2*ones(al,12) ones(al,al-24)]; %Age-related acquisition 

beta=(ptrans/100)*c2; 

immunity=[0.13 0.063; 0.03 0.077];
im=1;

h=pars(7);                   %proportion of severe diarrhea cases hospitalized
hosp1=immunity(1,im)*h*ones(1,al); 
hosp2=immunity(2,im)*h*ones(1,al); 
hosp3=pars(8)*zeros(1,al); 
delta1=0.41*h*ones(1,al); %proportion of primary infections that are symptomatic
delta2=0.35*h*ones(1,al); %proportion of secondary infections that are symptomatic
delta3=0.21*h*ones(1,al); %rate of detection of subsequent infection

reintro=0;

R0=max(eig(dur*beta.*(ones(al,1)*N)));


%%% VACCINATION PARAMETERS
v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al); %initialize vaccination rate across all ages
    
avacc=[3 4 0];                  %age of vaccine adminstration (2 months (1st dose) and 3months (2nd dose), 0 indicates no 3rd dose for 6/10 weeks schedule
sc1=sc1_par(icount);            %proportion that responded to the 1st dose  
sc2=sc2_par(icount);            %proportion that responded to the 2nd dose 
sc3=0*sc2_par(icount);          %proportion that responded to the 3rd dose 

if sc1 < sc2
res=sc1/sc2;                                            %probability of being a "responder"
sc2n=(1-sc2)*sc1/(1-sc1);                               %probability of responding to second dose given did not respond to first dose
sc3n=0*(1-sc2)^2*sc1/(1-res+res*(1-sc2)^2);             %probability of responding to third dose given did not respond to first or second dose
else
res=sc2/sc1;                                            %probability of being a "responder"
sc2n=(1-sc1)*sc2/(1-sc2);                               %probability of responding to second dose given did not respond to first dose
sc3n=0*(1-sc1)^2*sc2/(1-res+res*(1-sc1)^2);             %probability of responding to third dose given did not respond to first or second dose
end

%%% vaccine coverage
v1(:,avacc(1))=[zeros(t0+tvacc,1); vcov(1:end,1)];                %Vaccine coverage by month with 1 dose
if avacc(2)>0
v2(:,avacc(2))=[zeros(t0+tvacc+1,1); vcov(2:end,2)];              %Vaccine coverage by month with 2 doses
end
if avacc(3)>0
v3(:,avacc(3))=0*[zeros(t0+tvacc+2,1); vcov(3:end,2)];              %Vaccine coverage by month with 3 doses
end

wv=1/wv_par(icount);             %rate of waning of vaccine-induced immunity

%Initialize vector to keep track of the number of people in each state
St0=zeros(33*al,1);
St0(1:al,1)=[N(1) zeros(1,al-1)]; %Maternal immunity
St0(al+1:2*al,1)=[0 N(2:al)-[ones(1,al-11) zeros(1,10)]]; %Susceptible_0
St0(2*al+1:3*al,1)=[0 ones(1,al-11) zeros(1,10)]; %Infectious_1 (primary) 
St0(3*al+1:4*al,1)=zeros(1,al); %Recovered_1
St0(4*al+1:5*al,1)=zeros(1,al); %Susceptible_1
St0(5*al+1:6*al,1)=zeros(1,al); %Infectious_2 (2nd time)
St0(6*al+1:7*al,1)=zeros(1,al); %Recovered_2
St0(7*al+1:8*al,1)=zeros(1,al); %Susceptible-Resistant
St0(8*al+1:9*al,1)=zeros(1,al); %Asymptomatic Infectious_3 (subsequent)
St0(9*al+1:10*al,1)=zeros(1,al); %Temp Resistant
St0(10*al+1:11*al,1)=zeros(1,al); %Maternal immunity
St0(11*al+1:12*al,1)=zeros(1,al); %
St0(12*al+1:13*al,1)=zeros(1,al); %
St0(13*al+1:14*al,1)=zeros(1,al); %
St0(14*al+1:15*al,1)=zeros(1,al); %
St0(15*al+1:16*al,1)=zeros(1,al); %SV1
St0(16*al+1:17*al,1)=zeros(1,al); %IV0
St0(17*al+1:18*al,1)=zeros(1,al); %RV0
St0(18*al+1:19*al,1)=zeros(1,al); %SV1
St0(19*al+1:20*al,1)=zeros(1,al); %IV1
St0(20*al+1:21*al,1)=zeros(1,al); %RV1
St0(21*al+1:22*al,1)=zeros(1,al); %SV2
St0(22*al+1:23*al,1)=zeros(1,al); %IV2
St0(23*al+1:24*al,1)=zeros(1,al); %RV2
St0(24*al+1:25*al,1)=zeros(1,al); %MV1
St0(25*al+1:26*al,1)=zeros(1,al); %MV2
St0(26*al+1:27*al,1)=zeros(1,al); %MV3
St0(27*al+1:28*al,1)=zeros(1,al); %MV4
St0(28*al+1:29*al,1)=zeros(1,al); %MV5
St0(29*al+1:30*al,1)=zeros(1,al); %MV6


clear St lambda H

options=odeset('NonNegative',1:length(St0));
[time St]=ode45('rasisV1w',1:tmax,St0,options);


time(1:t0,:)=[];
St(1:t0,:)=[];


lambda=zeros(tmax-t0,al);
for t=1:size(St,1)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1*12)/12)+b2*cos(2*pi*(time(t)-phi2*12)/6))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al)+St(t,16*al+1:17*al)+ri2*St(t,19*al+1:20*al)+ri3*St(t,22*al+1:23*al))*beta)./sum(St(t,:)); %+b2*cos(2*pi*(time(t)-phi2)/26)
end


Cu=zeros(tmax-t0,al); Hu=zeros(tmax-t0,al); Cv=zeros(tmax-t0,al); Hv=zeros(tmax-t0,al); Incid=zeros(tmax-t0,al); Prev=zeros(tmax-t0,al);
for i=1:al
    Cu(:,i)=max(0,delta1(i)*St(:,al+i).*lambda(:,i)+delta2(i)*rr1*St(:,4*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,7*al+i).*lambda(:,i));
    Cv(:,i)=max(0,delta1(i)*St(:,15*al+i).*lambda(:,i)+delta2(i)*rr1*St(:,18*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,21*al+i).*lambda(:,i));
    Hu(:,i)=max(0,hosp1(:,i).*St(:,al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,4*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,7*al+i).*lambda(:,i));
    Hv(:,i)=max(0,hosp1(:,i).*St(:,15*al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,18*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,21*al+i).*lambda(:,i));
end
C=Cu+Cv;            %non-severe cases
H=Hu+Hv;            %moderate-to-severe


%age_distribution
A=(H(:,1:au5)./(sum(H(:,1:au5),2)*ones(1,au5)))*avgage';
A_prevacc=A(1:tvacc,1)'*sum(Hu(1:tvacc,1:au5),2)/sum(sum(Hu(1:tvacc,1:au5)));
A_postvacc=A(64:183,1)'*sum(H(64:183,1:au5),2)/sum(sum(H(64:288,1:au5)));
agedist_prevacc=sum(Hu(1:63,1:au5)./sum(sum(Hu(1:63,1:au5),2)));
agedist_postvacc=sum(H(64:183,1:au5))./sum(sum(H(64:183,1:au5),2));
agedistV=sum(Hu(64:183,1:au5))./sum(sum(Hu(64:183,1:au5),2));
agedistU=sum(Hu(1:tvacc,1:au5))./sum(sum(Hu(1:tvacc,1:au5),2));


% results -- moderate-to-severe  cases 
H=[sum(H(1:288,[1:7]),2) sum(H(1:288,[8:12]),2) sum(H(1:288,[13:24]),2) sum(H(1:288,[25:27]),2)];
Hn=sum(H,2); 
cases_U5Y(icount,1:120) = Hn(64:183);                       %cases u5 years
cases_U1Y(icount,1:120)=H(64:183,1)+H(64:183,2);              %cases under 1 year
cases_1Y(icount,1:120)=H(64:183,3);                           %cases 1Y  
cases_2_4Y(icount,1:120)=H(64:183,4);                         %cases 2-4Y  
age_dist(icount,1:3)=[sum(agedist_postvacc(1:12))  sum(agedist_postvacc(13:24)) sum(agedist_postvacc(25:27))];

%%- -  output file
save('navrongo_results.mat','cases_U5Y', 'cases_U1Y','cases_1Y', 'cases_2_4Y', 'age_dist');

end


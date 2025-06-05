% CCMV model with short selling, risk-neutral interest
% rate, transaction costs and option
% with historical and  forecasting data

clear all
clc
close all

%% data
n=28;
K=15;
rc=0.01;
c=0.1;
xl=-0.1;
xu= 0.1;
M=100;
W0=1000;


% k=0.01* rand(n,1);
% X= xlsread('PTAS2018.xlsx'); % forcasting data
% X(1:2:end, :)= X(1:2:end, :)-1;
% X(2:2:end, :)= X(2:2:end, :)+1;
% % X=X+6;
% Y= xlsread('TAS2018.xlsx') ; % real data

X= xlsread('PTAS2017.xlsx'); % forcasting data

% X(1:2:end, :)= X(1:2:end, :);
% X(2:2:end, :)= X(2:2:end, :);
% X=X+6;
Y= xlsread('TAS2017.xlsx') ; % real data
%  %%%%%%%%%%% 2017 %%%%%%%%%%%%%%%
a=[1 21 40 63 82 104 126 146 169 189 211 232];
b=[20 39 62 81 103 125 145 168 188 210 231 251];

%%%%%%%%%% 2018 %%%%%%%%%%%%%%%
% a=[1 22 41 62 83 105 126 147 170 189 212 233];
% b=[21 40 61 82 104 125 146 169 188 211 232 251];
kc=[ 0.0038
    0.0010
    0.0019
    0.0015
    0.0038
    0.0034
    0.0030
    0.0054
    0.0033
    0.0037
    0.0049
    0.0020
    0.0029
    0.0043
    0.0048
    0.0013
    0.0042
    0.0029
    0.0036
    0.0042
    0.0016
    0.0052
    0.0005
    0.0026
    0.0057
    0.0007
    0.0028
    0.0045];

%% data of assets prices forcasting
% %2017
% X= xlsread('PTA2017.xlsx');
% monthly risk
% a=[1 22 41 62 83 105 126 147 170 189 212 233];
% b=[21 40 61 82 104 125 146 169 188 211 232 251];

%% daily risk
% a1=linspace(1,250,250)
% b1=linspace(2,251,250)
% for i=1:n
%     for j=1:250
% Z(i,j)=(X(b1(j), i)- X(a1(j), i))/ X(a1(j), i);
%  end
% end
% %2018
% X= xlsread('PTA2018.xlsx');
% a=[1 20 39 62 81 103 125 145 168 188 210 231];
% b=[19 38 61 80 102 124 144 167 187 209 230 250];

%%
r=[];
for i=1:n
    r(i)= (X(b(12), i)- X(a(1), i))/ X(a(1), i);
end
mo1=r';
R=[];
%% return
%     for i=1:n
%         for  l=1:12
%             R(i,l)= (X(b(l), i)- X(a(l), i))/ X(a(l), i);
%         end
%     end
%     mo1=mean(R');
%     mo1=mo1';

%% covariance matrix
for i=1:n
    for j=1:12
        Z(i,j)=(X(b(j), i)- X(a(j), i))/ X(a(j), i);
    end
end
C1= cov(Z');

%% real data
r1=[];
for i=1:n
    r1(i)= (Y(b(12), i)- Y(a(1), i))/ Y(a(1), i);
end
mo2=r1';
R=[];

for i=1:n
    for j=1:12
        Z1(i,j)=(Y(b(j), i)- Y(a(j), i))/ Y(a(j), i);
    end
end
C2= cov(Z1');
%%
STR= Y(end,:);
S0=X(1,:); %S0 is first price for eny stock
ST= X(end,:); % ST is last price for eny stock
ST=ST';
S0=S0';
STR=STR';
for i=1:n
    if ST(i) < S0(i) % short
      Strike_price(i)=ST(i);%forcast price(finite price) *...
        option_price(i)=(max(STR(i)- Strike_price(i) , 0))*exp(-rc);%(forcast price - calculated price)*....
        % option_price(i)=(max(STR(i)- Strike_price(i) , 0))*exp(-rc);%(forcast price - calculated price)*....
    else %long
            Strike_price(i)=ST(i);%forcast price(finite price) *...
     option_price(i)=(max( Strike_price(i)- STR(i) , 0))*exp(-rc);%(forcast price - calculated price)*....
        
    end
end

Strike_price= Strike_price';
option_price= option_price';
O= option_price./S0; % option costs
% [Strike_price  option_price./S0];
% [ST  Strike_price];

%% program  without option
lambda= 0.5;
cvx_begin %quiet
cvx_solver mosek
cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME',20)
variables x(n) p(n)
variable z(n) binary
variable h(n) binary
% minimize(lambda*x1'*C1*x1 ...
%   -(1-lambda)*((x1'*mo1)- k'*x1- x1'*O -rc*sum(p)))
minimize(lambda*x'*C2*x ...
    -(1-lambda)*((x'*mo2)- kc'*abs(x)-rc*sum(p)))
subject to
sum(x)<=1;
sum(z)==K;
%short selling constraint
% x >= 0;
%%%%%%%%%%%%%%%
% risk neutral interest rate constraint
p <= M*(h);
c*x - p <= M*(1-h);
p<= c*x;
p>= -(h)*M;
x <= M*(1-h);
x >= -M*(h);
x >= xl*z;
x <= xu*z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  mo2.*x >= 0;
cvx_end
Ret= x'*mo2 - kc'*abs(x)-rc*sum(p);
out.Ret= x'*mo2 - kc'*abs(x)-rc*sum(p);
withoutop= W0*Ret;
out.withoutop= W0*Ret;

%% program  with option
lambda= 0.5;
cvx_begin %quiet
cvx_solver mosek
cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME',20)
variables x1(n) p1(n) %u(n)
variable z1(n) binary
variable h1(n) binary
% minimize(lambda*x1'*C1*x1 ...
%     -(1-lambda)*((x1'*mo1)- kc'*x1- x1'*O + x1'*pp))

minimize(lambda*x1'*C1*x1 ...
    -(1-lambda)*((x1'*mo1)- kc'*abs(x1)- O'*abs(x1)-rc*sum(p1)))
subject to
sum(x1)<=1;
sum(z1)==K;
% u <= x1
% -x1 <= u
x1 >= xl*z1;
x1 <= xu*z1;

%% short selling constraint
% x1 >= 0;

%% risk neutral interest rate constraint
p1 <= M*(h1);
c*x1 - p1 <= M*(1-h1);
p1 <= c*x1;
p1 >= -(h1)*M;
x1 <= M*(1-h1);
x1 >= -M*(h1);
x1 >= xl*z1;
x1 <= xu*z1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mo1.*x1 >= 0;
cvx_end
Retop= x1'*mo1- kc'*abs(x1) -rc*sum(p1)- abs(x1)'*O ;
out.Ret= x'*mo2 - kc'*abs(x)-rc*sum(p);
out.Retop= x1'*mo1- kc'*abs(x1) -rc*sum(p1) - abs(x1)'*O ;
% out.Sharp= ((1/n)*sum((x'*mo2 - kc'*abs(x))))/sqrt(x'*C2*x);
% out.Sharpop= ((1/n)*sum((x1'*mo1 - kc'*abs(x1)-rc*sum(p1) - abs(x1)'*O)))/sqrt(x1'*C1*x1);
out.Sharprc= ((1/n)*sum((x'*mo2 - kc'*abs(x)-rc*sum(p)))-rc)/sqrt(x'*C2*x);
out.Sharpoprc= ((1/n)*sum((x1'*mo1 - kc'*abs(x1) - abs(x1)'*O -rc*sum(p1)))-rc)/sqrt(x1'*C1*x1);

ind=find(x1> 0.0001);
ind1=find(x1< -0.0001);
innd=[];
innd=[ind' ind1'];
innd=innd';
% ZZ=S0(ind);
% SS=sum(ZZ);
% W=ZZ./SS;

%% put option for long position
% ind is for positive index
% DD= max((ST(ind)- Real(ind)),0)

Vput= max((Strike_price(ind)- STR(ind)),0);
weightput= (W0* abs(x1(ind)))./S0(ind);

% Vput= max((Strike_price(innd)- STR(innd)),0);
% weightput= (W0* abs(x1(innd)))./S0(innd);
 
%  Vput= max((Strike_price - STR ),0);
%  weightput= (W0* abs(x1))./S0;
 
Finalput= weightput'*Vput;
withop= W0*Retop;
out.withop= W0*Retop;

%% call option for short position
% ind1 is for negative index
Vcall= max( STR(ind1)- (Strike_price(ind1)),0);
weightcall= (W0* abs(x1(ind1)))./S0(ind1);
% 
% Vcall= max( STR(innd)- (Strike_price(innd)),0);
% weightcall= (W0* abs(x1(innd)))./S0(innd);

% Vcall= max( STR- (Strike_price),0);
% weightcall= (W0* abs(x1))./S0;

Finalcall= weightcall'*Vcall;

%%
TNoption=withop + (Finalcall)+ (Finalput);

out.TNoption=withop + (Finalcall)+ (Finalput);
% 
% out.sharpTN = (((withoutop/W0))-rc)/sqrt(x'*C2*x);
% out.sharpTNop = (((TNoption/W0))-rc)/sqrt(x1'*C1*x1);
Vput./S0(ind)
Vcall./S0(ind1)
% TNoption=withop- (weightput'*(option_price(ind))')...
%     - (weightcall'*(option_price(ind1))')...
%     + (Finalcall)+ (Finalput);
% out.TNoption=withop- (weightput'*(option_price(ind))')...
%     - (weightcall'*(option_price(ind1))')...
%     + (Finalcall)+ (Finalput);
out









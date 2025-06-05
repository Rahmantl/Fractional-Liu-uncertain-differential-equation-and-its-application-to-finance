% Paper title :
%
clc; clear all; close all;
warning off;
%% Real data     Load Dataset
Y = xlsread('2019realdata.xlsx') ;     % real data
%                Number of asset
n=23;
%
a=[1 22 43 64 85 106 127 148 169 190 211 232];
b=[21 42 63 84 105 126 147 168 189 210 231 252];
%%    Real data         Return
for i=1:n
    r(i)= (Y(b(12), i)- Y(a(1), i))/ Y(a(1), i);
end

%% protfolio with 15percent Real data
ninterval=5;
[col_ret_varY,Ynn,rYnn,nnn] = portfolio15percent(n,r,Y,ninterval);

%%    Real data            Covariance matrix
for i=1:nnn
    for j=1:12
        Z(i,j) = (Ynn(b(j), i) - Ynn(a(j), i)) / Ynn(a(j), i);
    end
end
Q= cov(Z');
var_Q=var(Q); %var_return

%% Set input parameters
%                The lower and upper bounds of the asset
xlow =0;
xup =0.1;
% cardinality parameter
K = 10;
lambda= 0.5; % Lambda in [0,1]
%%  Real data               CVX implemention
% if K>nnn
%     error('Cardinality parameter must less or equal than number of asset (K<=n)')
% end

 [x_real, zi] = solve_prob(lambda, Q, rYnn, K, xlow, xup, nnn);
 
 %format long
Column_number_return_variance_Realdata=[col_ret_varY var_Q']
%format short
%% *****************************
%% Predict               Load Dataset
X =   xlsread('2019predict.xlsx');  % predict
%%    Predict         Return
for i=1:n
    rp(i)= (X(b(12), i)- X(a(1), i))/ X(a(1), i);
end

%%
% % % % % % % index_X=1:n;
% % % % % % % rrp = [index_X; r];rrp=rrp';
% % % % % % % rrp=sortrows(rr,2);
% % % % % % % 
% % % % % % % [m,~] = size(rrp);
% % % % % % % m_per = floor(0.15*m);
% % % % % % % m1 = m_per+1; m_end = m-m_per;
% % % % % % % r_newp = rr(m1:m_end,:);
% % % % % % % 
% % % % % % % Xindex = r_newp(:,1);
% % % % % % % rnp=r_newp(:,2);
% % % % % % % Xn = X(:,Xindex);

[col_ret_varX,Xnn,rXnn,nnX] = portfolio15percent(n,rp,X,ninterval);
[~,nnp] = size(Xnn);

%%   Predict            Covariance matrix
for i=1:nnX
    for j=1:12
        Zp(i,j) = (Xnn(b(j), i) - Xnn(a(j), i)) / Xnn(a(j), i);
    end
end
Qp= cov(Zp');

% col_ret_varPrice_varRet = [col_ret_varPrice_varRet] 
% Column_number_return_variance_predict=col_ret_varPrice_varRet
%%  Predict               CVX implemention
lambda= 0.5; % Lambda in [0,1]
[x_pridict, zip] = solve_prob(lambda, Qp, rXnn, K, xlow, xup, nnX);


Column_number_return_variance_Realdata=[col_ret_varY var_Q']
Column_number_return_variance_predict=[col_ret_varX var_Q']

% 
% % 
% % 
% % %% Write answer in excell
% % 
% % %ans = []
% % 
% % %xlswrite('output.xls',ans)

function [valid_col_ret_var_final,Ynn,rYnn,nnn] = portfolio15percent(n,r,Y,n_interval)

%% remove %15 of lower and upper elements of returns
index_Y=1:n;
rr1 = [index_Y; r]';
rr=sortrows(rr1,2);

[m,~] = size(rr);
m_per = floor(0.05*m);   %%  percent
m_up = floor(0.10*m);   %%  percent
m1 = m_per+1; m_end = m-m_up;
r_new = rr(m1:m_end,:);

Yindex = r_new(:,1);
rn=r_new(:,2);
Yn = Y(:,Yindex);
[~,nn] = size(Yn);

%% divided return interval to 5 subinteval
%ninterval=5;
Lspace =linspace(rn(1),rn(end),n_interval+1);
Lspace(end) = Lspace(end)+0.002;
%% select one return in each subinterval with lower variance
elemnt=2;%%we need this elmnt in each interval
VarianceYn=var(Yn);
Lsize = size(Lspace);


i=1;
val_in_interval = zeros(n_interval*elemnt,3);
for indd=Lsize(end) :-1: Lsize(1)+1
    value_sign = (r_new(:,2) < Lspace(indd) & r_new(:,2) >= Lspace(indd-1));
    %compYn = value_sign .* r_new(:,2);
    indr = value_sign .* r_new(:,1);ss =size(indr);
    r_var_real = [indr'; r_new(:,2)'; VarianceYn]';
    r_var = ([indr'; r_new(:,2)'; VarianceYn]' .* [ones(ss) value_sign value_sign]);
    r_var( r_var == 0 ) = NaN;
    r_var = sortrows( r_var, 3 );
    r_var( isnan( r_var ) ) = 0;
     
     if nnz(r_var(:,1)) > elemnt
        nnz_each_interval(i) = 1;%2 means that mor than elemnt we have in this interval 
        %nn_element = nnz(r_var(:,1));
        i=i+1;
        %val_in_interval(j:j+,:) = r_var(1:elemnt,:);
        val_in_interval(indd*elemnt-2*elemnt+1:indd*elemnt-elemnt,:) = r_var(1:elemnt,:);
         
     elseif nnz(r_var(:,1)) == elemnt
        nnz_each_interval(i) = 0;%0 means that equal to elemnt we have in this interval 
        %nn_element = nnz(r_var(:,1));
        i=i+1;
         val_in_interval(indd*elemnt-2*elemnt+1:indd*elemnt-elemnt,:) = r_var(1:elemnt,:);
         % r_var_real(r_var_real(1:elemnt,:)) = NaN;
    elseif nnz(r_var(:,1)) < elemnt
        nnz_each_interval(i) = -1;%2 means that less than elemnt we have in this interval
        i=i+1;
        val_in_interval(indd*elemnt-2*elemnt+1:indd*elemnt-elemnt,:)  = r_var(1:elemnt,:);
        %remain=elemnt - nnz(r_var(:,1));
     end

end

col_ret_var=val_in_interval(1:end,1:end);

n_zero_element=numel(col_ret_var(:,1))-nnz(col_ret_var(:,1));


if n_zero_element>=1
    
valid_col_ret_var = col_ret_var;
valid_col_ret_var = col_ret_var;
valid_col_ret_var ( valid_col_ret_var  == 0 ) = NaN;
valid_col_ret_var  = sortrows( valid_col_ret_var , 2 );

valid_col_ret_var = valid_col_ret_var(1:end-n_zero_element,:);
we_need=n_zero_element

rr1(valid_col_ret_var(:,1),1:2) = zeros(numel(valid_col_ret_var(:,1)),2);
%rr1( rr1 == 0 ) = NaN;
rr1 = sortrows( rr1, 2 ,'descend');
%rr1( isnan( r_var ) ) = 0;
valid_col_ret_var_final = [valid_col_ret_var(:,1:2);rr1(n_zero_element,:)];
else 
    valid_col_ret_var_final = col_ret_var;
end


%% find companies column from data set where previous conditions are passed
Ynn=Y(:,valid_col_ret_var_final(:,1));
rYnn = valid_col_ret_var_final(:,2);
[~,nnn] = size(Ynn);

end


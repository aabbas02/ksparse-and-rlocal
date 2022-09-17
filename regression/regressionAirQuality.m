clc 
close all 
clear all
rng('default')
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
cd(newdir)
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altMinProposed'),...
        genpath('.\dataSets'))
A = readmatrix('air_quality_data.csv');
% retain rows of data from years 
idx1 = find(A(:,2) == 2013, 1 );
idx2 = find(A(:,2) == 2017, 1, 'last' );
A = A(idx1:idx2,:);
idx = 1:size(A,2);
% delete columsns 16, 18 two columns comprising NANs
idx = setdiff(idx,[size(A,2),size(A,2)-2]);
A = A(:,idx);
[row, col] = find(isnan(A));
A = A(setdiff(1:size(A,1),row),:);
% preprocess - remove outliers
Y = A(:,[6,7,8,9,11]);
[Y,TF] = rmoutliers(Y,'movmedian',128);
size(Y,1)
A = A(~TF,:);
Y = sqrt(Y);
X = zeros(size(A,1),27);
% WPSM is column 16, instead of column 17; column 16 delted above
X(:,1:6) = A(:,[10,12,13,14,15,16]); 
X(:,7:12) = X(:,1:6).^2;
t=13;
for i = 1 : 6
    for j = i+1:6
        X(:,t) = X(:,i).*X(:,j);
        t=t+1;
    end
end
X = X - mean(X,1);
Y = Y - mean(Y,1);

[U,S,V] = svd(X,'econ');
% improve conditioning
X = U;
% round temperature (col 12), air pressure (col 13) to nearest integers
A(:,12) = round(A(:,12));
A(:,13) = round(A(:,13));
% get block label
% A(:,2,3,4,5) - year,month,day,hour
%--------------------------------------
blkLabel = A(:,3) + 1e4*A(:,4);  %0.65,0.65,0.64,0.65
%blkLabel = A(:,2) + 1e4*A(:,3); 0.48,0.46,0.40
blkLabel = A(:,2) + 1e4*A(:,4); %0.62,0.60,0.44,0.60
%blkLabel = A(:,2) + 1e4*A(:,4);
%blkLabel = A(:,4) + 1e4*A(:,5); %0.62,0.59,0.17,0.59
%blkLabel = A(:,2) + 1e3*A(:,4) + 1e5*A(:,5);

%--------------------------------------
%blkLabel = A(:,2) + 1e5*A(:,5); % collapsed init better
%blkLabel = A(:,3) + 1e5*A(:,5); % LS init better '16 - '17, across years same
%blkLabel = A(:,3) + 1e5*A(:,4); % same
%blkLabel = A(:,2) + 1e5*A(:,3); % LS init much better
[blkLabelSorted,idx] = sort(blkLabel);
% order blockwise
Y = Y(idx,:);
X = X(idx,:);
% get lengths of blocks
temp = unique(blkLabelSorted);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blkLabelSorted == temp(i),1,'first');
    t2 = find(blkLabelSorted == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
n   = size(Y,1);
pi_ = get_permutation_r(n,r_); 
%pi_ = get_permutation_k(n,round(n/2));
rLocal = 1;
Y_permuted = Y(pi_,:);
%---------------- oracle -----------------------------------
Btrue = X \ Y;
R2_true  = 1 - norm(Y-X*Btrue,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%---------------- naive ------------------------------------
Bnaive = X \ Y_permuted;
R2_naive  = 1 - norm(Y-X*Bnaive,'fro')^2/norm(Y,'fro')^2
%---------------- proposed ----------------------------------
maxIter = 20;
lsInit = 0;
%---------------- w collapsed init --------------------------
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
BproErr  = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro   = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
%--------------- w least-squares init -----------------------
lsInit = 1;
[pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro = X(pi_hat,:) \ Y_permuted;
R2_proLS  = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
%------------------ slawski ---------------------------------
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:) \ Y_permuted;
beta_sls_err = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
%----------------- RLUS ---------------------------------------
tic
[pi_hat] = rlus(X,Y_permuted,r_,rLocal);
beta_RLUS = X(pi_hat,:) \ Y_permuted;
R2_rlus  = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
beta_rlus_err = norm(beta_RLUS - Btrue,2)/norm(Btrue,2);
tRlus = toc;
%----------------------------------------------------------------
%pi_icml     = icml_20(X,Y_permuted,r_);
%beta_icml   = X(pi_icml,:) \ Y_permuted;
%R2_icml     = 1 - norm(Y-X*beta_icml,'fro')^2/norm(Y,'fro')^2;
%BicmlErr    = norm(beta_icml - Btrue,2)/norm(Btrue,2); 

num_blocks = length(r_)
R2_true 
R2_naive
R2_pro
R2_proLS
R2_rlus
R2_sls
%R2_icml
%fVal
%R2_proLS
%fValLS
%R2_rlus
BproLSerr
BproErr
beta_sls_err
beta_rlus_err
%BicmlErr

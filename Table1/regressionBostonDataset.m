% dataset and loading code at https://www.mathworks.com/products/demos/machine-learning/boosted-regression.html
close all 
clear all
clc
rng('default')
dir = pwd;
% For linux, replace '\' with '/'
idcs   = strfind(dir,'\');
newdir = dir(1:idcs(end)-1);
cd (newdir)
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altMinProposed'),...
        genpath('.\dataSets'))

cd(dir)
[X,y] = loadBostonHousing();
R2Oracle = 1 - norm(y - X*pinv(X)*y,'fro')^2/norm(y,'fro')^2;
% zero-mean the columns
X = X - mean(X,1);
y = y - mean(y,1);% make one of the features the block label and round that feature to sf
% significant figures
sf = 1;
col = 1;
numBlocks = length(unique(round(X(:,col),sf)));
blkLabel = round(X(:,col),sf);
[blkLabelSorted,idx] = sort(blkLabel);
% order blockwise
y = y(idx,:);
X = X(idx,:);
% get lengths of blocks
temp = unique(blkLabelSorted);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blkLabelSorted == temp(i),1,'first');
    t2 = find(blkLabelSorted == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
% permute the data
n = size(X,1);
pi_ = get_permutation_r(n,r_);
rLocal = 1;
y_permuted = y(pi_,:);
%------------ oracle ---------------------------------------------------
Btrue = X\y;
yhat = X*Btrue;
R2_true =  1 - norm(y-yhat,'fro')^2/norm(y,'fro')^2;
%----------- naive -----------------------------------------------------
Bnaive = X\y_permuted;
yhat = X*Bnaive;
R2_naive =  1 - norm(y-yhat,'fro')^2/norm(y,'fro')^2;
beta_naive_err = norm(Bnaive - Btrue,2)/norm(Btrue,2);
%----------- proposed ----------------------------------
maxIter = 50;
lsInit = 0;
%---------- w collapsed init --------------------------
[pi_hat,fVal] = AltMin(X,y_permuted,r_,maxIter,rLocal,lsInit);
Bpro    = X(pi_hat,:) \ y_permuted;
beta_pro_err = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro       = 1 - norm(y-X*Bpro,'fro')^2/norm(y,'fro')^2;
%---------- w least-squares init -----------------------
lsInit       = 1;
[pi_hat,fValLS]   = AltMin(X,y_permuted,r_,maxIter,rLocal,lsInit);
Bpro         = X(pi_hat,:) \ y_permuted;
R2_proLS     = 1 - norm(y-X*Bpro,'fro')^2/norm(y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
%------------------ slawski ---------------------------------
noise_var    = norm(y_permuted-X*Bnaive,'fro')^2/(size(y,1)*size(y,2));
tic
[pi_hat,~]   = slawski(X,y_permuted,noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:)\y_permuted;
beta_sls_err = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(y-X*beta_sls,'fro')^2/norm(y,'fro')^2;
%----------------- RLUS ---------------------------------------
tic
[pi_hat] = rlus(X,y_permuted,r_,rLocal);
beta_RLUS = X(pi_hat,:) \ y_permuted;
R2_rlus  = 1 - norm(y-X*beta_RLUS,'fro')^2/norm(y,'fro')^2;
tRlus = toc;
beta_rlus_err = norm(beta_RLUS - Btrue,2)/norm(Btrue,2);

R2_true
R2_pro
R2_sls
R2_rlus
R2_proLS
function [X,y] = loadBostonHousing()
    filename = 'housing.txt';
    urlwrite('http://archive.ics.uci.edu/ml/machine-learning-databases/housing/housing.data',filename);
    inputNames = {'CRIM','ZN','INDUS','CHAS','NOX','RM','AGE','DIS','RAD','TAX','PTRATIO','B','LSTAT'};
    outputNames = {'MEDV'};
    housingAttributes = [inputNames,outputNames];

    formatSpec = '%8f%7f%8f%3f%8f%8f%7f%8f%4f%7f%7f%7f%7f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
    fclose(fileID);
    housing = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','VarName2','VarName3','VarName4','VarName5',...
        'VarName6','VarName7','VarName8','VarName9', ...
        'VarName10','VarName11','VarName12','VarName13','VarName14'});

% Delete the file and clear temporary variables
    clearvars filename formatSpec fileID dataArray ans;
    delete housing.txt

    housing.Properties.VariableNames = housingAttributes;
    X = housing{:,inputNames};  
    y = housing{:,outputNames};
end
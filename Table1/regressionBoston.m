clc 
close all 
clear all
rng('default')
dir = pwd;
cd .. 
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altGDMin'),...
        genpath('.\altMinProposed'),...
        genpath('.\dataSets'))
cd(dir)
% load Boston Housing Data set
[X, Y] = loadAndProcessHousingData();
% center the response and predictor variables
X = X - mean(X,1);
Y = Y - mean(Y,1);
% make one of the features the block label and round that feature to sf
% significant figures
sf = 0;
col = size(X,2) - 1;
numBlocks = length(unique(round(X(:,col),sf))); 
blkLabel = round(X(:,col),sf);
% sort blockwise
[blkLabelSorted,idx]  = sort(blkLabel);
X = X(idx,:);
Y = Y(idx,:);
% permute within block
temp = unique(blkLabel);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blkLabelSorted == temp(i),1,'first');
    t2 = find(blkLabelSorted == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
n   = size(Y,1);
pi_ = get_permutation_r(n,r_);
Y_permuted = Y(pi_,:);
% replace by SVD - often helps the collapsed initialization
[U,S,V] = svd(X,'econ');
X = U;
%X = U(:,1:size(X,2));
%------------ oracle ---------------------------------------------------
Btrue = X\Y;
Yhat = X*Btrue;
R2_true =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2;
%----------- naive -----------------------------------------------------
Bnaive = X\Y_permuted;
Yhat = X*Bnaive;
R2_naive =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2
beta_naive_err = norm(Bnaive - Btrue,2)/norm(Btrue,2);
%----------- proposed AltMin ----------------------------------
maxIter = 25;
%---------- w collapsed init --------------------------
tic
rLocal = 1;
lsInit = 0;
[pi_hat,~] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tAltMin = toc;
Bpro    = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro       = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
% %---------- w least-squares init -----------------------
lsInit       = 1;
[pi_hat,~] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro         = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
%---------- AltGDMin w collapsed init --------------------------
tic
rLocal = 1;
lsInit = 0;
[pi_hat,~] = altGDMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tAltMin = toc;
B_altGDMin = X(pi_hat,:) \ Y_permuted;
beta_pro_altGDMin_err = norm(B_altGDMin - Btrue,2)/norm(Btrue,2);
R2_altGDMin = 1 - norm(Y-X*B_altGDMin,'fro')^2/norm(Y,'fro')^2;
% %---------- AltGDMin w least-squares init -----------------------
lsInit = 1;
[pi_hat,~] = altGDMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
B_altGDMin         = X(pi_hat,:) \ Y_permuted;
R2_altGDMinLS     = 1 - norm(Y-X*B_altGDMin,'fro')^2/norm(Y,'fro')^2;
BproAltGDMinLSerr = norm(B_altGDMin - Btrue,2)/norm(Btrue,2);
%------------------ slawski ---------------------------------
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:)\Y_permuted;
beta_sls_err = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
%----------------- RLUS ---------------------------------------
tic
[pi_hat] = rlus(X,Y_permuted,r_,rLocal);
beta_RLUS = X(pi_hat,:) \ Y_permuted;
R2_rlus  = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
tRlus = toc;
beta_rlus_err = norm(beta_RLUS - Btrue,2)/norm(Btrue,2);
%---------------------------------------------------------------
R2_true 
R2_naive
R2_pro
R2_proLS
R2_sls
R2_rlus
R2_altGDMin
R2_altGDMinLS

beta_naive_err
beta_pro_err
BproLSerr
beta_sls_err
beta_rlus_err
beta_pro_altGDMin_err
BproAltGDMinLSerr
cd(dir)

function [X, Y] = loadAndProcessHousingData()
    % This function is taken entirely from the mathworks website
    % loadAndProcessHousingData Downloads and processes the housing dataset
    %
    % Input:
    %   url - URL of the dataset
    %   filename - Name of the file to save the downloaded data
    %
    % Output:
    %   X - Input attributes of the dataset
    %   Y - Output attribute of the dataset

    % Download the dataset from the provided URL
    filename = 'housing.txt';
    urlwrite('http://archive.ics.uci.edu/ml/machine-learning-databases/housing/housing.data',filename);
    % Define column names
    inputNames = {'CRIM','ZN','INDUS','CHAS','NOX','RM','AGE','DIS','RAD','TAX','PTRATIO','B','LSTAT'};
    outputNames = {'MEDV'};
    housingAttributes = [inputNames, outputNames];

    % Specify the format of the data to be read
    formatSpec = '%8f%7f%8f%3f%8f%8f%7f%8f%4f%7f%7f%7f%7f%f%[^\n\r]';
    
    % Open the file for reading
    fileID = fopen(filename, 'r');
    % Read data using the specified format
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
    % Close the file
    fclose(fileID);

    % Create a table from the data array
    housing = table(dataArray{1:end-1}, 'VariableNames', housingAttributes);
    
    % Extract input (X) and output (Y) data
    X = housing{:,inputNames};
    Y = housing{:,outputNames};

    % Clean up: Delete the file and clear temporary variables
    delete(filename);
    clearvars filename url formatSpec fileID dataArray;
end


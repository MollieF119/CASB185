%% Data Structure

% MEG/Restin/rmegpreproc/

% data (1x1 struct)
% data.trial -- MEG sig segmented into 2 second intervals, need to merge into one vector
% data.time -- 2 second time intervals, not very informative
% data.label -- label of each channel (bad channels removed)

%% Loading Data

data1 = load("Data/100307_MEG_5-Restin_rmegpreproc.mat");
data2 = load("Data/");
data3 = load("Data/");
data4 = load("Data/");
data5 = load("Data/");
data6 = load("Data/");
data7 = load("Data/");
data8 = load("Data/");
data9 = load("Data/");
data10 = load("Data/");

DATA = {data1, data2, data3, data4, data5, data6, data7, data8, data9, data10};

%% Variables

% sampling frequency
% Fs = data.fsample; 
Fs = 508.6275; % should be same for all data files

% each segment: 2 sec, 1018 samples
% seg = size(data.trial,2); 
seg = 100; % take the first 100 segments

% time vector (unit: seconds)
time = linspace(0,2*seg,1018*seg);
lenT = length(time);

% threshold of correlation
thres = 0.01;

% parameter p
p = 0:0.0001:1.5;
lenP = length(p);

%% Test Matrix

% Fs = 508.63;
% N = 50;
% time = linspace(0,2*50,1018*50);
% 
% MEG_raw = rand(N, length(time));

%% Data Processing

L2_Het_L_Tab = zeros(10,lenP);
L2_Het_P_Tab = zeros(10,lenP);
L2_Hom_L_Tab = zeros(10,lenP);
L2_Hom_P_Tab = zeros(10,lenP);

STD_Tab = zeros(10,4);
PME_Tab = zeros(10,4);

for i = 1:10

    data = DATA{i}.data; % for each individual:

    % number of channels
    N = size(data.label,1);

    % constructing meg_raw matrix 
    meg_raw = megraw(data, N, lenT, seg);

    % filtering meg_raw into alpha and beta frequencies --> meg_filt
    meg_filt = filt(meg_raw, N, lenT, Fs);

    % obtaining base adjacency matrix & remove below-threshold edges
    A = adj(meg_filt, N, thres);

    % constructing heterogeneous/homogeneous multilayer/multiplex networks
    [Ma, Mb, Pab_Het_L, Pba_Het_L, Pab_Het_P, Pba_Het_P, W, Pab_Hom_L, Pab_Hom_P] = prep(A,N);

    % calculating Lambda2 for each network
    L2_Het_L = calc(N, p, len, Ma, Mb, Pab_Het_L, Pba_Het_L);
    L2_Het_P = calc(N, p, len, Ma, Mb, Pab_Het_P, Pba_Het_P);
    L2_Hom_L = calc(N, p, len, Ma, Mb, Pab_Hom_L, Pab_Hom_L);
    L2_Hom_P = calc(N, p, len, Ma, Mb, Pab_Hom_P, Pab_Hom_P);

    L2_Het_L_Tab(i,:) = L2_Het_L;
    L2_Het_P_Tab(i,:) = L2_Het_P;
    L2_Hom_L_Tab(i,:) = L2_Hom_L;
    L2_Hom_P_Tab(i,:) = L2_Hom_P;
    
    % calculating std of interlayer edges for each network
    STD = interstd(Pab_Het_L, Pab_Het_P, Pab_Hom_L, Pab_Hom_P);
    STD_Tab(i,:) = STD;

    % calculating percent of missing interlayer edges for each network
    PME = pme(Pab_Het_L, Pab_Het_P, Pab_Hom_L, Pab_Hom_P);
    PME_Tab(i,:) = PME;

end

%% Average over 10 trials

AvgL2_Het_L = mean(L2_Het_L_Tab);
AvgL2_Het_P = mean(L2_Het_P_Tab);
AvgL2_Hom_L = mean(L2_Hom_L_Tab);
AvgL2_Hom_P = mean(L2_Hom_P_Tab);

AvgSTD = mean(STD_Tab);
AvgPME = mean(PME_Tab);

%% Plotting the graph

% X-axis: p; Y-axis: Lambda2
plot(p, AvgL2_Het_L, "bo");
hold on
plot(p, AvgL2_Het_P, "ro");
plot(p, AvgL2_Hom_L, Color="#EDB120", Marker="o",LineStyle="none");
plot(p, AvgL2_Hom_P, "ko");
hold off

legend('Multilayer', 'Multiplex', 'Homogeneous Multilayer', 'Homogeneous Multiplex');
xlabel('p');
ylabel('\lambda_2');


%% Function: MEG raw matrix

function meg_raw = megraw(data, N, lenT, seg)

    meg_raw = zeros(N,lenT); 
    
    for i = 1:seg
        meg_raw(:,(i-1)*1018+1 : i*1018) = data.trial{i};
    end

end

%% Function: Filtering (alpha & beta)

function meg_filt = filt(meg_raw, N, lenT, Fs)

    % MEG filtered signal matrix
    meg_filt = zeros(2*N, lenT);
    
    % Alpha bandpass filtering
    alpha_low = 8;
    alpha_high = 12;
    [b, a] = butter(2, [alpha_low, alpha_high]*2/Fs, 'bandpass');
    
    % Filling up rows 1 to N (alpha)
    for s = 1:N
        filteredAlpha = filtfilt(b, a, meg_raw(s,:)); %--> s-th row of the matrix
        meg_filt(s,:) = filteredAlpha;
    end
    
    % Beta bandpass filtering
    beta_low = 12;
    beta_high = 30;
    [b, a] = butter(2, [beta_low, beta_high]*2/Fs, 'bandpass');
    
    % Filling up rows N+1 to 2N (beta)
    for s = 1:N
        filteredBeta = filtfilt(b, a, meg_raw(s,:)); %--> (s+N)th row of the matrix
        meg_filt(s+N,:) = filteredBeta;
    end

end

%% Function: Base Adjacency Matrix (A) + Remove below-threshold edges

function adj = adj(meg_filt, N, thres)

    adj = zeros(2*N);
        
    for i = 1:2*N
        for j = 1:2*N
            C = corrcoef(meg_filt(i,:),meg_filt(j,:));
            if C > thres
                adj(i,j) = C(1,2);
            else
                adj(i,j) = 0;
            end
        end
    end
    
end

%% Function: Prepare the submatrices of each of the four networks

function [Ma, Mb, Pab_Het_L, Pba_Het_L, Pab_Het_P, Pba_Het_P, W, Pab_Hom_L, Pab_Hom_P] = prep(A,N)

    % 4 sub-matrices
    Ma = A(1:N, 1:N);
    Mb = A(N+1:2*N, N+1:2*N);
    Pab_Het_L = A(1:N, N+1:2*N);
    Pba_Het_L = A(N+1:2*N, 1:N);
        
    % Het Plex
    Pab_Het_P = Pab_Het_L .* eye(N); % element-wise multiplication with identity matrix
    Pba_Het_P = Pba_Het_L .* eye(N); % only diagonals of interlayer matrices remain
    
    % Mean Interlayer Edge Weight
    W = mean(Pab_Het_L, "all"); % --> for homogeneous networks

    % Hom Layer
    Pab_Hom_L = W * ones(N);
    
    % Hom Plex
    Pab_Hom_P = W * eye(N);

end

%% Function: calculate Lambda 2 for each network

function L2 = calc(N, p, len, Ma, Mb, Pab, Pba)
    
    L2 = ones(1,len);

    for k = 1:len
        
        % Adjacency Matrix
        A = [Ma, p(k).*Pab; p(k).*Pba, Mb];
        
        % Degree Matrix
        D = zeros(2*N);
        
        for i = 1:2*N
            deg = nnz(A(i,:));
            D(i,i) = deg;
        end
        
        % Laplacian Matrix
        L = D - A;
        
        % Lambda 2
        lambda2 = max(eigs(L,2,'smallestabs'));
        L2(k) = lambda2;
        
    end
end

%% Function: Calculate std of interlayer edge weight for each network

function stdTab = interstd(Pab_Het_L, Pab_Het_P, Pab_Hom_L, Pab_Hom_P)

    Pab_Het_L2 = Pab_Het_L;
    Pab_Het_P2 = Pab_Het_P;
    Pab_Hom_L2 = Pab_Hom_L;
    Pab_Hom_P2 = Pab_Hom_P;
    
    Pab_Het_L2(Pab_Het_L2 == 0) = NaN;
    Pab_Het_P2(Pab_Het_P2 == 0) = NaN;
    Pab_Hom_L2(Pab_Hom_L2 == 0) = NaN;
    Pab_Hom_P2(Pab_Hom_P2 == 0) = NaN;
    
    std_Het_L = std(Pab_Het_L2, [], "all", "omitnan");
    std_Het_P = std(Pab_Het_P2, [], "all", "omitnan");
    std_Hom_L = std(Pab_Hom_L2, [], "all", "omitnan"); % should be 0
    std_Hom_P = std(Pab_Hom_P2, [], "all", "omitnan"); % should be 0
    
    stdTab = [std_Het_L std_Het_P std_Hom_L std_Hom_P]; % tabulate

end

%% Function: Calculate % of missing edges for each network

function pmeTab = pme(Pab_Het_L, Pab_Het_P, Pab_Hom_L, Pab_Hom_P)

    pme_Het_L = (numel(Pab_Het_L)-nnz(Pab_Het_L))/numel(Pab_Het_L);
    pme_Het_P = (numel(Pab_Het_P)-nnz(Pab_Het_P))/numel(Pab_Het_P);
    pme_Hom_L = (numel(Pab_Hom_L)-nnz(Pab_Hom_L))/numel(Pab_Hom_L);
    pme_Hom_P = (numel(Pab_Hom_P)-nnz(Pab_Hom_P))/numel(Pab_Hom_P);
    
    pmeTab = [pme_Het_L pme_Het_P pme_Hom_L pme_Hom_P]; % tabulate

end

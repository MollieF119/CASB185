%% Data Structure

% MEG/Restin/rmegpreproc/

% data (1x1 struct)
% data.trial -- MEG sig segmented into 2 second intervals, need to merge into one vector
% data.time -- 2 second time intervals, not very informative
% data.label -- label of each channel (bad channels removed)

% one way to simplify: select only one trial per person, and only select trials w/o bad channels

%% Raw MEG Matrix

Fs = data.fsample; % sampling frequency
N = size(data.label,1); % number of channels
seg = size(data.trial,2); % each segment: 2 sec, 1018 samples
time = linspace(0,2*seg,1018*seg); % time vector (unit: seconds)

MEG_raw = zeros(N,length(time)); % MEG raw signal matrix
for i = 1:seg
    MEG_raw(:,(i-1)*1018+1 : i*1018) = data.trial{i};
end

%% Test Matrix

Fs = 508.63;
N = 50;
time = linspace(0,2*50,1018*50);

MEG_raw = rand(N, length(time));

%% Filter alpha waves

% Alpha bandpass filtering
alpha_low = 8;
alpha_high = 12;
[b, a] = butter(2, [alpha_low, alpha_high]*2/Fs, 'bandpass');

% MEG filtered signal matrix
MEG_filtered = zeros(2*N,length(time)); 

% Filling up rows 1 to N (alpha)
for s = 1:N
    filteredAlpha = filtfilt(b, a, MEG_raw(s,:)); %--> s-th row of the matrix
    MEG_filtered(s,:) = filteredAlpha;
end

%% Filter beta waves

% Beta bandpass filtering
beta_low = 12;
beta_high = 30;
[b, a] = butter(2, [beta_low, beta_high]*2/Fs, 'bandpass');

% Filling up rows N+1 to 2N (beta)
for s = 1:N
    filteredBeta = filtfilt(b, a, MEG_raw(s,:)); %--> (s+N)th row of the matrix
    MEG_filtered(s+N,:) = filteredBeta;
end

%% Plot for the k-th channel

% k = N-2;
% 
% figure
% plot(time,MEG_raw(k,:));
% hold on
% plot(time,MEG_filtered(k,:));
% plot(time,MEG_filtered(k+N,:));
% hold off
% legend("original","alpha","beta");

%% Base Adjacency Matrix (A)
% Also remove non-sifnificant edges

A = zeros(2*N);

thres = 0.01; % threshold of correlation coefficient

for i = 1:2*N
    for j = 1:2*N
        C = corrcoef(MEG_filtered(i,:),MEG_filtered(j,:));
        if C > thres
            A(i,j) = C(1,2);
        else
            A(i,j) = 0;
        end
    end
end

%% Prepare the submatrices of each of the four networks

% 4 sub-matrices
Ma = A(1:N, 1:N);
Mb = A(N+1:2*N, N+1:2*N);
Pab_Het_L = A(1:N, N+1:2*N);
Pba_Het_L = A(N+1:2*N, 1:N);

% Mean Interlayer Edge Weight
W = mean(Pab_Het_L, "all"); % --> for homogeneous networks

% Het Plex
Pab_Het_P = Pab_Het_L .* eye(N); % element-wise multiplication with identity matrix
Pba_Het_P = Pba_Het_L .* eye(N); % only diagonals of interlayer matrices remain

% Hom Layer
Pab_Hom_L = W * ones(N);

% Hom Plex
Pab_Hom_P = W * eye(N);

% parameter p
p = 0:0.0001:1.5;
len = length(p);

%% Calculate the lambda2 and Sp for each network

% A -- Adjacency Matrix
% D -- Degree Matrix
% L -- Laplacian Matrix
    % L = D - A;

% Heterogeneous Multi-Layer Network (Het_L)
[L2_Het_L, Sp_Het_L] = calc(N, p, Ma, Mb, Pab_Het_L, Pba_Het_L);

% Heterogeneous Multi-Plex Network (Het_P)
[L2_Het_P, Sp_Het_P] = calc(N, p, Ma, Mb, Pab_Het_P, Pba_Het_P);

% Homogeneous Multi-Layer Network (Hom_L)
[L2_Hom_L, Sp_Hom_L] = calc(N, p, Ma, Mb, Pab_Hom_L, Pab_Hom_L);

% Homogeneous Multi-Plex Network (Hom_P)
[L2_Hom_P, Sp_Hom_P] = calc(N, p, Ma, Mb, Pab_Hom_P, Pab_Hom_P);

%% Calculate standard deviation of all interlayer edge weights

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

stdTab = [std_Het_L std_Het_P; std_Hom_L std_Hom_P]; % tabulate

%% Calculate the percentage of missing interlayer edges (PercentME) for each network

PercentME_Het_L = (numel(Pab_Het_L)-nnz(Pab_Het_L))/numel(Pab_Het_L);
PercentME_Het_P = (numel(Pab_Het_P)-nnz(Pab_Het_P))/numel(Pab_Het_P);
PercentME_Hom_L = (numel(Pab_Hom_L)-nnz(Pab_Hom_L))/numel(Pab_Hom_L);
PercentME_Hom_P = (numel(Pab_Hom_P)-nnz(Pab_Hom_P))/numel(Pab_Hom_P);

percentMETab = [PercentME_Het_L PercentME_Het_P; PercentME_Hom_L PercentME_Hom_P]; % tabulate

%% plotting the graph

plot(p, L2_Het_L, p, L2_Het_P, p, L2_Hom_L, p, L2_Hom_P) % p as x axis
legend('Heterogeneous Multilayer', 'Heterogeneous Multiplex', 'Homogeneous Multilayer', 'Homogeneous Multiplex')
xlabel('p');
ylabel('Lambda2 values');

%% Function that calculates the Laplacian matrix for each network and loops through Sp

function [L2, Sp] = calc(N, p, Ma, Mb, Pab, Pba)
    
    len = length(p);
    L2 = ones(1,len);
    Sp = ones(1,len);

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
        
        % Sp
        sp = sum(p(k)*Pab, "all");
        Sp(k) = sp;

    end
end


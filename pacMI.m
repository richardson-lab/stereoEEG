function pacMI
% function pacMI
%   Phase-amplitude coupling (PAC) significance for sEEG data, following
%   methods of Mukamel et al 2014 J Neurosci. Computes and evaluates
%   significance of Tort's modulation index (MI) for specified channel x
%   phase frequency x amplitude frequency x time step.
%
%   DR 04/2023

% parameters
base = '/Users/drew/Library/CloudStorage/Box-Box/sEEG/';
subj = 'HUP241_RID890'; % session name
chnm = 'RB01'; % channel name
fP = logspace(log10(1),log10(20),21); % phase: frequencies (Hz)
fA = logspace(log10(5),log10(200),21); % amplitude: frequencies (Hz)
delta = pi/8; % phase bin size (rad)
cph = 1; % correct phase for waveform asymmetry? (0 or 1)
tbin = 120; % time bin length (s)
tstp = 30; % time step size (s)
Nperm = 200; % number of random permutations
tperm = [-60 60]; % interval of random permutation time shift (s)

% load data
cd(fullfile(base,subj));
load([subj '_Induction.mat']);
    
% remove bad channels and white matter channels
ibad = notSEEGchannels(channel_labels,1,bad_channels);
data(:,ibad) = [];
channel_labels(ibad,:) = [];
Nch = size(data,2);

% line noise filter
for ich = 1:Nch
    ln = mtmlinenoise(data(:,ich),3,sample_rate,sample_rate,60:60:max(fA));
    data(:,ich) = data(:,ich) - ln;
end

% re-reference to common average of each lead
ind = 1;
while ind < Nch
    l = channel_labels{ind,1};
    lead = l(isstrprop(l,'alpha'));
    ilead = find(strncmp(channel_labels,lead,length(lead)));
    if length(ilead) > 1
        disp(['re-referencing ' lead ' x ' num2str(length(ilead))]);
        data(:,ilead) = data(:,ilead) - mean(data(:,ilead),2);
    end
    ind = ilead(end) + 1; % assumes labels are sorted
end

% channel to analyze
ich = find(strcmp(channel_labels(:,1),chnm), 1);
if isempty(ich), error('no channel found'); end

% phase-amplitude coupling histograms
sbin = round(tbin*sample_rate); % samples per bin
sstp = round(tstp*sample_rate); % samples per step
rT = 1:sstp:size(data,1)-sbin; % start sample of each bin
Nt = length(rT);
t = (rT+sbin/2)/sample_rate; % s
rP = [fP(1:end-1)' fP(2:end)']; % frequency bin ranges
rA = [fA(1:end-1)' fA(2:end)'];
NP = size(rP,1);
NA = size(rA,1);
edges = -pi:delta:pi; % phase bins
x = edges(1:end-1)+delta/2; % radians
Nx = length(x);
MI = zeros(NP,NA,Nt);
MIpval = zeros(NP,NA,Nt);
for iP = 1:NP
    tic
    [bP,aP] = butter(2,rP(iP,:)/(sample_rate/2));
    P = angle(hilbert(filtfilt(bP,aP,data(:,ich)))); % phase x sample
    if cph
        ECDFt = sort(P);
        ECDFx = (1:length(ECDFt))/length(ECDFt);
        cx = interp1(ECDFt,ECDFx,P,'linear');
        P = 2*pi*cx-pi; % phase corrected for waveform asymmetry (Siapas et al 2005)
    end
    for iA = 1:NA
        [bA,aA] = butter(2,rA(iA,:)/(sample_rate/2));
        A = abs(hilbert(filtfilt(bA,aA,data(:,ich)))); % amplitude x sample
        for iT = 1:Nt
            t1 = rT(iT);
            t2 = rT(iT)+sbin-1;
            cP = P(t1:t2);
            cA = A(t1:t2);
            [~,bin] = histc(cP,edges);
            cAm = zeros(1,Nx);
            for jj = 1:Nx
                cAm(jj) = trimmean(cA(bin==jj),10); % robust mean amplitude in each phase bin
            end
            cAm = cAm/sum(cAm);
            cAm(cAm==0) = 1e-10; % to avoid -Inf for log(0)
            MI(iP,iA,iT) = (log(Nx)+sum(cAm.*log(cAm)))/log(Nx); % modulation index (Tort et al 2010)
        end
        MIperm = zeros(Nperm,Nt);
        for iperm = 1:Nperm
            tshift = tperm(1) + diff(tperm)*rand; % random time shift (s)
            nshift = round(tshift*sample_rate); % random time shift (samples)
            Ashift = circshift(A,nshift);
            for iT = 1:Nt
                t1 = rT(iT);
                t2 = rT(iT)+sbin-1;
                cP = P(t1:t2);
                cA = Ashift(t1:t2);
                [~,bin] = histc(cP,edges);
                cAm = zeros(1,Nx);
                for jj = 1:Nx
                    cAm(jj) = trimmean(cA(bin==jj),10); % robust mean amplitude in each phase bin
                end
                cAm = cAm/sum(cAm);
                cAm(cAm==0) = 1e-10; % to avoid -Inf for log(0)
                MIperm(iperm,iT) = (log(Nx)+sum(cAm.*log(cAm)))/log(Nx); % modulation index (Tort et al 2010)
            end
        end
        n = sum(MIperm >= ones(Nperm,1)*squeeze(MI(iP,iA,:))');
        MIpval(iP,iA,:) = (n+1)/(Nperm+1);
    end
    et = toc;
    disp([num2str(iP) '/' num2str(NP) ' : ' num2str(et) 's']); % progress report
end
eval(['PACmi.' chnm ' = MI;']);
eval(['PACmi.' chnm 'pval = MIpval;']);

% save
save([subj '_Induction.mat'],'PACmi','-append');

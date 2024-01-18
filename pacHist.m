function pacHist
% function pacHist
%   Phase-amplitude coupling (PAC) analysis for sEEG data. Computes PAC
%   histograms for each channel x phase frequency x amplitude frequency x
%   time step. Also computes Tort's modulation index and, optionally, its
%   significance based on permutation test of Mukamel et al 2014 J
%   Neurosci.
%
%   DR 04/2023

% parameters
fP = logspace(log10(1),log10(20),21); % phase: frequencies (Hz)
fA = logspace(log10(5),log10(200),21); % amplitude: frequencies (Hz)
delta = pi/8; % phase bin size (rad)
cph = 1; % correct phase for waveform asymmetry? (0 or 1)
tbin = 120; % time bin length (s)
tstp = 30; % time step size (s)
Nperm = 0; % number of permutations for MI significance test (leave empty or set to zero to not compute MIpval which takes a long time)
tperm = [-60 60]; % interval of random permutation time shift (s)

%
parpool(maxNumCompThreads);

% load data
selpath = uigetdir([],'Select top directory of sEEG data');
datfile = dir(fullfile(selpath,'**','*_Induction.mat'));
for ifile = 1:length(datfile)
    cd(datfile(ifile).folder);
%     load(datfile(ifile).name,'bad_channels','channel_labels','data','sample_rate');
    S = load(datfile(ifile).name,'bad_channels','channel_labels','data','sample_rate');
    [bad_channels,channel_labels,data,sample_rate] = v2struct(S);
    disp(datfile(ifile).name);
    
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

    % diagnostics
    % pwelchPlot(data,sample_rate,[0.5 500],5,1,500,channel_labels); % PSD of cleaned data (takes a long time)

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
    PAC = zeros(Nch,NP,NA,Nt,Nx);
    PACmi = zeros(Nch,NP,NA,Nt,2);
    for ich = 1:Nch
        tic
        for iP = 1:NP
            [bP,aP] = butter(2,rP(iP,:)/(sample_rate/2));
            P = angle(hilbert(filtfilt(bP,aP,data(:,ich)))); % phase x sample
            if cph
                ECDFt = sort(P);
                ECDFx = (1:length(ECDFt))/length(ECDFt);
                cx = interp1(ECDFt,ECDFx,P,'linear');
                P = 2*pi*cx-pi; % phase corrected for waveform asymmetry (Siapas et al 2005)
            end
            Pbin = zeros(sbin,Nt);
            for iT = 1:Nt
                t1 = rT(iT);
                t2 = rT(iT)+sbin-1;
                cP = P(t1:t2);
                [~,~,Pbin(:,iT)] = histcounts(cP,edges);
            end
            for iA = 1:NA
                [bA,aA] = butter(2,rA(iA,:)/(sample_rate/2));
                A = abs(hilbert(filtfilt(bA,aA,data(:,ich)))); % amplitude x sample
                for iT = 1:Nt
                    t1 = rT(iT);
                    t2 = rT(iT)+sbin-1;
                    cA = A(t1:t2);
                    cAm = zeros(1,Nx);
                    for jj = 1:Nx
                        cAm(jj) = mean(cA(Pbin(:,iT)==jj)); % mean amplitude in each phase bin
                    end
                    PAC(ich,iP,iA,iT,:) = cAm;
                    cAm = cAm/sum(cAm);
                    cAm(cAm==0) = 1e-10; % to avoid -Inf for log(0)
                    PACmi(ich,iP,iA,iT,1) = (log(Nx)+sum(cAm.*log(cAm)))/log(Nx); % modulation index (Tort et al 2010)
                end
                if Nperm>1 % permutation test for significance of MI
                    MIperm = zeros(Nperm,Nt);
                    parfor iperm = 1:Nperm
                        tshift = tperm(1) + diff(tperm)*rand; % random time shift (s)
                        nshift = round(tshift*sample_rate); % random time shift (samples)
                        Ashift = circshift(A,nshift);
                        for iT = 1:Nt
                            t1 = rT(iT);
                            t2 = rT(iT)+sbin-1;
                            cA = Ashift(t1:t2);
                            cAm = zeros(1,Nx);
                            for jj = 1:Nx
                                cAm(jj) = mean(cA(Pbin(:,iT)==jj));
                            end
                            cAm = cAm/sum(cAm);
                            cAm(cAm==0) = 1e-10;
                            MIperm(iperm,iT) = (log(Nx)+sum(cAm.*log(cAm)))/log(Nx);
                        end
                    end
                    n = sum(MIperm >= ones(Nperm,1)*squeeze(PACmi(ich,iP,iA,:,1))');
                    PACmi(ich,iP,iA,:,2) = (n+1)/(Nperm+1);
                end
            end
        end
        et = toc;
        disp([num2str(ifile) '/' num2str(length(datfile)) ' : ' num2str(ich) '/' num2str(Nch) ' : ' num2str(et) 's']); % progress report
    end

    % save
    PACparam.rP = rP;
    PACparam.rA = rA;
    PACparam.delta = delta;
    PACparam.cph = cph;
    PACparam.tbin = tbin;
    PACparam.tstp = tstp;
    PACparam.Nperm = Nperm;
    PACparam.tperm = tperm;
    PACparam.t = t;
    PACparam.x = x;
    PACparam.channel_labels = channel_labels;
    % save(datfile(ifile).name,'PAC','PACmi','PACparam','-append');
end

delete(gcp);
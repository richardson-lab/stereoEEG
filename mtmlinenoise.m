function ln = mtmlinenoise(x,NW,Nwin,Fs,Fsine)
% function ln = mtmlinenoise(x,NW,Nwin,Fs,Fsine)
%   Estimate line noise (sine waves of known frequency, e.g. 60 Hz) from a
%   data sequence using Thomson's multitaper method (see Mitra and Pesaran
%   1999) 
%
%   x       data sequence
%   NW      time-bandwidth product
%   Nwin    number of window samples
%   Fs      sampling frequency (Hz)
%   Fsine   sine wave frequencies (Hz)
%   ln      linear superposition of fitted sine waves
%
%   DR 01/2009

N = length(x);
x = x(:);
x = x - mean(x); % dc offsets can cause problems with fit
nfft = Fs;
if (Nwin > nfft) || (Nwin > N),
    error('parameters must satisfy: Nwin <= Fs and Nwin <= length(x)');
end
E = dpss(Nwin,NW);
U = fft(E,Nwin);
k = min(round(2*NW),Nwin); 
k = max(k-1,1);
U = U(:,1:k);
f = Fs/2*linspace(0,1,ceil(nfft/2)+1)';
if size(Fsine,2) == 1, Fsine = Fsine'; end
Nsine = length(Fsine);
indsine = zeros(Nsine,1);
for ii = 1:length(Fsine)
    [~,indsine(ii)] = min(abs(f-Fsine(ii)));
end
ln = zeros(size(x));
blk = 1:Nwin:N;
for ii = 1:length(blk)
    ind = blk(ii):min(blk(ii)+Nwin-1,N);
    if (length(ind) < Nwin) % what to do about last segment
        if (length(ind)/2 > NW) 
            E = dpss(length(ind),NW);
            U = fft(E,length(ind)); U = U(:,1:k);
        else continue;
        end
    end
%     x(ind) = x(ind) - median(x(ind)); % dc offsets can cause problems with fit
    Xx = fft(E(:,1:k).*x(ind,ones(1,k)),nfft); % windowed DFTs
    AP = sum(Xx(indsine,:).*(ones(Nsine,1)*conj(U(1,:))),2)/sum(abs(U(1,:)).^2); % amplitude, phase @ Fsine
    ln(ind) = exp(1i*2*pi*(0:length(ind)-1)'*Fsine/Fs)*AP + exp(-1i*2*pi*(0:length(ind)-1)'*Fsine/Fs)*conj(AP); % linear superposition of fitted sine waves
end
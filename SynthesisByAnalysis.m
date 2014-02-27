function [y yh ys]=SynthesisByAnalysis(pitch_file, wavfile, win, NFFT, nHarm, thsld, maxhd, stocf)
%=> analysis/synthesis of a sound using the harmonic plus stochastic model
%
% This is basically taken from the AMP course from Xavier Serra with small
% modifications
%
% pitch_file = pitch of the source which has to be synthesized
% wavfile = wave file which has to be analysed for extracting harmonic
%           amplitudes
% win = window that has to be used in analysis (odd size)
% NFFT = number of FFT points to be used in analysis (min 512)
% nHarm = number of harmonic to be resynthesized
% thshld = threshold below which peaks wont be considerd in spectrum (negative in dB)
% maxhd: max. relative deviation in harmonic detection (ex: .2)
% stocf: decimation factor of mag spectrum for stochastic analysis
% y: output sound, yh: harmonic component, ys: stochastic component
%
% NOTE: if you are performing synthesis for listening dont use higher
% number of harmonics because it casuses distortion. This method is not
% very sophisticated, its either ok for loudness computation of single
% source or casual synthesis of lead instrument (with lesser number of harmonics)


%%%%% Reading audio %%%%
[audio fs]=wavread(wavfile);
x=audio(:,1);

%%%%%% Reading the information in the pitch file %%%%%%%%%%
[time_pitch]=dlmread(pitch_file);

% ----------Common params ------------
H = 256;                                 % hop size for analysis and synthesis
N2 = NFFT/2+1;                              % half-size of spectrum
soundlength = length(x);                 % length of input sound

% ----------Analysis related Params ------------
M = length(win);                           % analysis window size
hM = (M-1)/2;                            % half analysis window size
fftbuffer = zeros(NFFT,1);                  % initialize buffer for FFT
win = win/sum(win);                            % normalize analysis window

% ----------Synthesis related Params ------------
Ns = 1024;                               % FFT size for synthesis
hNs = Ns/2;                              % half synthesis window size
yh = zeros(soundlength+Ns/2,1);       % output sine component
ys = zeros(soundlength+Ns/2,1);       % output residual component
sw = zeros(Ns,1);
ow = triang(2*H-1);                      % overlapping window
ovidx = Ns/2+1-H+1:Ns/2+H;               % overlap indexes
sw(ovidx) = ow(1:2*H-1);
bh = blackmanharris(Ns);                 % synthesis window
bh = bh ./ sum(bh);                      % normalize synthesis window
wr = bh;                                 % window for residual
sw(ovidx) = sw(ovidx) ./ bh(ovidx);
sws = H*hanning(Ns)/2;               % synthesis window for stochastic
lastyhloc = zeros(nHarm,1);             % initialize synthesis harmonic locations
yhphase = 2*pi*rand(nHarm,1);           % initialize synthesis harmonic phases


frame_cnt=1;
pin = max(hNs+1,1+hM);   % initialize sound pointer to middle of analysis window
pend = soundlength-max(hM,hNs);            % last sample to start a frame

% ----------Interpolating the pitch samples to suite the hop of the
% anlaysis----------
pitch_resampled = interp1(time_pitch(:,1), time_pitch(:,2), (pin:H:pend)/fs);
pitch_resampled(pitch_resampled<50)=0;


while pin<pend
    %-----analysis-----%
    xw = x(pin-hM:pin+hM).*win(1:M);         % window the input sound
    fftbuffer(:) = 0;                      % reset buffer
    fftbuffer(1:(M+1)/2) = xw((M+1)/2:M);  % zero-phase window in fftbuffer
    fftbuffer(NFFT-(M-1)/2+1:NFFT) = xw(1:(M-1)/2);
    X = fft(fftbuffer);                    % compute the FFT
    mX = 20*log10(abs(X(1:N2)));           % magnitude spectrum
    pX = unwrap(angle(X(1:NFFT/2+1)));        % unwrapped phase spectrum
    ploc = 1 + find((mX(2:N2-1)>thsld) .* (mX(2:N2-1)>mX(3:N2)) ...
        .* (mX(2:N2-1)>mX(1:N2-2)));     % find peaks
    [ploc,pmag,pphase] = peakinterp(mX,pX,ploc);     % refine peak values
    yinws = round(fs*0.0125);          % using approx. a 12.5 ms window for yin
    yinws = yinws+mod(yinws,2);        % make it even
    yb = pin-yinws/2;
    ye = pin+yinws/2+yinws;
    %   if (yb<1 || ye>length(x))          % out of boundaries
    %      f0 = 0;
    %   else
    %      f0 = f0detectionyin(x(yb:ye),fs,yinws,minf0,maxf0); % compute f0
    %   end
    f0 = pitch_resampled(frame_cnt);
    hloc = zeros(nHarm,1);                    % initialize harmonic locations
    hmag = zeros(nHarm,1)-100;                % initialize harmonic magnitudes
    hphase = zeros(nHarm,1);                  % initialize harmonic phases
    hf = (f0>0).*(f0.*(1:nHarm));             % initialize harmonic frequencies
    hi = 1;                                % initialize harmonic index
    npeaks = length(ploc);                 % number of peaks found
    while (f0>0 && hi<=nHarm && hf(hi)<fs/2)  % find harmonic peaks
        [dev,pei] = min(abs((ploc(1:npeaks)-1)/NFFT*fs-hf(hi)));   % closest peak
        if ((hi==1 || ~any(hloc(1:hi-1)==ploc(pei))) && dev<maxhd*hf(hi))
            hloc(hi) = ploc(pei);              % harmonic locations
            hmag(hi) = pmag(pei);              % harmonic magnitudes
            hphase(hi) = pphase(pei);          % harmonic phases
        end
        hi = hi+1;                           % increase harmonic index
    end
    hloc(1:hi-1) = (hloc(1:hi-1)~=0).*((hloc(1:hi-1)-1)*Ns/NFFT);   % synth. locs
    ri= pin-hNs;                     % input sound pointer for residual analysis
    
    if(stocf>0)        
        xr = x(ri:ri+Ns-1).*wr(1:Ns);           % window the input sound
        Xr = fft(fftshift(xr));                 % compute FFT for residual analysis
        Xh = genspecsines(hloc(1:hi-1),hmag,hphase,Ns);            % generate sines
        Xr = Xr-Xh;                             % get the residual complex spectrum
        mXr = 20*log10(abs(Xr(1:Ns/2+1)));      % magnitude spectrum of residual
        mXsenv = decimate(max(-200,mXr),stocf); % decimate the magnitude spectrum
    end                             % and avoid -Inf
    %-----synthesis data-----%
    yhloc = hloc;                           % synthesis harmonics locs
    yhmag = hmag;                           % synthesis harmonic amplitudes
    yf0 = f0;                               % synthesis f0
    if(stocf>0)
        mYsenv = mXsenv;                        % synthesis residual envelope
    end
    %-----transformations-----%
    %-----synthesis-----%
    yhphase = yhphase + 2*pi*(lastyhloc+yhloc)/2/Ns*H;   % propagate phases
    lastyhloc = yhloc;
    Yh = genspecsines(yhloc,yhmag,yhphase,Ns);    % generate sines
    yhw = fftshift(real(ifft(Yh)));       % sines in time domain using IFFT
    yh(ri:ri+Ns-1) = yh(ri:ri+Ns-1)+yhw(1:Ns).*sw;  % overlap-add for sines
    
    if(stocf>0)
        mYs = interp(mYsenv,stocf);                   % interpolate to original size
        roffset = ceil(stocf/2)-1;                    % interpolated array offset
        mYs = [ mYs(1)*ones(roffset,1); mYs(1:Ns/2+1-roffset) ];
        mYs = 10.^(mYs/20);                           % dB to linear magnitude
        if (f0>0)
            mYs = mYs .* cos(pi*[0:Ns/2]'/Ns*fs/yf0).^2;  % filter residual
        end
        fc = 1+round(500/fs*Ns);                          % 500 Hz
        mYs(1:fc) = mYs(1:fc) .* ([0:fc-1]'/(fc-1)).^2;   % HPF
        pYs = 2*pi*rand(Ns/2+1,1);      % generate phase spectrum with random values
        mYs1 = [mYs(1:Ns/2+1); mYs(Ns/2:-1:2)]; % create complete magnitude spectrum
        pYs1 = [pYs(1:Ns/2+1); -1*pYs(Ns/2:-1:2)];  % create complete phase spectrum
        Ys = mYs1.*cos(pYs1)+1i*mYs1.*sin(pYs1);    % compute complex spectrum
        ysw = fftshift(real(ifft(Ys)));       % stochastic in time domain using IFFT
        ys(ri:ri+Ns-1) = ys(ri:ri+Ns-1)+ysw(1:Ns).*sws; % overlap-add for stochastic
    end
    pin = pin+H;                                  % advance the sound pointer
    frame_cnt = frame_cnt +1;
end
y= yh+ys;                                         % sum sines and stochastic
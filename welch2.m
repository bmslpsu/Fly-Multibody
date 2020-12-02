function varargout = welch2(x,funcName,varargin)
%#codegen
%WELCH Welch spectral estimation method.
%   [Pxx,F] = WELCH(X,'pwelch',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,TRACE,'psd')
%   [Pxx,F] = WELCH(X,'pwelch',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,TRACE,'ms')
%   [Pxx,F] = WELCH(X,'pwelch',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,TRACE,'power')
%   [Pxy,F] = WELCH({X,Y},'cpsd',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE)
%   [Txy,F] = WELCH({X,Y},'tfe',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE)
%   [Txy,F] = WELCH({X,Y},'tfeh2',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE)
%   [Cxy,F] = WELCH({X,Y},'mscohere',WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE)
%   [...,F] = WELCH({X,Y},'mscohere',WINDOW,NOVERLAP,NFFT,Fs,...'mimo')
%   [Pxx,F,Pxxc] = WELCH(...)
%   [Pxx,F,Pxxc] = WELCH(...,'ConfidenceLevel',P)
%
%   Inputs:
%      see "help pwelch" for complete description of all input arguments.
%      ESTTYPE - is a string specifying the type of estimate to return, the
%                choices are: psd, cpsd, tfe, tfeh2, and mscohere.
%
%   Outputs:
%      Depends on the input string ESTTYPE:
%      Pxx - Power Spectral Density (PSD) estimate, or
%      MS  - Mean-square spectrum, or
%      Pxy - Cross Power Spectral Density (CPSD) estimate, or
%      Txy - Transfer Function Estimate (TFE), or
%      Cxy - Magnitude Squared Coherence.
%      F   - frequency vector, in Hz if Fs is specified, otherwise it has
%            units of rad/sample

%   Copyright 1988-2020 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and
%         Modeling, John Wiley & Sons, 1996.

narginchk(2,11);
nargoutchk(0,3);

% for codegen build: make sure that the C/C++ code is generated only for
% one computation function
assert(coder.internal.isConst(funcName));

if strcmp(funcName,'pwelch')
    % look for psd, power, and ms window compensation flags for pwelch
    [esttype, args] = signal.internal.psdesttype({'psd','power','ms'},'psd',varargin);
else
    esttype = funcName;
    args = varargin;
end

% Parse input arguments.
[x1,~,~,y1,~,win1,winName,winParam,noverlap,k1,L,options] = ...
    signal.internal.spectral.welchparse(x,esttype,args{:});
% Cast to enforce precision rules
options.nfft = signal.internal.sigcasttofloat(options.nfft,'double',...
  'WELCH','NFFT','allownumeric');
noverlap = signal.internal.sigcasttofloat(noverlap,'double','WELCH',...
  'NOVERLAP','allownumeric');
options.Fs = signal.internal.sigcasttofloat(options.Fs,'double','WELCH',...
  'Fs','allownumeric');
k = double(k1);

if any([signal.internal.sigcheckfloattype(x1,'single')...
    signal.internal.sigcheckfloattype(y1,'single'),...
    isa(win1,'single')])
  x = single(x1);
  y = single(y1);
  win = single(win1);
else
    x = x1;
    y = y1;
    win = win1;
end

% Frequency vector was specified, return and plot two-sided PSD
freqVectorSpecified = false; 
nrow = 1;
if length(options.nfft) > 1
    freqVectorSpecified = true;
    [~,nrow] = size(options.nfft);
end

% Compute the periodogram power spectrum of each segment and average always
% compute the whole power spectrum, we force Fs = 1 to get a PS not a PSD.

% Initialize
if freqVectorSpecified
    if strcmpi(options.range,'onesided')
        coder.internal.warning('signal:welch:InconsistentRangeOption');
    end
    options.range = 'twosided';
end
% Cast to enforce precision rules
zeroPrototype = zeros(1,class(x)); %#ok<*ZEROLIKE>

LminusOverlap = L-noverlap;
xStart = 1:LminusOverlap:k*LminusOverlap;
xEnd   = xStart+L-1;

method = "plus";
cmethod = @plus;

% maximum cache to compute MIMO cpsd and mscohere
cacheSize = 2e9; % 2 GB

isyEmpty = isempty(y);

switch funcName
    case {'pwelch'}
        % Combining method
        if options.maxhold
            method = "max";
            cmethod1 = @(seg,nextseg) max(seg,real(k*nextseg)); % k will be removed below
            % real is required for coder to infer it as double and not complex double
            [Pxx,w,units] = localComputeSpectra(-Inf(1,class(x)),x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod1,freqVectorSpecified);
        elseif options.minhold
            method = "min";
            cmethod2 = @(seg,nextseg) min(seg,real(k*nextseg)); % k will be removed below
            % real is required for coder to infer it as double and not complex double
            [Pxx,w,units] = localComputeSpectra(Inf(1,class(x)),x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod2,freqVectorSpecified);
        else
            [Pxx,w,units] = localComputeSpectra(zeroPrototype,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
        end
        
    case 'cpsd'
        [numChan, chindx, chindy] = localChannelIdx(x,y,options.MIMO);
        % get memory usage in bytes per channel for x, y
        xMemPerCh = getDataMemory(x);
        yMemPerCh = getDataMemory(y);
        
        % number of channel for computation under the cache size
        maxCol = getNumCh(cacheSize,xMemPerCh,yMemPerCh,numChan);      
        
        if maxCol == numChan || ~options.MIMO
            % if all channel can be calculated
            [Pxx,w,units] = localComputeSpectra(zeroPrototype,x(:,chindx),y(:,chindy),xStart,xEnd,win,...
                options,esttype,k,cmethod,freqVectorSpecified);
        else
            [Pxx,w,units] = computeMIMOForLoop(x,y,chindx,chindy,maxCol,numChan,xStart,...
                    xEnd,win,options,esttype,k,cmethod,freqVectorSpecified);           
        end
        
        % y will never be empty for 'tfe','tfe2','mscohere'. isempty(y) is
        % a guard for coder size inference.
      
    case 'tfe'
        if options.MIMO
            [~, chindy, chindx] = localChannelIdx(y,x,options.MIMO);
            [~, chindx1, chindx2] = localChannelIdx(x,x,options.MIMO);

            [Pxx,w,units] = localComputeSpectra(zeroPrototype,x(:,chindx1),x(:,chindx2),xStart,xEnd,win,...
                options,esttype,k,cmethod,freqVectorSpecified);
            % Cross PSD.  The frequency vector and xunits are not used.
            if ~isyEmpty
                Pxy = localComputeSpectra(zeroPrototype,y(:,chindy),x(:,chindx),xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
            else
                Pxy = zeros(0,class(x));
            end
        else
            [Pxx,w,units] = localComputeSpectra(zeroPrototype,x,[],xStart,xEnd,win,...
                options,esttype,k,cmethod,freqVectorSpecified);
            % Cross PSD.  The frequency vector and xunits are not used.
            if ~isyEmpty
                Pxy = localComputeSpectra(zeroPrototype,y,x,xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
            else
                Pxy = zeros(0,class(x));
            end
        end
          
        if options.MIMO && size(x,2) > 1
            Pxx = localComputeMIMO(x,y,Pxx,Pxy,[],'tfe');
        else
            Pxx = bsxfun(@rdivide, Pxy, Pxx);
        end
        
    case 'tfeh2'
        if options.MIMO
            [~, chindx, chindy] = localChannelIdx(x,y,options.MIMO);
            [~, chindy1, chindy2] = localChannelIdx(y,y,options.MIMO);

            if ~isyEmpty
                [Pyy,w,units] = localComputeSpectra(zeroPrototype,y(:,chindy1),y(:,chindy2),xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
                % Cross PSD.  The frequency vector and xunits are not used.
                Pxy = localComputeSpectra(zeroPrototype,x(:,chindx),y(:,chindy),xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
            else
                Pyy = zeros(0,class(x));
                Pxy = zeros(0,class(x));
                w = [];
                units = '';
            end
        else
            if ~isyEmpty
                [Pyy,w,units] = localComputeSpectra(zeroPrototype,y,[],xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
                % Cross PSD.  The frequency vector and xunits are not used.
                Pxy = localComputeSpectra(zeroPrototype,x,y,xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
            else
                Pyy = zeros(0,class(x));
                Pxy = zeros(0,class(x));
                w = [];
                units = '';
            end
        end
        
        if options.MIMO && size(x,2) > 1
            Pxx = localComputeMIMO(x,y,[],Pxy,Pyy,'tfeh2');
        else
            Pxx = bsxfun(@rdivide, Pyy, Pxy);
        end
        
    case 'mscohere'
        % Note: (Sxy1+Sxy2)/(Sxx1+Sxx2) != (Sxy1/Sxy2) + (Sxx1/Sxx2)
        % ie, we can't push the computation of Cxy into computeperiodogram.
        if options.MIMO
            [numChanxy, chindx, chindy] = localChannelIdx(x,y,options.MIMO);
            [numChanxx, chindx1, chindx2] = localChannelIdx(x,x,options.MIMO);
            
            % get memory usage in bytes per channel for x
            xMemPerCh = getDataMemory(x);

            % number of channel for computation under the cache size
            maxCol = getNumCh(cacheSize,xMemPerCh,xMemPerCh,numChanxx);  
            
            if maxCol == numChanxx
                [Pxx,w,units] = localComputeSpectra(zeroPrototype,x(:,chindx1),x(:,chindx2),xStart,xEnd,win,...
                options,esttype,k,cmethod,freqVectorSpecified);
            else
                [Pxx,w,units] = computeMIMOForLoop(x,x,chindx1,chindx2,maxCol,numChanxx,xStart,...
                    xEnd,win,options,esttype,k,cmethod,freqVectorSpecified);
            end           
            
            if ~isyEmpty
                Pyy = localComputeSpectra(zeroPrototype,y,[],xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
                % Cross PSD.  The frequency vector and xunits are not used.
                
                % get memory usage in bytes per channel for x
                yMemPerCh = getDataMemory(y);

                % number of channel for computation under the cache size
                maxCol = getNumCh(cacheSize,xMemPerCh,yMemPerCh,numChanxy);  
                
                if maxCol == numChanxy
                    Pxy = localComputeSpectra(zeroPrototype,x(:,chindx),y(:,chindy),xStart,xEnd,win,...
                        options,esttype,k,cmethod,freqVectorSpecified);
                else
                    Pxy = computeMIMOForLoop(x,y,chindx,chindy,maxCol,numChanxy,xStart,...
                        xEnd,win,options,esttype,k,cmethod,freqVectorSpecified);
                end 
            else
                Pyy = zeros(0,class(x));
                Pxy = zeros(0,class(x));
            end
        else
            [Pxx,w,units] = localComputeSpectra(zeroPrototype,x,[],xStart,xEnd,win,...
                options,esttype,k,cmethod,freqVectorSpecified);
            if ~isyEmpty
                Pyy = localComputeSpectra(zeroPrototype,y,[],xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
                % Cross PSD.  The frequency vector and xunits are not used.
                Pxy = localComputeSpectra(zeroPrototype,x,y,xStart,xEnd,win,...
                    options,esttype,k,cmethod,freqVectorSpecified);
            else
                Pyy = zeros(0,class(x));
                Pxy = zeros(0,class(x));
            end
        end
        
        if options.MIMO && size(x,2) > 1
            Pxx = localComputeMIMO(x,y,Pxx,Pxy,Pyy,'mscohere'); 
        else
            Pxx = (abs(Pxy).^2)./bsxfun(@times,Pxx,Pyy); % Cxy
        end 
        
    otherwise
           Pxx = zeros(0,class(x));
           w = [];
           method = "";
        
end

% Compute confidence intervals if needed.
if ~isnan(options.conflevel)
    if any(strcmpi(esttype,{'ms','power','psd'})) && method ~= "plus"
        % Always compute the confidence interval around the average spectrum
        [avgPxx,w] = localComputeSpectra(zeroPrototype,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
    else
        avgPxx = Pxx;
    end
    % Cast to enforce precision rules
    avgPxx = double(avgPxx);
    Pxxc = signal.internal.spectral.confInterval(options.conflevel, avgPxx, isreal(x), w, options.Fs, k);
elseif nargout>2
    if any(strcmpi(esttype,{'ms','power','psd'})) && method ~= "plus"
        % Always compute the confidence interval around the average spectrum
        [avgPxx,w] = localComputeSpectra(zeroPrototype,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
    else
        avgPxx = Pxx;
    end
    % Cast to enforce precision rules
    avgPxx = double(avgPxx);
    Pxxc = signal.internal.spectral.confInterval(0.95, avgPxx, isreal(x), w, options.Fs, k);
else
    Pxxc = [];
end

if nargout==0 && coder.target('MATLAB')
    signal.internal.spectral.plotWelch(Pxx,w,Pxxc,esttype,noverlap,L,winName,winParam,units,options);
else
    if options.centerdc
        [Pxx1, w1, Pxxc1] = signal.internal.spectral.psdcenterdc(Pxx, w, Pxxc, options, esttype);
    else
        Pxx1 = Pxx;
        w1 = w;
        Pxxc1 = Pxxc;
    end
    
    % If the input is a vector and a row frequency vector was specified,
    % return output as a row vector for backwards compatibility.
    if nrow > 1 && isvector(x)
        Pxx2 = Pxx1.'; 
        w2 = w1.';
    else
        Pxx2 = Pxx1;
        w2 = w1;
    end
    
    % Cast to enforce precision rules   
    % Only cast if output is requested, otherwise, plot using double
    % precision frequency vector.
    if isa(Pxx2,'single')
        w3 = single(w2);
    else
        w3 = w2;
    end
    
    % Reshape output if MIMO
    if options.MIMO
        switch lower(esttype)
            case {'tfe','tfeh2'}
                Pxx3 = reshape(Pxx2,size(Pxx2,1),size(y,2),size(x,2));
            case 'cpsd'
                Pxx3 = reshape(Pxx2,size(Pxx2,1),size(x,2),size(y,2));
            otherwise
                Pxx3 = Pxx2;
        end
    else
        Pxx3 = Pxx2;
    end
    
    if nargout < 2 
        varargout = {Pxx3,w3}; % Pxx=PSD, MEANSQUARE, CPSD, or TFE
    else
        varargout = {Pxx3,w3,Pxxc1};
    end       
    
end

end

function [Pxx,w,units] = localComputeSpectra(Sxx,x,y,xStart,xEnd,win,options,esttype,k,cmethod,freqVectorSpecified)
 
 w=[]; 
 if ~isempty(y)
      Sxx1 = zeros(0,'like',1i*Sxx); %initialization
     for ii = 1:k
         [Sxxk,w] = computeperiodogram({x(xStart(ii):xEnd(ii),:),...
             y(xStart(ii):xEnd(ii),:)},win,options.nfft,esttype,options.Fs);
         if ii == 1
             Sxx1 = cmethod(Sxx,Sxxk);
         else
             Sxx1 = cmethod(Sxx1,Sxxk);
         end
     end
     
 else
	Sxx1 = zeros(0,class(Sxx)); %initialization
    for ii = 1:k
         [Sxxk,w] = computeperiodogram2(x(xStart(ii):xEnd(ii),:),win,...
             options.nfft,esttype,options.Fs);
         if ii == 1
             % use Sxx for applying cmethod in first iteration
             Sxx1 = cmethod(Sxx,Sxxk);
         else
             % accumulate from the second iteration
             Sxx1  = cmethod(Sxx1,Sxxk);
         end
    end
 end

Sxx1 = Sxx1./k; % Average the sum of the periodograms

% Generate the freq vector directly in Hz to avoid roundoff errors due to
% conversions later.
if ~freqVectorSpecified
    w = psdfreqvec('npts',options.nfft, 'Fs',options.Fs);
end

% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, corresponding freq vector and freq units.
[Pxx,w,units] = computepsd(Sxx1,w,options.range,options.nfft,options.Fs,esttype);
end

function [numChan, chindx, chindy] = localChannelIdx(x,y,MIMO)
% Find the number of cross-spectra of the columns of x and y to compute,
% numChan, and the indices into the inputs x and y to compute each
% cross-spectrum. If MIMO is specified, every combination of the columns of
% x and y is needed. Otherwise, compute cross-spectra of column-wise pairs
% of column vectors in x and y.
Nx = size(x,2);
Ny = size(y,2);
if MIMO
    numChan = Nx*Ny;
    chindx = repmat(1:Nx,1,Ny);
    chindy1 = repmat(1:Ny,Nx,1);
    chindy = chindy1(:)';
else
    numChan = max(size(x,2),size(y,2));
    chindx = 1:Nx;
    chindy = 1:Ny;
end
end

function Pout = localComputeMIMO(x,y,Pxx,Pxy,Pyy,esttype)
% Compute the transfer function for a MIMO system. At each frequency, the
% transfer function H is given by H(f) = MPxy(f)/MPxx(f), where MPxy and
% MPxx are matrices reshaped from the columns of Pxy and Pxx. Note that Pxy
% is a misnomer here, since it is formed from the DFT of y times the
% complex conjugate of the DFT of x.
nIn = size(x,2);
nOut = size(y,2);

if coder.target('MATLAB')
    % Disable warnings for the case that Pxx is singular. This will avoid
    % repeated warnings.
    [msg0, id0] = lastwarn('');
    state(1) = warning('off','MATLAB:nearlySingularMatrix'); %#ok<*EMVDF>
    state(2) = warning('off','MATLAB:singularMatrix');
    cleanupObj = onCleanup(@()warning(state));
end

switch esttype
    case 'tfe'
        Pout = zeros(size(Pxy),'like',Pxy); %Txy
        for i = 1:size(Pxy,1)
            P = reshape(Pxy(i,:),nOut,nIn)/reshape(Pxx(i,:),nIn,nIn);
            Pout(i,:) = P(:);
        end
    case 'tfeh2'
        % Number of input and output channels are equal
        Pout = zeros(size(Pxy),'like',Pxy); %Txy
        for i = 1:size(Pxy,1)
            P = reshape(Pyy(i,:),nOut,nOut)/reshape(Pxy(i,:),nOut,nOut);
            Pout(i,:) = P(:);
        end
    case 'mscohere'  
        Pout = zeros(size(Pxx,1),nOut,'like',Pxy); %Cxy
        Pxx0 = reshape(Pxx,[],nIn,nIn);
        Pxx0 = permute(Pxx0,[2,3,1]);
        for iOut = 1:nOut
            Pyx0 = Pxy(:,(iOut-1)*nIn+1:(iOut-1)*nIn+nIn);
            Pxy0 = Pyx0.';
            Pyx0 = conj(Pyx0);
            Pyy0 = Pyy(:,iOut);
            for fNum = 1:size(Pxx,1)
                Pout(fNum,iOut) = real(Pyx0(fNum,:)/Pxx0(:,:,fNum)*Pxy0(:,fNum)/Pyy0(fNum));
            end
        end
end

if coder.target('MATLAB')
    % Warn if a singular matrix warning occured. Reset lastwarn if no warnings
    % occured.
    [~,msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:nearlySingularMatrix') || strcmp(msgid,'MATLAB:singularMatrix')
        coder.internal.warning('signal:welch:SingularMatrixMIMO');
    elseif isempty(msgid)
        lastwarn(msg0,id0);
    end
end

end

function mem = getDataMemory(x)
% Get the memory usage of each channel in x
if isa(x,'double')
    mem = 8;
end

if isa(x,'single')
    mem = 4;
end

if ~isreal(x)
    mem = mem*2;
end

mem = mem*size(x,1);
end

function numCol = getNumCh(cacheSize,xMem,yMem,totNumCh)
% Get the number of channels to be computed based on cacheSize and data
% memory usage for one channel

% Always set numCol as 1 if cacheSize<(xMem+yMem)
numCol = ceil(cacheSize/(xMem+yMem));
if totNumCh <= numCol
    numCol(1) = totNumCh;
end

end

function [P,w,units] = computeMIMOForLoop(x,y,chindx,chindy,maxCol,numChan,xStart,xEnd,win,options,esttype,k,cmethod,freqVectorSpecified)
% For-loop to calculate localComputeSpectra. In each iteration, numCol
% channels are computed.

% calculate the first channel (for codegen)
Stmp = zeros(1,class(x));
[Ptmp,w,units] = localComputeSpectra(Stmp,x(:,chindx(1)),...
    y(:,chindy(1)),xStart,xEnd,win,options,esttype,k,cmethod,freqVectorSpecified);

P = zeros(size(Ptmp,1),numChan,'like',Ptmp);
P(:,1) = Ptmp(:,1);

numGroups = ceil((numChan-1)/maxCol);

% for-loop
for idx = 1:numGroups
    if idx == numGroups
        iCol = (idx-1)*maxCol+2:min(idx*maxCol+1,numChan);
    else             
        iCol = (idx-1)*maxCol+2:idx*maxCol+1;
    end
    P(:,iCol) = localComputeSpectra(Stmp,x(:,chindx(iCol)),y(:,chindy(iCol)),xStart,xEnd,win,...
        options,esttype,k,cmethod,freqVectorSpecified);
end
end
% [EOF]

% LocalWords:  Pxx NOVERLAP NFFT Fs SPECTRUMTYPE ESTTYPE Pxy Txy tfeh Cxy mimo
% LocalWords:  Pxxc Petre Stoica Monson allownumeric xunits Sxy Sxx
% LocalWords:  computeperiodogram conflevel MEANSQUARE cmethod periodograms
% LocalWords:  roundoff npts Manolakis Ingle Kagon Graw MPxy MPxx DFT occured

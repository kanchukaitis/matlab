function [sd99,sd95,sd90,sd5,psig]=gsbtest(seriesLength,iters,window,seriesCorrelation,corrType,genMethod,obs)

% GSBTEST  'Gershunov' Test of Moving Correlation Significance Levels
%
%       [sd99,sd95,sd90,sd5,psig]=gsbtest(seriesLength,iters,window,seriesCorrelation,corrType,genMethod,obs)
%
% Implementation of a Monte Carlo type simulation of moving correlation functions
% between two weakly correlated noise series for the purposes of constructing
% significance tests for moving correlation functions (or traces) based on their standard
% deviation (Gershunov et al. 2001). Updated version now supports both Pearson and Spearman moving
% correlation via CORRTYPE.
%
% Inputs:
%   seriesLength         = length of the time series used in the running correlation analysis
%   iters                = number of Monte Carlo iterations to perform (should be >= 1,000)
%   window               = width (in units of the time step for the data) of the correlation window
%   seriesCorrelation    = desired overall correlation of the simulated time series (simulations will be within 0.5% of this value)
%   corrType [optional]  = 'pearson' (default) or 'spearman' -- type of moving correlation to test
%   genMethod [optional] = which method to use to create the correlated random series [default = 1], see below
%   obs [optional]       = find the significance level for a specific (observed) value of the std of the moving correlation function
%
% Series-Generation Methods (genMethod)
%   1 = Cholesky decomposition of the covariance matrix [default]
%   2 = Matrix square root
%
% Outputs:
%  sd99  = bootstrapped 99% confidence level
%  sd95  = bootstrapped 95% confidence level
%  sd90  = bootstrapped 90% confidence level
%  sd5   = bootstrapped 5% confidence level
%  psig  = significance for specific (observed) standard deviation (if input argument 'obs' is given)
%
% Example:
%  > [sd99,sd95]=gsbtest(100,1000,11,0.30,'pearson',2)
%
% sd99 =
%
%    0.4002
%
% sd95 =
%
%    0.3624
%
% NOTE: Although the solution should converge for iters >= 1000, small differences (due to the
% correlation between the two artificial time series generated in the script) may still result,
% especially for a small number (<<10000) of simulations; so it is recommended that you verify
% convergence with several runs and be cautious of values too close to your chosen test statistic
% (particularly at 90% or 95% CI). You can compare the above result for 'sd95' to Table 1 from
% Gershunov et al. 2001 (for interseries r = 0.30, window length = 11, their 95% confidence level = 0.36).
%
% References
% ----------
% Gershunov, A., N. Schneider, and T. Barnett, 2001. Low-frequency modulation of the ENSO-Indian
% monsoon rainfall relationship: Signal or noise?, Journal of Climate, 14: 2486-2492.
%

% Update History:
%
% 07.26b -- modernized several earlier era MATLAB patterns:
%             (1) randn('state', sum(100*clock)) (deprecated since R2011a)
%               replaced with rng('shuffle').
%             (2) manual repmat()-based mean broadcasting in pairwise() and
%               pairwiseSpearman() replaced with implicit expansion
%               (available since R2016b; repmat was the standard
%               workaround before then).
%             (3) disp() console message on rejection-sampling failure
%               replaced with a proper identified warning()
%               ('gsbtest:noConvergence'), which is filterable/suppressible
%               via warning('off', ...) and visible in more contexts than
%               a bare disp().
%             (4) legacy scalar '|' replaced with 'nargout <= 1' in
%               findnearest()'s output-count check.
%           Note: this version targets MATLAB R2016b+ for implicit
%           expansion, and will not work in Octave or prior to MATLAB R2026b
% 07.26a -- merged gsbtest.m and gsbtest_spearman.m into a single file with a
%          corrType switch ('pearson'/'spearman'); fixed three bugs found in
%          the process:
%            (1) gsbtest_spearman.m's Cholesky branch had a stray '%' that
%                silently commented out '*seriesCorrelation' in the
%                rejection-sampling tolerance check, so its calibration
%                tolerance was 0.01 (absolute) instead of the intended
%                0.005*seriesCorrelation (relative) -- much looser than the
%                Pearson version, especially at low target correlations.
%            (2) gsbtest_spearman.m allowed only 10 rejection-sampling
%                retries per iteration (forgive=10) vs. gsbtest.m's 1000,
%                making it far more likely to silently accept a poorly
%                calibrated surrogate pair. Both now use 1000 by default.
%            (3) gsbtest_spearman.m's matrix-square-root branch (genMethod=2)
%                checked convergence with plain Pearson corrcoef() even
%                though the eventual moving-correlation statistic was
%                Spearman -- so it was rejection-sampling toward the wrong
%                target. Convergence checks now always match the corrType.
%          Also switched pairwise() to nanmean/nanstd (matching
%          pairwiseSpearman(), and the movingcorrelation.m this is meant to
%          be consistent with) for robustness to missing data.
% 09 13 -- minor code improvements for clarity; added 'findnearest' as subfunction
% 03.11 -- added calculation of 'psig' using input 'obs'
% 06.08 -- posted online
% 10.07 -- initial public version
% 09.03 -- initial version
%
% Function gsbtest.m written by Kevin Anchukaitis
% Subfunction build_win is a modification of the function by David Meko 
% Subfunction pairwise is a modification of corrpair.m by David Meko 
% Included subfunction 'findnearest' by Tom Benson

if nargin<7; obs = [];  end
if nargin<6 || isempty(genMethod); genMethod = 1; end
if nargin<5 || isempty(corrType);  corrType = 'pearson'; end
if nargin<4; error('Too Few Input Parameters'); end

corrType = lower(corrType);
if ~any(strcmp(corrType, {'pearson','spearman'}))
    error('corrType must be ''pearson'' or ''spearman''.');
end
if ~any(genMethod == [1 2])
    error('Invalid Series Creation Method (genMethod must be 1 or 2)');
end

if nargout==5 && isempty(obs)
    error('Need to provide desired (observed) value in ''obs'' to output significance level')
end

% Create a correlation matrix from the desired correlation value
Rmat = [1 seriesCorrelation; seriesCorrelation 1];

rng('shuffle');  % seed the global random stream from the current time

% Rejection-sampling tolerance for the realized correlation of each simulated pair
tol     = 0.005;  % relative tolerance (fraction of seriesCorrelation)
forgive = 1000;   % max rejection-sampling retries per iteration before giving up

mcf_var = zeros(iters,1);
mcf_std = zeros(iters,1);

for i = 1:iters
    rcheck = 1;
    count_repeated_checks = 0;
    while abs(rcheck) > tol*seriesCorrelation
        X = generate_pair(Rmat, seriesLength, genMethod);
        rcheck = sample_corr(X, corrType) - seriesCorrelation;
        count_repeated_checks = count_repeated_checks + 1;
        if count_repeated_checks >= forgive
            warning('gsbtest:noConvergence', ...
                'stuck on recheck = %.4f, giving up ...', rcheck)
            break
        end
    end

    % moving correlation function on the (accepted) series in X
    x1 = build_win(X(:,1),window);
    x2 = build_win(X(:,2),window);
    s  = moving_corr(x1, x2, corrType);

    mcf_var(i) = var(s);
    mcf_std(i) = std(s);
end

% sort the vector of standard deviations
sorted_mcf_std = sort(mcf_std);

% get the desired values
sd99 = sorted_mcf_std(round(0.99 * iters));
sd95 = sorted_mcf_std(round(0.95 * iters));
sd90 = sorted_mcf_std(round(0.90 * iters));
sd5  = sorted_mcf_std(round(0.05 * iters));

% calculate the significance level for the observed, if passed as input
if ~isempty(obs)
    flag_pval = findnearest(obs,sorted_mcf_std);
    psig = (flag_pval)/length(sorted_mcf_std);
end


%% subfunctions

function X = generate_pair(Rmat, n, genMethod)
% GENERATE_PAIR  Draw one n-sample pair of Gaussian series with covariance Rmat,
%                using either a Cholesky decomposition (genMethod=1) or a
%                matrix square root (genMethod=2).
switch genMethod
    case 1
        C = chol(Rmat);          % Cholesky decomposition of the correlation matrix
        Z = randn(2,n);
        X = (C'*Z)';
    case 2
        X = (sqrtm(Rmat)*randn(2,n))';
end


function r = sample_corr(X, corrType)
% SAMPLE_CORR  Realized overall correlation of a 2-column series pair X,
%              computed with whichever correlation type (Pearson or
%              Spearman) matches the eventual test statistic -- this is
%              what the rejection-sampling loop checks against the target.
switch corrType
    case 'pearson'
        c = corrcoef(X);
    case 'spearman'
        c = corr(X,'Type','Spearman','rows','pairwise');
end
r = c(2,1);


function s = moving_corr(X, Y, corrType)
% MOVING_CORR  Dispatch to the Pearson or Spearman windowed-correlation
%              subfunction, matching corrType.
switch corrType
    case 'pearson'
        s = pairwise(X,Y);
    case 'spearman'
        s = pairwiseSpearman(X,Y);
end


function X=build_win(x,m)
% BUILD_WIN  Build windowed matrix from time series vector
%            Adapted from code originally written by David Meko
[mx,nx]=size(x);
if nx~=1;
    error('x must be a column vector');
end;

if m>=mx;
    error('window width must be shorter than time series!');
end;

% Make index to end times
i2 = fliplr(mx:-1:m);

% Expand to matrix
A = repmat(i2,m,1);
[mA,nA]=size(A);

% increment vector
b = (0:1:(m-1))';
B=repmat(b,1,nA);

C=A-B;
C=flipud(C);
X=x(C);


function s = pairwise(X,Y)
% PAIRWISE  Pearson correlations between paired columns of two matrices,
%           modified from David Meko's corrpair.m. Uses nanmean/nanstd (as
%           pairwiseSpearman below already does) for robustness to missing
%           data, and relies on implicit expansion (R2016b+) for the
%           mean-centering step rather than manually building repmat'd
%           mean matrices, which was the standard workaround before
%           implicit expansion existed.
[mX,nX]=size(X);
[mY,nY]=size(Y);
x_bar = nanmean(X);
y_bar = nanmean(Y);

% Matrices of std devs; these are computed with "N-1" in denominator
xstd = nanstd(X);
ystd = nanstd(Y);

% Matrices of departures (implicit expansion: x_bar/y_bar are 1-by-n row
% vectors, broadcast against the mX-by-n / mY-by-n data matrices)
Dx = X - x_bar;
Dy = Y - y_bar;

% Product of departures
D2 =  Dx .* Dy;

% Sample covariance
D = sum(D2)/(mX-1); % rv

% Std of x times std of y
xystd = xstd .* ystd; % a rv

% Correlation
s = D ./ xystd;


function s = pairwiseSpearman(X, Y)
% PAIRWISESPEARMAN  Spearman correlations between paired columns of two
%                    matrices. Uses implicit expansion (R2016b+) for
%                    mean-centering, same as pairwise() above.
[mX, nX] = size(X);
[mY, nY] = size(Y);

if mX ~= mY || nX ~= nY
    error('Matrices X and Y must have the same dimensions.');
end

% Rank the data
X_ranked = tiedrank(X);
Y_ranked = tiedrank(Y);

% Calculate mean of the ranks
x_bar = nanmean(X_ranked);
y_bar = nanmean(Y_ranked);

% Calculate the deviations of ranks from the mean ranks (implicit expansion)
Dx = X_ranked - x_bar;
Dy = Y_ranked - y_bar;

% Calculate the product of deviations
Dxy = Dx .* Dy;

% Sum up the products for covariance calculation
cov_xy = sum(Dxy) / (mX - 1);

% Calculate standard deviations of ranked data
xstd = nanstd(X_ranked);
ystd = nanstd(Y_ranked);

% Calculate Spearman correlation
s = cov_xy ./ (xstd .* ystd);


function [r,c,V] = findnearest(srchvalue,srcharray,bias)

% Usage:
% Find the nearest numerical value in an array to a search value
% All occurances are returned as array subscripts
%
% Output:
%
% For 2D matrix subscripts (r,c) use:
%
%       [r,c] = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% To also output the found value (V) use:
%
%       [r,c,V] = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% For single subscript (i) use:
%
%         i   = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% Inputs:
%
%    srchvalue = a numerical search value
%    srcharray = the array to be searched
%    bias      = 0 (default) for no bias
%                -1 to bias the output to lower values
%                 1 to bias the search to higher values
%                (in the latter cases if no values are found
%                 an empty array is ouput)
%
%
% By Tom Benson (2002)
% University College London
% t.benson@ucl.ac.uk

if nargin<2
    error('Need two inputs: Search value and search array')
elseif nargin<3
    bias = 0;
end

% find the differences
srcharray = srcharray-srchvalue;

if bias == -1   % only choose values <= to the search value

    srcharray(srcharray>0) =inf;

elseif bias == 1  % only choose values >= to the search value

    srcharray(srcharray<0) =inf;

end

% give the correct output
if nargout <= 1

    if all(isinf(srcharray(:)))
        r = [];
    else
        r = find(abs(srcharray)==min(abs(srcharray(:))));
    end

elseif nargout>1
    if all(isinf(srcharray(:)))
        r = [];c=[];
    else
        [r,c] = find(abs(srcharray)==min(abs(srcharray(:))));
    end

    if nargout==3
        V = srcharray(r,c)+srchvalue;
    end
end

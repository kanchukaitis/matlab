function [sd99,sd95,sd90,sd5,psig]=gsbtest(seriesLength,iters,window,seriesCorrelation,method,obs)

% GSBTEST  'Gershunov' Test of Moving Correlation Significance Levels
%              
%      [sd99,sd95,sd90,sd5,pcrit]=gsbtest(sL,sN,window,R,method,pval)
%
% Implementation of a Monte Carlo type simulation of moving correlation functions 
% between two weakly correlated white noise series for the purposes of constructing
% significance test for moving correlation functions based on their standard deviation (Gershunov et al. 2001).
%
% Inputs:
%   seriesLength = length of the time series used in the running correlation analysis
%   iters = number of Monte Carlo iterations to perform (should be >= 1,000)
%   window = width (in units of the time step for the data) of the correlation window
%   seriesCorrelation = desired overall correlation of the simulated time series (simulations will be within 1% of this value)
%   method [optional] = Which method to use to create correlated random series [default = 1], see below
%   obs [optional] = find the significance level for a specific (observed) value of the std of the moving correlation function
%
% Correlation Methods
%   1 = Cholsky decomposition of the covariance matrix [default]
%   2 = Matrix square root
%
% Outputs:
%  sd99  = bootstrapped 99% confidence level
%  sd95  = bootstrapped 95% confidence level
%  sd90  = bootstrapped 90% confidence level
%  sd5   = bootstrapped 5% confidence level  
%  psig  = significance for specific (observed) standard deviation (obs)
%

% Example:
%  > [sd99,sd95]=gsbtest(100,1000,11,0.30,2)
%
% sd99 =
%
%    0.4002
%
% sd95 =
%
%    0.3624
%           
% NOTE: Although the solution should converge for sN >= 1000, small differences (due to the correlation between
% the two artificial time series generated in the script), may still result, especially for a small number (<<10000)
% of simulations; so it is recommended that you verify convergence with several runs and be cautious of values too
% close to your chosen test statistic (particularly at 90% or 95% CI).  You can compare the above result for 'sd95'
% to Table 1 from Gershunov et al. 2001 (for interseries r = 0.30, window length = 11, their 95% confidence level = 0.36).
%
% References
% ----------
% Gershunov, A., N. Schneider, and T. Barnett, 2001. Low-frequency modulation of the ENSO-Indian 
% monsoon rainfall relationship: Signal or noise?, Journal of Climate, 14: 2486-2492.
%

% Update History:
%
% 09 13 -- minor code improvements for clarity; added 'findnearest' as subfunction
% 03.11 -- added calculation of 'psig' using input 'obs' 
% 06.08 -- posted online
% 10.07 -- initial public version
% 09.03 -- initial version
% 
% Function gsbtest.m written by Kevin Anchukaitis
% Subfunction build_win is a modification of the function by David Meko (dmeko@ltrr.arizona.edu)
% Subfunction pairwise is a modification of corrpair.m by David Meko (dmeko@ltrr.arizona.edu)
% Included subfunction 'findnearest' by Tom Benson

if nargin<6; obs = [];  end;
if nargin<5; method = 1; end;
if nargin<4; error('Too Few Input Parameters'); end;

if nargout==5 && isempty(obs)
    error('Need to provide desired p value to output significance level')
    return
end

% Create a correlation matrix from the desired correlation value for the white noise series
Rmat=[1 seriesCorrelation;seriesCorrelation 1];

randn('state',sum(100*clock));  % re-seed the random number generator (outdated)

% this should probably be vectorized, and other things, to improve speed
tol     = 0.005; % hardcoded tolerance for the simulated correlation
forgive = 1000;   % how many times will the loop run to get below tol?
if method==1 % Cholsky decomposition
    	for i=1:iters
            C=chol(Rmat);  % Cholsky decomposition of the correlation matrix
            rcheck = 1;
            count_repeated_checks = 0; 
             while abs(rcheck) > tol*seriesCorrelation; 
              Z=randn(2,seriesLength);
              X=C'*Z; X=X';
              check = corrcoef(X);
              rcheck = check(2,1) - seriesCorrelation;
              count_repeated_checks = count_repeated_checks + 1;
              if count_repeated_checks >= forgive
                  disp(['stuck on recheck = ',num2str(rcheck),', giving up ...'])
                  break
              end
             end 
            % now do the moving correlation function on the series in X
            x1=build_win(X(:,1),window); x2=build_win(X(:,2),window);
            s = pairwise(x1,x2);
            mcf_var(i) = var(s); 
            mcf_std(i) = std(s);
        end
		
 elseif method==2 % Matrix square root
		for i=1:iters
		    rcheck = 1;
             while abs(rcheck) > tol*seriesCorrelation; 
              X=(sqrtm(Rmat)*randn(2,seriesLength))';
              check = corrcoef(X);
              rcheck = check(2,1) - seriesCorrelation;
             end
            % now do the moving correlation function on the series in X
            x1=build_win(X(:,1),window); x2=build_win(X(:,2),window);
            s = pairwise(x1,x2);
            mcf_var(i) = var(s); 
            mcf_std(i) = std(s);
		end   

else 
error('Invalid Series Creation Method'); 
end
    
% sort the matrix of variance or standard deviation 
sorted_mcf_var = sort(mcf_var'); 
sorted_mcf_std = sort(mcf_std');

% get the desired values
sd99 = sorted_mcf_std(round(0.99 * seriesLength)); 
sd95 = sorted_mcf_std(round(0.95 * seriesLength)); 
sd90 = sorted_mcf_std(round(0.90 * seriesLength)); 
sd5 =  sorted_mcf_std(round(0.05 * seriesLength));   

% calculate the significance level for the observed, if passed as input
if ~isempty(obs)
 flag_pval = findnearest(obs,sorted_mcf_std);
 psig = (flag_pval)/length(sorted_mcf_std);
end



%% subfunctions

function X=build_win(x,m)
% Build windowed matrix from time series vector
% Adapted from code originally written by David Meko
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
% modified from David Meko's corrpair.m
[mX,nX]=size(X);
[mY,nY]=size(Y);
x_bar = mean(X);
X_bar = repmat(x_bar,mX,1);
y_bar = mean(Y);
Y_bar = repmat(y_bar,mY,1);

% Matrices of std devs; these are computed with "N-1" in denominator
xstd = std(X);
Xstd = repmat(xstd,mX,1);
ystd = std(Y);
Ystd = repmat(ystd,mY,1);

% Matrices of departures
Dx = X - X_bar;
Dy = Y - Y_bar;

% Matrices of squared departures
D2x = Dx .* Dx;
D2y = Dy .* Dy;

% Product of departures
D2 =  Dx .* Dy;

% Sample covariance
D = sum(D2)/(mX-1); % rv

% Std of x time std of y
xystd = xstd .* ystd; % a rv

% Correlation 
s = D ./ xystd;


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
if nargout==1 | nargout==0
    
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


    


function [output,newx,beta,debug] = quantile_mapping(rec,obs,nlls,fullSeries,nboot,span,prct)

% QUANTILE MAPPING  Bias correction of a time series compared to the corresponding observations using quantile mapping
%
%  [output,newx,beta,debug] = quantile_mapping(rec,obs,nlls,fullSeries,nboot,span,prct)
%
%
% The inputs to the function are:
%  rec = [n x 2] reconstruction, with time in the first column, and values in the second
%  obs = [p x 2] observations, with time in the first column, and values in the second
%  nlls = the number of samples to consider for the neighborhood - default is to use the nearest 10 neighbors
%  fullSeries - [m x 2] specify an alternative series to apply to mapping to for the output, default is to use the full span of the reconstruction e.g. rec(:,2)
%  nboot = number of iterations to run a bootstrap that samples (with replacement) from the reconstruction and observation
%  span = a vector of years (to be comparable to rec(:,1) and obs(:,1)) to use the fit the mapping
%  prct = specify a set of quantiles to use, defauly is from 0 to 1 by 0.01 (for 101 total steps)
%
% The outputs are:
%  output - this are the quantile mapping bias corrected values, will have a length of either length(rec(:,2) or length(fullSeries)
%  newx - these are the values at each quantile of the original series from rec(:,2)
%  beta - the weights from the local regression ~ the mapping and the slope of the local relationship
%  debug - structure for storing diagnostics and other things you might want to look at to make sure things are running as expected

% version information:
% This function is based on the R code from Gudmundsson (qmap) and the application to paleoclimate time series in Robeson et al (2020)
% v0.1 02/14/2021 - initial version by Kevin Anchukaitis

% Citations:
%
% Gudmundsson et al. (2012) Technical Note: Downscaling RCM precipitation to the station scale using statistical transformations - a comparison of methods. Hydrology and Earth System Sciences, 16, 3383-3390
% Teutschbein and Seibert (2012). Bias correction of regional climate model simulations for hydrological climate?change impact studies: Review and evaluation of different methods. Journal of Hydrology, 456, 12?29.
% Maraun (2013). Bias correction, quantile mapping, and downscaling: Revisiting the inflation issue. Journal of Climate, 26(6), 2137-2143.
% Cannon et al. (2015). Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes?. Journal of Climate, 28(17), 6938-6959.
% Robeson et al. (2020). Bias correction of paleoclimatic reconstructions: A new look at 1,200+ years of Upper Colorado River flow. Geophysical Research Letters, 47(1), e2019GL086689.
%

% Notes (2/14): It is possible (likely) that this function provides different output than the R code in qmap, most likely due to different behavior in the interp1() in MATLAB vs
% approx() from R.

% a little bit of error checking
if nargin <2
    disp('At least rec and obs must be provided to run this code')
    return
end

% if size(rec,2)~= 2 | size(obs,2) ~= 2
%     disp('Size of rec or obs not as expected. Please be sure you provide a 2 column matrix with time in the first column and the data in the second')
%     return
% end

% now set some defaults
if nargin < 7
    prct = 0:0.01:1.0; % quantiles we'll evaluate for local regression
end

if nargin < 6 % automated common interval setting
    [commonYear,indx,jndx] = intersect(rec(:,1),obs(:,1));
    span = commonYear;
    full_rec = rec(:,2:end);
    rec = rec(indx,:);
    obs = obs(jndx,:);
else
    full_rec = rec(:,2:end);
    [~,~,indx] = intersect(span,rec(:,1));
    rec = rec(indx,:);
    [~,~,jndx] = intersect(span,obs(:,1));
    obs = obs(jndx,:);
end

if nargin < 5
    nboot = 1; % no bootstrapping by default
end

if nargin<3
    nlls = 10; % by default use 10 values for the neighborhood, ~10%
end

%% start quantile mapping operations here
% extract values from input
y = obs(:,2:end);
x = rec(:,2:end);

% sort the observations and reconstructions
ys = sort(y); % observations
xs = sort(x); % reconstruction

% these are the values at each quantile in the reconstruction
newx = quantile(xs,prct); % value at every quantile in prct

% preallocate the beta matrix
beta = NaN(size(newx,1),size(newx,2),2,nboot);

for j = 1:nboot % loop over bootstrap iterationss
    
    for i = 1:length(newx) % loop over quantile steps
        
        for k = 1:size(newx,2)
            
            if j==1 % no bootstrap or first bootstrap, use series as normal
                xss = sort(xs(:,k));
                yss = sort(ys(:,k));
            else % bootstrap, draw with replacement from reconstruction                
                [xss,idx] = datasample(xs(:,k),length(xs),'Replace',true);
                yss = ys(idx,k); % associated values from the observations
                xss = sort(xss); % re-sort for safety
                yss = sort(yss); % re-sort for safety
            end
            
            % this is a straight port of the R code from qmap, one way of identifying the k neighbors for the local regression
            xc = xss - newx(i,k); mdist = sort(abs(xc)); mdist = mdist(nlls); ik = abs(xc) <= mdist;
            xcs = [ones(length(xc(ik)),1) xc(ik)]; % these are the values in the neighborhood plus a column of ones
            beta(i,k,1:2,j) = regress(yss(ik),xcs); % here, the regression of the sorted local values gives us an estimate of x in beta(:,1), a slope in beta(:,2)
        end
    end
end

% average over the bootstrapped values
beta = squeeze(nanmean(beta,4));

if nargin<4
    fullSeries = full_rec; % you can use the series you put in for the fit, or you can use another series for the estimate
end

for k = 1:size(newx,2)
    output(:,k) = interp1(newx(:,k),beta(:,k,1),fullSeries,'linear','extrap');
end
% for testing - extrapolation can actually miss extreme values in the paleoclimate record - in theory this could correct for that
%outsideBounds = (output < min(newx) | output > max(newx));
%output(outsideBounds) = fullSeries(outsideBounds);


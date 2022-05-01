%% REDFIT
% Following the REDFIT algorithm of Mudelsee 2002, with code modified from
% the original FORTRAN 90 program, the dlpR R package and the astroChron R
% package
%
% References:
% *REDFIT*:
% Authors of REDFIT
% Authors: Michael Schulz, MARUM and Faculty of Geosciences, Univ. Bremen
% -------- Klagenfurter Str., D-28334 Bremen
%          mschulz@marum.de
%          www.geo.uni-bremen.de/~mschulz
%
%          Manfred Mudelsee, Inst. of Meteorology, Univ. Leipzig
%          Stephanstr. 3, D-04103 Leipzig
%          Mudelsee@rz.uni-leipzig.de
% Reference: Schulz, M. and Mudelsee, M. (2002) REDFIT: Estimating
% ---------- red-noise spectra directly from unevenly spaced paleoclimatic
%            time series. Computers and Geosciences, 28, 421-426.
% *dlpR:*
% Andy Bunn, Mikko Korpela, Franco Biondi, Filipe Campelo, Pierre M ?erian,
% Fares Qeadan, Christian Zang, Darwin Pucha-Cofrep and Jakob Wernicke
% (2018). dplR: Dendrochronology Program Library in R. R package version
% 1.6.7. https://r-forge.r-project.org/projects/dplr/
%
%astroChron and 'Seeing Red':
% Meyers, S.R. (2014). Astrochron: An R Package for Astrochronology.
% http://cran.r-project.org/package=astrochron

%% Preamble
%[~,n]=unique(t); %ensure t has no repeating values
%t=t(n); 
%x=x(n);
x=irdmar(~isnan(irdmar));
t=age(~isnan(irdmar));
%x=detrend(x); %remove first order linear trend from x
npts=length(x); %number of points
nsim=1000; %number of simulations, default 1000 but higher for publication / thesis
nseg = length(t);

%% window weights 
% Welch 1 - commented out as it is not currently used in the code, possibly
% integrate at a later time, for now, ignore
% nseg = length(t);
% tlen=t(nseg)-t(1);
% tlenfull=nseg*tlen/(nseg-1);
% tpeak=t(nseg)-tlenfull/2;
% ww = (t - tpeak) / (t(nseg) - tpeak);
% ww = 1 - ww * ww;
% ww = ww * sqrt(nseg / sum(ww * ww)); %scale window weights

%% get spectra for data (x)
%firstly, estimate gxx, the lomb-scargle power spectrum for x
ofac = 4; %oversampling factor

[gxx,f]=plomb(x,t,[],ofac,'normalized'); 
varx = mean(diff(f)) * sum(gxx); %estimate area under gxx, to get data variance idea

%% TAUEST
%persistence estimation for unevenly spaced time-series, from Manfred
%Mudelsee
n50=3; %number of segements for tauestimation
nseg = round(npts / (n50 + 1) * 2); %indexing
segskip = (npts-nseg)/(n50-1); %more indexing
rhovec = n50; %pointless renaming of variable
    for i = 1:length(n50)
        if i == 1
        iseg = [1,nseg];
        end
        if i == 2
                iseg = [nseg,nseg+segskip];
        end
        if i == 3
            iseg = [nseg+segskip,length(t)];
        end
        twk = t(iseg(:,1):iseg(:,2)); 
        twkM(:,i) = twk;
        xwk = x(iseg(:,1):iseg(:,2));
        %detrend segment
        xwk = detrend(xwk); 
       
    % estimate rho for each segment
    xscal=xwk/std(xwk); %scale segment
    dt=mean(diff(twk));
    lag0=xscal(1:end-1);
    lag1=xscal(2:end);
    rho=corr(lag0,lag1); %find autocorrelation for each segment
    %now correct and average
    rhovec(i)=(rho*(nseg-1)+1)/(nseg-4); %bias correction for rho (Kendall & Stuart, 1967; Vol. 3))
    end
rho = mean(rhovec);
tau = -dt/log(rho); 
%we don't actually need tau, thanks to the dlpR workaround for finding rho,
%but it's including to check that all is well, i.e. positive value

%% generate red noise spectrum
sdev=std(x);
dt=mean(diff(t));
%initialize matrix for sim results
red=zeros(npts,nsim);
%dim(red) = c(npts,nsim);

for i = 1:nsim
    %generate normal deviates 
    white=randn(npts,1);
% generate AR1 red noise
red(1,i)=white(1);
for ii=2:npts
    red(ii,i)=rho*red(ii-1,i)+white(ii,1);
end

    red(:,i)=red(:,i)-mean(red(:,i));
    red(:,i)=red(:,i)*sdev/std(red(:,i));
    red(:,i)=red(:,i)+mean(x);
end
for i= 1:nsim
    [grx(i,:),ff]=plomb(red(:,i),t,[],4,'normalized'); %compute spectra for red noise
    varr1 = mean(diff(ff))*sum(grx(i,:)); %find variance and scale
    grx(i,:)=varx/varr1*grx(i,:);
end

redav=mean(grx,1); %average spectra
varr2 = mean(diff(ff)) * sum(redav);
redav = varx / varr2 * redav; %scale to match area under curve of gxx
rhosq=rho*rho;
%set theoretical spectrum (e.g., Mann and Lees, 1996, Eq. 4)
gredth=(1-rhosq) ./ (1+rhosq-2*rho*cos(0:pi/(length(f)-1):pi));
cor = redav/gredth;
invcorr = gredth / redav;
 gxxc = gxx*invcorr;
 
%% Red noise confidence levels
%false-alarm levels from percentiles of MC simulation
        ci99=prctile(grx,99)*invcorr;
         ci995=prctile(grx,99.5)*invcorr;
          ci95=prctile(grx,95)*invcorr;
          ci90=prctile(grx,90)*invcorr;
%% plot
plot(f,smooth(gxx,15),'color','k','linewidth',2);
hold on
plot(f,smooth(ci95,length(ci95)/2),'linestyle','--')
plot(f,smooth(ci90,length(ci95)/2),'linestyle','--')

% legend('?^1^3C Spectrum','99.9th red noise percentile','95th red noise percentile','90th red noise percentile','location','best')
% legend('boxoff')
xlim([0 0.15])
xlabel('Frequency (kyr^-^1)');
ylabel('Normalized Power')
xticks([0.0025 0.01 0.025 1/23])
xticklabels([1/0.0025, 1/0.01, 1/0.025, 23]);
xtickangle(45)
legend('IRD MAR','95% CI','90% CI')
box off
set(gca,'tickdir','out','color','none')
title('Lomb-Scargle Periodogram with AR(1) noise model')
set(gca,'fontsize',13)

%% interp and stats
fs=prctile(diff(age),75);
t=[1:fs:age(end)];
x1=interp1(age(~isnan(irdmar)),irdmar(~isnan(irdmar)),t);
t=t(~isnan(x1));
x1=x1(~isnan(x1));
[P,s,ci]=pmtmPH(x1,fs,3,1,length(x1));

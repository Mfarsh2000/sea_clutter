function [gammamatrix,gaussmatrix,CurrentScan]=gammageneratortime(Radar,gaussmatrix)
% Clutter Generator
% Masoud Farshchian, NRL 2010
% Clutter Params

SS=Radar.Clutter.SS;  % Sea State - Douglass
height=Radar.Tx.Height; % Height in Meters
Ts=Radar.Clutter.CorrelationTime;  % Correlation Time in secnds
azBW=Radar.Tx.BEAM_WIDTH;  % Beam Width in Degrees


NumbPRIs=Radar.Simulation.NumbPRIs;  % number of PRIs per Scan to simulate
RangeCells=Radar.Simulation.RangeCells; % number of Range Cells
RPM=Radar.Simulation.RPM;

polarization=Radar.Tx.Polarization;  % Polarization
rangeres=Radar.Tx.res;
Az_Filter=Radar.Tx.AzFilter;
numbrangecells=Radar.Simulation.RangeCells;
frequencyagility=Radar.Tx.frequencyagility;


Aw=Radar.Clutter.Aw;
samplerange=Radar.Target.SampleRange; % nm 


Ts=RPM/60*Ts;  % convert time correlation constant to scan correlation constant
Az_Filter=Az_Filter./norm(Az_Filter);
AzFilterSize=max(size(Az_Filter));

if SS==0  % Sea State Zero - No range correlation
    Tr=0;
else

    W=[2.5 4.5 6 8.5 11 14];
    
    gravity=9.8;
    Tr=(pi/2)*(W(SS)^2/gravity)*(3*((cos(Aw*pi/180)).^2)+1).^(0.5);
    Tr=Tr/rangeres;             % Correlation Range
    % Cell-average CFAR gain in spatially
    % correlated K-distributed Clutter
    % S. Watts
    % choose the appropriate correlation size
    % in cells
end
%Tr=repmat(Tr.',1,max(size(RangeCells)));

alpha=shapefactor(polarization,samplerange,height,rangeres,azBW,Aw);
alpha=alpha*ones(NumbPRIs,numbrangecells);
T=(1+.15./alpha.^(0.7));

% range
a = exp(-1./(Tr.*T));
b = sqrt(1-a.^2);

Ts=abs(Ts*cos(Aw*pi/180));  % time correlation 
%Ts=repmat(Ts.',1,max(size(RangeCells)));
% time
A = exp(-1./(Ts.*T));
B = sqrt(1-A.^2);


if  isempty(gaussmatrix)
    gaussmatrix=randn(NumbPRIs,numbrangecells);
    for i=2:numbrangecells, gaussmatrix(:,i) = a(:,i-1).*gaussmatrix(:,(i-1))+b(:,i-1).*randn(NumbPRIs,1); end    % range
else
    gaussnew=randn(NumbPRIs,numbrangecells);
    for i=2:numbrangecells, gaussnew(:,i) = a(:,i-1).*gaussnew(:,(i-1))+b(:,i-1).*randn(NumbPRIs,1); end % range
    gaussmatrix=A.*gaussmatrix+B.*gaussnew;     % time
end


GaussGamma=filter(Az_Filter,1,gaussmatrix);


%gammamatrix = rt_gaminv(0.5 *erfc(-GaussGamma / sqrt(2)),1,alpha);
%gammamatrix = gaminv(0.5 *erfc(-GaussGamma / sqrt(2)),alpha,1./alpha);
gammamatrix = icdf('gam',0.5 *erfc(-GaussGamma / sqrt(2)),alpha,1./alpha);
% xgam = RT_acorr(gam_noise-1,lag_time); % auto correlation of Gamma
gammamatrix((isnan(gammamatrix)))=0;  % R-T Toolbox breaks down for very small nu's

if frequencyagility
   specklecorr=1/sqrt(2)*(randn(NumbPRIs,numbrangecells)+sqrt(-1)*randn(NumbPRIs,numbrangecells));
else
    speckfilterAz=ones(1,AzFilterSize);
    speckfilterAz=speckfilterAz./norm(speckfilterAz);
    %speckfilterR=ones(1,RangeCells);
    %speckfilterR=speckfilterR./norm(speckfilterR);
    specklecorr=1/sqrt(2)*(randn(NumbPRIs,numbrangecells)+sqrt(-1)*randn(NumbPRIs,numbrangecells));
    %specklecorr=filter(speckfilterAz,speckfilterR,specklecorr);
    specklecorr=filter(speckfilterAz,1,specklecorr);
end
% 
 % AdditiveNoise=1/sqrt(2)*(randn(NumbPRIs,numbrangecells)+sqrt(-1)*randn(NumbPRIs,numbrangecells));
% 
% CurrentScan=sqrt(1/(1+CNR)).*AdditiveNoise+sqrt(CNR/(1+CNR)).*sqrt(gammamatrix).*specklecorr;

%CurrentScan=sqrt(Pc).*sqrt(gammamatrix).*specklecorr;
CurrentScan=sqrt(gammamatrix).*specklecorr;


function  [vH] =shapefactor(polarization,range_cells,height,res,BW_deg,Aw)
%
% Purpose: Calculates K-distribution shape parameter
%
% Input: polarization 0 for VV and 1 for HH
%        range_cells = the  range cells in(nm)
%        height  = platform height (meters)
%        res, BW = radar resolution (meter, deg)
%        BW_deg = the azimuth Bandwidth
%        theta_value = look angle
%        
% comments from previous version:
% For up or down swell direction, a_wind = 0, 180: f_wind = -1/3
% For across direction, a_wind = 90: f_wind = 1/3
% For no swell or a_wind = 45: f_wind = 0
% Masoud Farshchian
% see : Ward, Tough, Watts "Sea Clutter: Scattering, the K Distribution and
% Radar Performance", IET , 2006. pg 237


d2r = pi/180;


%f_wind = -cos(2*(theta_value-az_wind))/3;
f_wind=-cos(2*Aw*pi/180)/3;  %upwind/crosswind only for now

% --- CLUTTER PATCH: AREA and GRAZING ANGLE -------------------------
beamwidth = BW_deg * d2r;              % rad



% create a PRI_SCANxrangecells matrix
r=range_cells*1852;  % range cells in meters


re = 8495124;                          % m  (4/3 Earth radius)
h = height;

angle=asin(h./r+h^2./(2*re.*r)-r./(2*re));
area = beamwidth.*r*res.*sec(angle); % m2

f_wind=repmat(f_wind.',1,max(size(range_cells)));
area=repmat(area,max(size(Aw)),1);
angle=repmat(angle,max(size(Aw)),1);

if polarization == 1
    vH = 10 .^( 2/3*log10(angle*180/pi) + 5/8*log10(area) - 2.09 + f_wind);
else
    vH = 10 .^( 2/3*log10(angle*180/pi) + 5/8*log10(area) - 1.39 + f_wind);
end

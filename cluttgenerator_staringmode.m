% Clutter generatore R-T

%Simulation Parameters
Radar.Simulation.NumbPRIs=256; % numb PRIs per Scan per the azimuth sector
Radar.Simulation.RangeCells=256; % number of range cells to simulatie

% Transmitter Properties
Radar.Tx.Fc=9.5e9;   % Transmit Frequency
Radar.Tx.Lambda=3e8/Radar.Tx.Fc; % Wavelength
Radar.Tx.Polarization=1;  % Polarization 0 (ver) 1 (Horizontal)
Radar.Tx.Height=100; % Height in Meters
Radar.Tx.PRF=1800;  % PRF 
Radar.Tx.BW=400e6; % Bandwidth
Radar.Tx.res=3e8/(2*Radar.Tx.BW); %Range Resolution - For now assume Chirp
Radar.Tx.Pulsewidth=40e-6;
Radar.Tx.frequencyagility=0;  %0 off 1 on
Radar.Tx.BEAM_WIDTH=2;  % Azimuth Beam Width in Degrees

Radar.Tx.frequencyagility=0;  %0 off 1 on

Radar.Target.SampleRange=6;  %nm


Radar.Tx.PRF=2880;

% Clutter Parameters
Radar.Clutter.SS=3;  % Sea State - Douglass
Radar.Clutter.Aw=0;     % 0 upwind, 90 cross-wind, 180 downwind  (deg)
Radar.Clutter.CorrelationTime=2;  % Correlation Time in seconds (Upwind) for Underlying mean
Radar.Clutter.CorrelationTimeSpeckle=10e-3;  % Correlation Time for speckle set at 10 milli-second per experimental data


% Power CNR
CNR=30;  % CNR in dB
CNRv=10^(CNR/10);
% Initialization varitables
gammamatrix=[];  % For Clutter Correlation
gaussmatrix=[];  % For Clutter Correlation
CurrentScan=[];  % Store the data for all Scans before PP

    
      
%      % generate Clutter
    [CurrentScan]=gammageneratortime_staringmode(Radar);
    
%      % add Noise
      AdditiveNoise=1/sqrt(2)*(randn(Radar.Simulation.NumbPRIs,Radar.Simulation.RangeCells)+sqrt(-1)*randn(Radar.Simulation.NumbPRIs,Radar.Simulation.RangeCells));
      CurrentScan=sqrt(1/(1+CNRv)).*AdditiveNoise+sqrt(CNRv/(1+CNRv)).*CurrentScan;
    
      TotScans=[CurrentScan];    
      
TotScans=TotScans';
     
imagesc(abs(TotScans))
xlabel('time')
ylabel('range')
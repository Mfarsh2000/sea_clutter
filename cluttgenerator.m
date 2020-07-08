% Clutter generatore R-T

%Simulation Parameters
Radar.Simulation.NumbPRIs=256; % numb PRIs per Scan per the azimuth sector
Radar.Simulation.RangeCells=256; % number of range cells to simulatie
Radar.Simulation.RPM=16;  % 16 rev per minute
Radar.Simulation.numbscans= 16%  16 scans to make 256x256
% Transmitter Properties
Radar.Tx.Fc=9.5e9;   % Transmit Frequency
Radar.Tx.Lambda=3e8/Radar.Tx.Fc; % Wavelength
Radar.Tx.Polarization=1;  % Polarization 0 (ver) 1 (Horizontal)
Radar.Tx.Height=100; % Height in Meters
Radar.Tx.PRF=1800;  % PRF 
Radar.Tx.BW=25e6; % Bandwidth
Radar.Tx.res=3e8/(2*Radar.Tx.BW); %Range Resolution - For now assume Chirp
Radar.Tx.Pulsewidth=40e-6;
Radar.Tx.frequencyagility=0;  %0 off 1 on
Radar.Tx.BEAM_WIDTH=2;  % Azimuth Beam Width in Degrees
Radar.Tx.frequencyagility=0;  %0 off 1 on
Radar.Target.SampleRange=6;  %nm
Radar.Tx.PRF=2880;
Radar.Tx.AzFilter=ones(1,16);  % Azimuth beamwidth shape  % Could use Gaussian or Sinc^2 beamwidth


%Radar.Tx.AzFilter=ones(1,16);  % Azimuth beamwidth shape  % Could use Gaussian or Sinc^2 beamwidth
% Clutter Parameters
Radar.Clutter.SS=3;  % Sea State - Douglass
Radar.Clutter.Aw=0;     % 0 upwind, 90 cross-wind, 180 downwind  (deg)
Radar.Clutter.CorrelationTime=2;  % Correlation Time in seconds (Upwind)

% Power CNR
CNR=30;  % CNR in dB
CNRv=10^(CNR/10);
% Initialization varitables
gammamatrix=[];  % For Clutter Correlation
gaussmatrix=[];  % For Clutter Correlation
CurrentScan=[];  % Store the data for all Scans before PP
TotScans=[];

 for ScanNumber=1:Radar.Simulation.numbscans
%      
      
%      % generate Clutter
      [gammamatrix,gaussmatrix,CurrentScan]=gammageneratortime(Radar,gaussmatrix);
%      % add Noise
      AdditiveNoise=1/sqrt(2)*(randn(Radar.Simulation.NumbPRIs,Radar.Simulation.RangeCells)+sqrt(-1)*randn(Radar.Simulation.NumbPRIs,Radar.Simulation.RangeCells));
      CurrentScan=sqrt(1/(1+CNRv)).*AdditiveNoise+sqrt(CNRv/(1+CNRv)).*CurrentScan;
    
      TotScans=[TotScans; CurrentScan];    
       
 end
      
TotScans=TotScans';
     
imagesc(abs(TotScans))
xlabel('time')
ylabel('range')

%***************************************************************************************************************
%                                 PROGRAMME FOR Plate Media
%***************************************************************************************************************
% MAIN BODY OF THE PROGRAMME TO CALCULATE SOURCE STRENGTH OF EACH LAYER OF SOURCES INCLUDING TRANSDUCER
% THIS IS TO CALCULATE THE PRESSURE DISTRIBUTION IN MULTILAYER SOLIDS DUE TO TRANSDUCERS KEPT INCLINED TO
% INTERFACES IN THE FLUID. IN THIS PROGRAMME 2 FLUIDS AND 2 SOLIDS HAS BEEN CONSIDERED BUT FOR MORE SOLIDS PROGRAMME                                                                                                    *
% RUNNER NEEDS TO ADD FEW LINES IN THE "Calculation Of M Matrix" PART OF THE PROGRAMME ACCORDINGLY BY SEEING
% THE FORMULATION WITH SOUND UNDERSTANDING AND NEED TO CHANGE THE CONSECUENCES AFTERWARDS.
%
% SPECIAL NOTE: INCLUDE SUBPROGRAMMES
% 1. DiscretizeInterfaces.m      % TO GET THE GEOMETRY OF SOURCES SINGLE POINT SOURCE AT FLUID-SOLID INTERFACE
% 2. Christoffel_Sol.m           % TO GET THE PHASE VELOCITIES AND VECTORS BY SOLVING CHRISTOFFEL'S EQUATION
% 3. GeomSensCircl.m             % FOR GEOMETRY OF SENSOR
% 4. PressureFluidMatt.m         % PRESSURE VELOCITY MATRIX CALCULATOR IN FLUID
% 5. DisFluidmatt.m              % DISPLACEMENT CALCULATOR IN FLUID
% 6. Aniso_Green.m               % STRESS  AND DISPLACEMENT CALCULATOR IN SOLID
% 7. CompPressureFluid.m         % SUB TO SUBPROGRAMME - EQUATIONS TO CALCULATE PRESSURE GREEN FUNCTION IN FLUID
% 8. CompWaveFieldSolid.m        % SUB TO SUBPROGRAMME - EQUATIONS TO CALCULATE STRESS GREEN FUNCTION AND DISP GREEN FUNCTION IN SOLID
% 9. solid_green.m               % SUBPROGRAM TO CALCULATE DISPLACEMENT AND STRESS GREEN'S FUNCTION AT ANY POINT IN SOLID
% 10. SphereChristofel.m         % SUBPROGRAM TO CALCULATE THE PHASE VELOCITIES AND VECTORS BY SOLVING CHRISTOFFEL'S EQUATION
%
%***************************************************************************************************************
%%

tic;

format long;
%%  INPUT DATA

%***************************************************************************************************************
%% Geometrical data

NumSolidLay=1;                                                   % Number of solid layers in multilayer solids
NumFluidLay=2;                                                   % Number of Fluid layers on both side of solid faces
NumIntrFc=NumSolidLay+1;                                         % Number of interfaces
NumSolidFluidIntrFc=2;                                           % Number of fluid solid interface
NumTrans=2;                                                      % Number of Transducers

NumSourcePt_IntrFc_x=115;%115;%                                 % Number source points (on either side) on the surface along X axis preferably odd number
NumSourcePt_IntrFc_y=3;                                         % Number source points (on either side) on the surface along Y axis preferably odd number
Length_IntrFc_x=20;%20;%                                        % Length of the interface along X axis%20;%
Length_IntrFc_y=(Length_IntrFc_x*NumSourcePt_IntrFc_y)/NumSourcePt_IntrFc_x; % Length of the interface along Y axis

% NumSourcePt_IntrFc_x=1;                                 % Number source points (on either side) on the surface along X axis preferably odd number
% NumSourcePt_IntrFc_y=115;                                     % Number source points (on either side) on the surface along Y axis preferably odd number
% Length_IntrFc_y=20;                                             % Length of the interface along X axis
% Length_IntrFc_x=(Length_IntrFc_y*NumSourcePt_IntrFc_x)/NumSourcePt_IntrFc_y;      % Length of the interface along Y axis

%% Transducer Input Parameters

InnerR_Trans=0;                          % Inner radius of transducer
OuterR_Trans=2;%3.5;                          % Outer radius of transducer
NumSourcePt_Trans=100;                   % Number of transducer source point

% Origin of the transducer
% TransCoord_z=[0,0,0];                  % Coordinate of the centers of the transducers in interface axis system,

% shifting the origin of the transducer
IntrFcShift=0;%.50;                           % Percentage of shift of interface coordinate towards left.

% Rotation of the trasnsducers (anticlockwise is positive)
Rotation_Trans1=0;%8.9743;%30;%50.7377;%59.1701;%0;%52.4241;%59.37;                       % Angle of rotattion of transducer 1, i.e inclination of the transducer anticlockwise
Rotation_Trans2=0;%-8.9743;%30;%50.7377;%59.1701;%0;%-52.4241;%59.37;                      % Angle of rotattion of transducer 2, i.e inclination of the transducer clockwise

Dist_IntrFc1=5;                                            % Distance between transducer and Interface 1 (fluid-solid interface)
Dist_IntrFc2=8;                                            % Distance between transducer and Interface 2 (solid-solid interface)

Dist_2ndTrans=Dist_IntrFc2+Dist_IntrFc1;                   % Distance of the 2nd transducer from the first one

ReceiverShift=-2*(IntrFcShift/100)*Length_IntrFc_x;         % Shift of Receiver transducer along X axis
Dist_IntrFc=[0 Dist_IntrFc1 Dist_IntrFc2 Dist_2ndTrans];    % Distance between transducer and Interfaces

%% Frequency Signal Input from the experiment

% Time Domain Data
% f=[2];% 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5];                        % Actuation frequency in MHz
% FtSignal = [2.25]; % Frequency Amplitude Signal in time domain (MHz)
SampRate = 10;        % Sampling Rate (MS/sec)
delTime = 1/(SampRate); % dt for time domain in micro sec.
NumSampPt = 256;   % No. of sample point

% Pre allocation of memory
% TimeStamp=zeros(NumSampPt,1);
% ForceTimeSignal=zeros(NumSampPt,1);
signal_tb=zeros(NumSampPt,1);         % Signal with a single frequency content
% Signal=zeros(NumSampPt,2);

F=1;%(f_start:f_increment:f_end);  % in MHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Signal data
H=1e-3;%str2num(get(findobj('tag','H'),'string'));                     % Unit is in m.
% Signal Parameters
k=5000;%str2num(get(findobj('tag','k'),'string'));                     % A signal shape factor
NumCycles=5;%str2num(get(findobj('tag','NumCycles'),'string'));                % Number of cycles to be generated
% Central Frequency in rad / sec.

CentFreq=F;
OmegaC=2*pi*CentFreq; % Central Frequency in MHz
Tao=NumCycles/CentFreq;      % Total band of pulse
Tao0=Tao/2;
p1=(k*H*CentFreq/NumCycles)^2; % A factor
TimeStamp=(0:NumSampPt-1)*delTime;
ForceTimeSignal=exp(-p1.*((TimeStamp-Tao0).^2)./2).*sin(OmegaC.*TimeStamp);

% Signal(:,1)=TimeStamp;
% Signal(:,2)=ForceTimeSignal;

% Frequency Transform of the Signal
FwSignal = fft(ForceTimeSignal); % Frequency Amplitude Signal in frequency domain
delW = (0:NumSampPt-1)/(NumSampPt*delTime); % dw for frequency domain

FreqStamp=delW(1:NumSampPt/2); % Store half of the generated 'Actuation frequency' vector due to symmetricity about Nyquist Freq
FreqAmp=FwSignal(1:NumSampPt/2); % Store half of the generated 'Amplitude' vector due to symmetricity about Nyquist Freq

%% Input Signal Plot
%  figure(1)
% plot(Signal(:,1),Signal(:,2));
% plot(TimeStamp,ForceTimeSignal);
%  figure(2)
% plot(f,abs(frAmp));

%%                            Discretization CALCULATION
%*************************************************************************************************************************
NumIntfcSourceTot=NumSourcePt_IntrFc_x*NumSourcePt_IntrFc_y;   % Total no. of sources
IntrFcArea=Length_IntrFc_x*Length_IntrFc_y;               % Area of the interface
Source_EqivR=(IntrFcArea/(NumIntfcSourceTot*2*pi))^0.5;	      % Equivalent radius of the sources
DistSource_x=Length_IntrFc_x/NumSourcePt_IntrFc_x;        % Distance between two point sources along x-axis
DistSource_y=Length_IntrFc_y/NumSourcePt_IntrFc_y;        % Distance between two point sources along x-axis
IntrFcShift=-(IntrFcShift/100)*Length_IntrFc_x;           % Percentage of shift of interface coordinate towards right
IntrFcCoord_z=[0,0,Dist_IntrFc1;0,0,Dist_IntrFc2];        % Include coordinate of all interfaces
TransCoord_z=[0,0,0;0,0,Dist_2ndTrans];                   % Z-coordinate of the centers of the transducers
%% Discretization of Transducers and the Interface

% Calling of subprogramme DiscretizeInterfaces and GeomSensCircl to get point source coordinates for both interfaces and transducers

[IntrFcCoord_Cent,IntrFcCoord_Top,IntrFcCoord_Btm,Sw_IntrFcCoord_Cent,Sw_IntrFcCoord_Top,Sw_IntrFcCoord_Btm]...
    =DiscretizeInterfaces(Length_IntrFc_x,Length_IntrFc_y,NumSourcePt_IntrFc_x,NumSourcePt_IntrFc_y,NumSolidFluidIntrFc,IntrFcCoord_z,IntrFcShift);

[TransCoord_Cent,TransCoord_Top,TransCoord_Btm,NumSourcePt_Trans]=DiscretizeTransducer(InnerR_Trans,OuterR_Trans,NumSourcePt_Trans,NumTrans,TransCoord_z,Rotation_Trans1,Rotation_Trans2,Dist_2ndTrans);

%% Line Discretization along which Time Domain PE Signal to be plotted
ZLine=(Dist_IntrFc1+(Dist_IntrFc2-Dist_IntrFc1)):0.125:Dist_2ndTrans;
xyAxis = zeros(1,size(ZLine,2));

PE_TargetLine = [xyAxis; xyAxis; ZLine]; % keeping zero in X and Y but discretize Z

% NumSourcePt_Trans_New=TransCoord_x_Cent(1,:,1);   % (NumSourcePt_Trans_New?)

%% Wave Field Computation Parameters
WaveFieldMode=3;                       % Mode of Scan for wave field
% 1 = Point Calculation A-Scan
% 2 = Line Calculation B-Scan
% 3 = Plane calculation C-Scan

PlotMode=3;                            % Mode of Plotting for wave field
% 1 = XY plane plot
% 2 = YZ plane plot
% 3 = XZ plane plot
% 4 = All plane plot

%% Parameters required for Plotting the Results
% for WaveFieldMode=3
NumTarget_x=NumSourcePt_IntrFc_x;          % Number of points taken As target (to plot) along X axis must be odd number
NumTarget_y=NumSourcePt_IntrFc_y;          % Number of points taken As target (to plot) along Y axis must be odd number
NumTarget_z=31;                            % Number of points taken As target (to plot) along Z axis must be odd number
X_PlaneCoord=0.0;                          % Coordinate of X plane where surface plot is intended to find out
Y_PlaneCoord=0.0;                          % Coordinate of Y plane where surface plot is intended to find out
Z_PlaneCoord=Dist_IntrFc1+0.1;                          % Coordinate of Z plane where surface plot is intended to find out

% For calulating the wave field at a specific Point for WaveFieldMode=1
% TrgtX=10;                          % x-coordinate for A scan
% TrgtY=0;                           % y-coordinate for A scan
% TrgtZ=0;                           % z-coordinate for A scan

% For calculating the wave field on a line joining the given points for WaveFieldMode=2;
% StartPt=[5 6 0];                   % Start point for B scan
% EndPt=[0 0 3];                     % End point for B scan

%% Material Properties of Fluid
Fluid_rho=1;                       % Density of fluid layers
WaveVel_P=1.48;                    % P-wave speed


%% Material Properties of Solid [GPa] : Constitutive stiffness matrix

% Material: Transversely Isotropic
C = [143.8 6.2 6.2 0 0 0; 6.2 13.3 6.5 0 0 0; 6.2 6.5 13.3 0 0 0; 0 0 0 3.4 0 0; 0 0 0 0 5.7 0; 0 0 0 0 0 5.7];
Solid_rho=1.56;                    % Density of solid layers

% Fully Orthtropic
% C = [70 23.9 6.2 0 0 0;23.9 33 6.8 0 0 0;6.2 6.8 14.7 0 0 0;0 0 0 4.2 0 0;0 0 0 0 4.7 0;0 0 0 0 0 21.9];
% Solid_rho=1.5;                    % Density of solid layers

% Monoclinic
% C = [102.6 24.1 6.3 0 0 40;24.1 18.7 6.4 0 0 10;6.3 6.4 13.3 0 0 -0.1;0 0 0 3.8 0.9 0;0 0 0 0.9 5.3 0;40 10 -0.1 0 0 23.6];
% Solid_rho=1.56;
%% Christoffel Solution Parameters

AngTestPt=15;                       % Angle in degrees for discretization of Sphere
dTheta = AngTestPt;
dPhi = AngTestPt;

%% End of Input Values
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Transducer Boundary Condition
Vso=0;                         % Velocity of transducer (Source or Transmitter) face enter magnitude only BOTTOM Transducer
Vto=1;                         % Velocity of transducer (Source or Receiver) face enter magnitude only

%% Initialization
% Source Strengths
As = zeros(NumSourcePt_Trans,size(FreqStamp,2));
A1 = zeros(NumSourcePt_IntrFc_x*NumSourcePt_IntrFc_y,size(FreqStamp,2));
A_1=zeros(NumSourcePt_IntrFc_x*NumSourcePt_IntrFc_y*3,size(FreqStamp,2));
A2=zeros(NumSourcePt_IntrFc_x*NumSourcePt_IntrFc_y*3,size(FreqStamp,2));
A_2 = zeros(NumSourcePt_IntrFc_x*NumSourcePt_IntrFc_y,size(FreqStamp,2));
Ar = zeros(NumSourcePt_Trans,size(FreqStamp,2));

% Pressure and Displacement Along Line
% T for Total, PR for Pressure, D for Displacement,
% L for Line, S for Bottom Transducer, R for Top Transducer

% Pressure Due to unit force
PR_L = zeros(size(ZLine,2),size(FreqStamp,2));
D_L = zeros(size(ZLine,2),size(FreqStamp,2));

% Pressure and Displacement At Transducer Surface
PR_S = zeros(NumSourcePt_Trans,size(FreqStamp,2));
D_S = zeros(NumSourcePt_Trans,size(FreqStamp,2));
PR_R = zeros(NumSourcePt_Trans,size(FreqStamp,2));
D_R = zeros(NumSourcePt_Trans,size(FreqStamp,2));

% Pressure Due to actual force
T_PR_L = zeros(size(ZLine,2),size(FreqStamp,2));
T_D_L = zeros(size(ZLine,2),size(FreqStamp,2));

% Pressure and Displacement At Transducer Surface
T_PR_S = zeros(NumSourcePt_Trans,size(FreqStamp,2));
T_D_S = zeros(NumSourcePt_Trans,size(FreqStamp,2));
T_PR_R = zeros(NumSourcePt_Trans,size(FreqStamp,2));
T_D_R = zeros(NumSourcePt_Trans,size(FreqStamp,2));

%% Loop to store source strength for different frequencies
for freqIndex = 1: size(FreqStamp,2)
    
    
    freq=FreqStamp(freqIndex);
    if freq==0
        freq=0.01;
    end
    
    omega = 2*pi*freq;
    
    %% Calculated Material Properties of Fluid
    WaveNum_P=(omega)/WaveVel_P;   % P-Wave number
    
    %% Calculating Eigen vector(FI) and eigen value(CV) using Christoffel Solution
    
    [CV,FI,Theta,Phi,R,NumTestPt]=SphereChristofel(C,Solid_rho,omega,AngTestPt);
    
    %% Calculation of the DPSM Green's Function Matrix
    
    % Calling of subprogramme PressureFluidMatt to calculate the Velocity and Pressure Green's Function Matrix for Fluid Interface
    [MSS,~,MSi,~,QSS,QiS,QSi,Q11] = PressureFluidMatt(TransCoord_Cent{1}(:,:),TransCoord_Btm{1}(:,:),IntrFcCoord_Cent{1}(:,:),IntrFcCoord_Top{1}(:,:),...
        WaveNum_P,Fluid_rho,WaveVel_P,NumIntfcSourceTot,NumSourcePt_Trans,Rotation_Trans1,1); % At Interface 1
    
    [MRR,~,MRi,~,QRR,QiR,QRi,Q22] = PressureFluidMatt(TransCoord_Cent{2}(:,:),TransCoord_Top{2}(:,:),IntrFcCoord_Cent{2}(:,:),IntrFcCoord_Btm{2}(:,:),...
        WaveNum_P,Fluid_rho,WaveVel_P,NumIntfcSourceTot,NumSourcePt_Trans,Rotation_Trans2,2); % At Interface 2
    
    % S for Source , R for Received, i for interface 1 for bottom interface 2
    % for top interface
    
    
    % Calling of subprogramme DisFluidMatt to calculate the Displacement Green's Function Matrix for Fluid Interface
    [DF3_SS,DF3_1S,DF3_S1,DF3_11]=DisFluidMatt(TransCoord_Cent{1}(:,:),TransCoord_Btm{1}(:,:),IntrFcCoord_Cent{1}(:,:),IntrFcCoord_Top{1}(:,:),WaveNum_P,freq,Fluid_rho,NumIntfcSourceTot,NumSourcePt_Trans);
    [DF3_33,DF3_R3,DF3_3R,DF3_RR]=DisFluidMatt(IntrFcCoord_Cent{2}(:,:),IntrFcCoord_Btm{2}(:,:),TransCoord_Cent{2}(:,:),TransCoord_Top{2}(:,:),WaveNum_P,freq,Fluid_rho,NumSourcePt_Trans,NumIntfcSourceTot);
    
    % clear DF3_SS DF3_S1;
    % clear DF3_RR DF3_R3;
    
    %% Calculating Green's Function
    
    SrcIndex = (NumIntfcSourceTot+1)/2;                                                        % To determine the central Point Source to impose symmetry
    CoordCentSourcePt_x=IntrFcCoord_Btm{1}(1,SrcIndex);                                   % x coordinate for Central point source
    CoordCentSourcePt_y=IntrFcCoord_Btm{1}(2,SrcIndex);                                   % y coordinate for Central point source
    CoordCentSourcePt_z=IntrFcCoord_Btm{1}(3,SrcIndex);                                   % z coordinate for Central point source
    CoordCentSourcePt_Btm=[CoordCentSourcePt_x CoordCentSourcePt_y CoordCentSourcePt_z];  % Central point source coordinate Vector for bottom interface
    
    CoordCentSourcePt_x=IntrFcCoord_Top{2}(1,SrcIndex);                                   % x coordinate for Central point source
    CoordCentSourcePt_y=IntrFcCoord_Top{2}(2,SrcIndex);                                   % y coordinate for Central point source
    CoordCentSourcePt_z=IntrFcCoord_Top{2}(3,SrcIndex);                                   % z coordinate for Central point source
    CoordCentSourcePt_Top=[CoordCentSourcePt_x CoordCentSourcePt_y CoordCentSourcePt_z];  % Central point source coordinate Vector for top interface
    
    % Calling of subprogramme Aniso_green to calculate the Displacement and Stress Green's Function Matrix for Solid Interface
    %Green's Function Matrix due to bottom point source
    [~, ~, DS3_11, ~, ~, S33_11, ~, S31_11, S32_11]...
        =Aniso_green(Sw_IntrFcCoord_Cent{1}(:,:), CoordCentSourcePt_Btm, C, Theta, Phi, CV, FI, NumTestPt, NumSourcePt_IntrFc_x, NumSourcePt_IntrFc_y, Solid_rho, omega, dTheta, dPhi);
    [~, ~, DS3_21, ~, ~, S33_21, ~, S31_21, S32_21]...
        =Aniso_green(Sw_IntrFcCoord_Cent{2}(:,:), CoordCentSourcePt_Btm, C, Theta, Phi, CV, FI, NumTestPt, NumSourcePt_IntrFc_x, NumSourcePt_IntrFc_y, Solid_rho, omega, dTheta, dPhi);
    
    %Green's Function Matrix due to top point source
    
    % [DS1_12,DS2_12,DS3_12, S11_12, S22_12, S33_12, S12_12, S31_12, S32_12]...
    %     =Aniso_green(Sw_IntrFcCoord_Cent{1}(:,:), CoordCentSourcePt_Top, C, Theta, Phi, CV, FI, NumTestPt, NumSourcePt_IntrFc_x, NumSourcePt_IntrFc_y, Solid_rho, w, dTheta, dPhi);
    % [DS1_22,DS2_22,DS3_22, S11_22, S22_22, S33_22, S12_22, S31_22, S32_22]...
    %     =Aniso_green(Sw_IntrFcCoord_Cent{2}(:,:), CoordCentSourcePt_Top, C, Theta, Phi, CV, FI, NumTestPt, NumSourcePt_IntrFc_x, NumSourcePt_IntrFc_y, Solid_rho, w, dTheta, dPhi);
    
    %By the Application of Symmetry :-->
    
    
    %For Force along 1 direction
    DS3_12(:,:,1) = -rot90(DS3_21(:,:,1),2);
    S33_12(:,:,1) = rot90(S33_21(:,:,1),2);
    S31_12(:,:,1) = -rot90(S31_21(:,:,1),2);
    S32_12(:,:,1) = -rot90(S32_21(:,:,1),2);
    
    DS3_22(:,:,1) = -rot90(DS3_11(:,:,1),2);
    S33_22(:,:,1) = rot90(S33_11(:,:,1),2);
    S31_22(:,:,1) = -rot90(S31_11(:,:,1),2);
    S32_22(:,:,1) = -rot90(S32_11(:,:,1),2);
    
    %For Force along 2 direction
    DS3_12(:,:,2) = -rot90(DS3_21(:,:,2),2);
    S33_12(:,:,2) = rot90(S33_21(:,:,2),2);
    S31_12(:,:,2) = -rot90(S31_21(:,:,2),2);
    S32_12(:,:,2) = -rot90(S32_21(:,:,2),2);
    
    DS3_22(:,:,2) = -rot90(DS3_11(:,:,2),2);
    S33_22(:,:,2) = rot90(S33_11(:,:,2),2);
    S31_22(:,:,2) = -rot90(S31_11(:,:,2),2);
    S32_22(:,:,2) = -rot90(S32_11(:,:,2),2);
    
    %For Force along 3 direction
    DS3_12(:,:,3) = rot90(DS3_21(:,:,3),2);
    S33_12(:,:,3) = -rot90(S33_21(:,:,3),2);
    S31_12(:,:,3) = rot90(S31_21(:,:,3),2);
    S32_12(:,:,3) = rot90(S32_21(:,:,3),2);
    
    DS3_22(:,:,3) = rot90(DS3_11(:,:,3),2);
    S33_22(:,:,3) = -rot90(S33_11(:,:,3),2);
    S31_22(:,:,3) = rot90(S31_11(:,:,3),2);
    S32_22(:,:,3) = rot90(S32_11(:,:,3),2);
    
    % %For Force along 1 direction
    % DS3_12(:,:,1) = rot90(DS3_21(:,:,1),2);
    % S33_12(:,:,1) = -rot90(S33_21(:,:,1),2);
    % S31_12(:,:,1) = rot90(S31_21(:,:,1),2);
    % S32_12(:,:,1) = rot90(S32_21(:,:,1),2);
    %
    % DS3_22(:,:,1) = rot90(DS3_11(:,:,1),2);
    % S33_22(:,:,1) = -rot90(S33_11(:,:,1),2);
    % S31_22(:,:,1) = rot90(S31_11(:,:,1),2);
    % S32_22(:,:,1) = rot90(S32_11(:,:,1),2);
    %
    % %For Force along 2 direction
    % DS3_12(:,:,2) = rot90(DS3_21(:,:,2),2);
    % S33_12(:,:,2) = -rot90(S33_21(:,:,2),2);
    % S31_12(:,:,2) = rot90(S31_21(:,:,2),2);
    % S32_12(:,:,2) = rot90(S32_21(:,:,2),2);
    %
    % DS3_22(:,:,2) = rot90(DS3_11(:,:,2),2);
    % S33_22(:,:,2) = -rot90(S33_11(:,:,2),2);
    % S31_22(:,:,2) = rot90(S31_11(:,:,2),2);
    % S32_22(:,:,2) = rot90(S32_11(:,:,2),2);
    %
    % %For Force along 3 direction
    % DS3_12(:,:,3) = -rot90(DS3_21(:,:,3),2);
    % S33_12(:,:,3) = rot90(S33_21(:,:,3),2);
    % S31_12(:,:,3) = -rot90(S31_21(:,:,3),2);
    % S32_12(:,:,3) = -rot90(S32_21(:,:,3),2);
    %
    % DS3_22(:,:,3) = -rot90(DS3_11(:,:,3),2);
    % S33_22(:,:,3) = rot90(S33_11(:,:,3),2);
    % S31_22(:,:,3) = -rot90(S31_11(:,:,3),2);
    % S32_22(:,:,3) = -rot90(S32_11(:,:,3),2);
    
    % save check_relation_TI_5deg.mat
    %% Creating Green's Function Matrix and Finding Source Stregth
    zero1=zeros(size(MSS));
    zero2=zeros(size(DF3_3R));
    zero3=zeros(size(DF3_33));
    zero4=zeros(size(MSi));
    zero5=[zero3 zero3 zero3];
    zero6=[zero4 zero4 zero4];
    
    MAT1=[MSS MSi zero6 zero6 zero4 zero1];
    MAT2=[DF3_1S DF3_11 -DS3_11(:,:,1) -DS3_11(:,:,2) -DS3_11(:,:,3) -DS3_12(:,:,1) -DS3_12(:,:,2) -DS3_12(:,:,3) zero3 zero2];
    MAT3=[QiS Q11 S33_11(:,:,1) S33_11(:,:,2) S33_11(:,:,3) S33_12(:,:,1) S33_12(:,:,2) S33_12(:,:,3) zero3 zero2];
    MAT4=[zero2 zero3 S31_11(:,:,1) S31_11(:,:,2) S31_11(:,:,3) S31_12(:,:,1) S31_12(:,:,2) S31_12(:,:,3) zero3 zero2];
    MAT5=[zero2 zero3 S32_11(:,:,1) S32_11(:,:,2) S32_11(:,:,3) S32_12(:,:,1) S32_12(:,:,2) S32_12(:,:,3) zero3 zero2];
    
    MAT6=[zero2 zero3 DS3_21(:,:,1) DS3_21(:,:,2) DS3_21(:,:,3) DS3_22(:,:,1) DS3_22(:,:,2) DS3_22(:,:,3) -DF3_33 -DF3_3R];
    MAT7=[zero2 zero3 S33_21(:,:,1) S33_21(:,:,2) S33_21(:,:,3) S33_22(:,:,1) S33_22(:,:,2) S33_22(:,:,3) Q22 QiR];
    MAT8=[zero2 zero3 S31_21(:,:,1) S31_21(:,:,2) S31_21(:,:,3) S31_22(:,:,1) S31_22(:,:,2) S31_22(:,:,3) zero3 zero2];
    MAT9=[zero2 zero3 S32_21(:,:,1) S32_21(:,:,2) S32_21(:,:,3) S32_22(:,:,1) S32_22(:,:,2) S32_22(:,:,3) zero3 zero2];
    MAT10=[zero1 zero4 zero6 zero6 MRi MRR];
    
    
    
    MATRIX=[MAT1 ; MAT2 ; MAT3 ; MAT4 ; MAT5 ; MAT6 ; MAT7 ; MAT8 ; MAT9 ; MAT10];    % Green’s Function Matrix
    
    
    % Calculation of Boundary Velocity Vector(Vs)
    Vs=Vso*ones(1,NumSourcePt_Trans);                   % Boundary velocity vector
    Vt=Vto*ones(1,NumSourcePt_Trans);                   % Boundary velocity vector
    Vzero=zeros(1,NumIntfcSourceTot);
    
    V=[Vs Vzero Vzero Vzero Vzero Vzero Vzero Vzero Vzero Vt]';                         % Boundary condition vector
    
    clear MAT1 MAT2 MAT3 MAT4 MAT5 MAT6 MAT7 MAT8 MAT9 MAT10;
    %% Inversion of the Green's function Matrix
    %MINV=inv(MATRIX);
    AS=MATRIX\V;                                                     % Source Strengths for all the distributed point sources
    disp ('Matrix inversion is done');
    
    %% Computation of the Source Strength and their distribution among different layers accordingly
    % Required Source Strengths: As,A1,A_2,Ar
    As(:,freqIndex)=AS(1:NumSourcePt_Trans);                                                                    % Source Strengths for the distributed point sources of transducer
    A1(:,freqIndex)=AS((NumSourcePt_Trans+1):(NumSourcePt_Trans+NumIntfcSourceTot));                                 % Source Strengths for the distributed point sources at the bottom of interface
    A_1(:,freqIndex)=AS((NumSourcePt_Trans+NumIntfcSourceTot+1):(NumSourcePt_Trans+(4*NumIntfcSourceTot)));               % Source Strengths for the distributed point sources at the top of interface
    A2(:,freqIndex)=AS((NumSourcePt_Trans+(4*NumIntfcSourceTot)+1):(NumSourcePt_Trans+(7*NumIntfcSourceTot)));
    A_2(:,freqIndex)=AS((NumSourcePt_Trans+(7*NumIntfcSourceTot)+1):(NumSourcePt_Trans+(8*NumIntfcSourceTot)));
    Ar(:,freqIndex)=AS((NumSourcePt_Trans+(8*NumIntfcSourceTot)+1):((2*NumSourcePt_Trans)+(8*NumIntfcSourceTot)));
    
    clear MATRIX;
    clear V AS;
    
    disp ('Aquisition of Source Strength is done');
    
    
    %% Calculation of the wave field at different layers and over the entire domain
    if WaveFieldMode==3
        
        if PlotMode==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For X-Y plot
            cell_TrgCoord=cell(1,1);              % Target Coordinate point for plotting
            Length_IntrFc_x=Length_IntrFc_x-DistSource_x;                   % Length of the interface along X axis
            Length_IntrFc_y=Length_IntrFc_y-DistSource_y;
            TrgCoord_x = -Length_IntrFc_x/2+IntrFcShift:Length_IntrFc_x/(NumTarget_x-1):Length_IntrFc_x/2+IntrFcShift;    % Target coordinate along x-axis
            TrgCoord_y = -Length_IntrFc_y/2:Length_IntrFc_y/(NumTarget_y-1):Length_IntrFc_y/2;                            % Target coordinate along y-axis
            
            for yin=1:NumTarget_y
                for xin=1:NumTarget_x
                    indice=xin+((yin-1)*(NumTarget_x));
                    cell_TrgCoord{1}(:,indice)=[TrgCoord_x(xin) TrgCoord_y(yin) Z_PlaneCoord ];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        if PlotMode==2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For Y-Z plot
            cell_TrgCoord=cell(1,(NumFluidLay+NumSolidLay));              % Target Coordinate point for plotting
            Length_IntrFc_y=Length_IntrFc_y-DistSource_y;                 % Length of the interface along Y axis
            TrgCoord_y = -Length_IntrFc_y/2:Length_IntrFc_y/(NumTarget_y-1):Length_IntrFc_y/2;    % Target coordinate along y-axis
            
            for ifl=1:(NumFluidLay+NumSolidLay)
                TrgCoord_z = Dist_IntrFc(ifl):(Dist_IntrFc(ifl+1)-Dist_IntrFc(ifl))/(NumTarget_z-1):Dist_IntrFc(ifl+1);    % Target coordinate along z-axis
                for zin=1:NumTarget_z
                    for yin=1:NumTarget_y
                        indice=yin+((zin-1)*(NumTarget_y));
                        cell_TrgCoord{ifl}(:,indice)=[X_PlaneCoord+IntrFcShift TrgCoord_y(yin) TrgCoord_z(zin)];
                    end
                end
            end
            
            disp ('Fluid Target point is done');
            TotNumTarget=NumTarget_y*NumTarget_z;            % Total number of target points
            
            %% Calculation of the Pressure Wave Field in the Fluid
            
            % Calling of subprogramme CompPressureFluid to calculate the Pressure Matrix for Fluid Domain
            [PR1]=CompPressureFluid(As(:,freqIndex),A1(:,freqIndex),TransCoord_Btm{1}(:,:),IntrFcCoord_Top{1}(:,:),cell_TrgCoord{1}(:,:),WaveNum_P,NumTarget_y,NumTarget_z,NumSourcePt_Trans,NumIntfcSourceTot);
            [PR4]=CompPressureFluid(A_2(:,freqIndex),Ar(:,freqIndex),IntrFcCoord_Btm{2}(:,:),TransCoord_Top{2}(:,:),cell_TrgCoord{3}(:,:),WaveNum_P,NumTarget_y,NumTarget_z,NumIntfcSourceTot,NumSourcePt_Trans);
            disp ('Pressure Calculation is done');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        if PlotMode==3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For Pressure and Displacement Calculation at the transducer surfaces
            % Transducer Parameters: TransCoord_Cent,TransCoord_Top,TransCoord_Btm,NumSourcePt_Trans
            % 1 : Bottom Transducer ; 2: Top transducer
            
            % Required Source Strengths: As,A1,A_2,Ar
            
            [~,~,QLR,QL2] = PressureFluidMatLine(PE_TargetLine,TransCoord_Top{2}(:,:),IntrFcCoord_Btm{2}(:,:),...
                WaveNum_P,Fluid_rho,WaveVel_P,NumIntfcSourceTot,NumSourcePt_Trans,Rotation_Trans2,2,size(ZLine,2)); % Along Line
            [DF3_LR,DF3_L2]=DisFluidMatLine(PE_TargetLine,TransCoord_Top{2}(:,:),IntrFcCoord_Btm{2}(:,:),...
                WaveNum_P,freq,Fluid_rho,NumIntfcSourceTot,NumSourcePt_Trans,size(ZLine,2));
            
            
            %% Calculation of the Pressure and Displacement Wave Field at Transducer "S"
            PR_S(:,freqIndex) = QSS*As(:,freqIndex) + QSi*A1(:,freqIndex);
            D_S(:,freqIndex) = DF3_SS*As(:,freqIndex) + DF3_S1*A1(:,freqIndex);
            
            %% Calculation of the Pressure and Displacement Wave Field at Transducer "R"
            PR_R(:,freqIndex) = QRR*Ar(:,freqIndex) + QRi*A_2(:,freqIndex);
            D_R(:,freqIndex) = DF3_RR*Ar(:,freqIndex) + DF3_R3*A_2(:,freqIndex);
            
            %% Calculation of the Pressure and Displacement Wave Field Along Line Below Transducer "R"
            PR_L(:,freqIndex) = QLR*Ar(:,freqIndex) + QL2*A_2(:,freqIndex);
            D_L(:,freqIndex) = DF3_LR*Ar(:,freqIndex) + DF3_L2*A_2(:,freqIndex);
            
            disp ('Pressure and Displacement Calculation is done');
            
            %% With Conjugate
            %Actual Displacement and Pressure at Transducer "S" and "R"
            %             T_PR_S(:,freqIndex) = conj(PR_S(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            %             T_D_S(:,freqIndex) = conj(D_S(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            %             T_PR_R(:,freqIndex) = conj(PR_R(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            %             T_D_R(:,freqIndex) = conj(D_R(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            %
            %             %Actual Displacement and Pressure Along Line Below Transducer "R"
            %             T_PR_L(:,freqIndex) = conj(PR_L(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            %             T_D_L(:,freqIndex) = conj(D_L(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            
            %% Without Conjugate
            %Actual Displacement and Pressure at Transducer "S" and "R"
            T_PR_S(:,freqIndex) = (PR_S(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            T_D_S(:,freqIndex) = (D_S(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            T_PR_R(:,freqIndex) = (PR_R(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            T_D_R(:,freqIndex) = (D_R(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            
            %Actual Displacement and Pressure Along Line Below Transducer "R"
            T_PR_L(:,freqIndex) = (PR_L(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            T_D_L(:,freqIndex) = (D_L(:,freqIndex))*FreqAmp(freqIndex)/NumSampPt;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if PlotMode==4
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For X-Y plot
            A_1s = zeros(3,NumIntfcSourceTot);                                                                                                        % Source Strengths for the distributed point sources at the top of interface
            A2s = zeros(3,NumIntfcSourceTot);                                                                                                          % Source Strengths for the distributed point sources at the top of interface
            for dir=1:3
                for indice=1:NumIntfcSourceTot
                    index=indice+(dir-1)*NumIntfcSourceTot;
                    A_1s(dir,indice)=A_1(index);
                    A2s(dir,indice)=A2(index);
                end
            end
            
            cell_TrgCoord=cell(1,1);              % Target Coordinate point for plotting
            Length_IntrFc_xp=Length_IntrFc_x-DistSource_x;                   % Length of the interface along X axis
            Length_IntrFc_yp=Length_IntrFc_y-DistSource_y;
            TrgCoord_x = -Length_IntrFc_xp/2+IntrFcShift:Length_IntrFc_xp/(NumTarget_x-1):Length_IntrFc_xp/2+IntrFcShift;    % Target coordinate along x-axis
            TrgCoord_y = -Length_IntrFc_yp/2:Length_IntrFc_yp/(NumTarget_y-1):Length_IntrFc_yp/2;                            % Target coordinate along y-axis
            
            for yin=1:NumTarget_y
                for xin=1:NumTarget_x
                    indice=xin+((yin-1)*(NumTarget_x));
                    cell_TrgCoord{1}(:,indice)=[TrgCoord_x(xin) TrgCoord_y(yin) Z_PlaneCoord ];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Calculation for XY plot is done');
            save ForXYPLot.mat
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For Y-Z plot
            
            A_1_cell = cell(1,3);                                                                                            % Source Strengths for the distributed point sources at the top of interface
            A2_cell = cell(1,3);
            for dir=1:3
                for i=1:NumSourcePt_IntrFc_y
                    for j=1:NumSourcePt_IntrFc_x
                        index=j+(i-1)*NumSourcePt_IntrFc_x+(dir-1)*NumIntfcSourceTot;
                        A_1_cell{dir}(j,i)=A_1(index);
                        A2_cell{dir}(j,i)=A2(index);
                    end
                end
            end
            
            A_1s = zeros(3,NumIntfcSourceTot);                                                                                            % Source Strengths for the distributed point sources at the top of interface
            A2s = zeros(3,NumIntfcSourceTot);
            for dir=1:3
                for i=1:NumSourcePt_IntrFc_x
                    for j=1:NumSourcePt_IntrFc_y
                        index=j+(i-1)*NumSourcePt_IntrFc_y;
                        A_1s(dir,index) = A_1_cell{dir}(i,j);
                        A2s(dir,index) = A2_cell{dir}(i,j);
                    end
                end
            end
            
            clear A_1_cell A2_cell ;
            
            cell_TrgCoord=cell(1,(NumFluidLay+NumSolidLay));              % Target Coordinate point for plotting
            Length_IntrFc_y=Length_IntrFc_y-DistSource_y;                 % Length of the interface along Y axis
            TrgCoord_y = -Length_IntrFc_y/2:Length_IntrFc_y/(NumTarget_y-1):Length_IntrFc_y/2;    % Target coordinate along y-axis
            
            for ifl=1:(NumFluidLay+NumSolidLay)
                TrgCoord_z = Dist_IntrFc(ifl):(Dist_IntrFc(ifl+1)-Dist_IntrFc(ifl))/(NumTarget_z-1):Dist_IntrFc(ifl+1);    % Target coordinate along z-axis
                for zin=1:NumTarget_z
                    for yin=1:NumTarget_y
                        indice=yin+((zin-1)*(NumTarget_y));
                        cell_TrgCoord{ifl}(:,indice)=[X_PlaneCoord+IntrFcShift TrgCoord_y(yin) TrgCoord_z(zin)];
                    end
                end
            end
            
            disp ('Fluid Target point is done');
            TotNumTarget=NumTarget_y*NumTarget_z;            % Total number of target points
            
            %% Calculation of the Pressure Wave Field in the Fluid
            
            % Calling of subprogramme CompPressureFluid to calculate the Pressure Matrix for Fluid Domain
            [PR1]=CompPressureFluid(As,A1,TransCoord_Btm{1}(:,:),IntrFcCoord_Top{1}(:,:),cell_TrgCoord{1}(:,:),WaveNum_P,NumTarget_y,NumTarget_z,NumSourcePt_Trans,NumIntfcSourceTot);
            [PR4]=CompPressureFluid(A_2,Ar,IntrFcCoord_Btm{2}(:,:),TransCoord_Top{2}(:,:),cell_TrgCoord{3}(:,:),WaveNum_P,NumTarget_y,NumTarget_z,NumIntfcSourceTot,NumSourcePt_Trans);
            disp ('Pressure Calculation is done');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Calculation for YZ plot is done');
            save ForYZPLot.mat
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For X-Z plot
            
            A_1s = zeros(3,NumIntfcSourceTot);                                                                                            % Source Strengths for the distributed point sources at the top of interface
            A2s = zeros(3,NumIntfcSourceTot);                                                                                            % Source Strengths for the distributed point sources at the top of interface
            for dir=1:3
                for indice=1:NumIntfcSourceTot
                    index=indice+(dir-1)*NumIntfcSourceTot;
                    A_1s(dir,indice)=A_1(index);
                    A2s(dir,indice)=A2(index);
                end
            end
            
            cell_TrgCoord=cell(1,(NumFluidLay+NumSolidLay));              % Target Coordinate point for plotting
            Length_IntrFc_x=Length_IntrFc_x-DistSource_x;                   % Length of the interface along X axis
            TrgCoord_x = -Length_IntrFc_x/2+IntrFcShift:Length_IntrFc_x/(NumTarget_x-1):Length_IntrFc_x/2+IntrFcShift;    % Target coordinate along x-axis
            
            for ifl=1:(NumFluidLay+NumSolidLay)
                TrgCoord_z = Dist_IntrFc(ifl):(Dist_IntrFc(ifl+1)-Dist_IntrFc(ifl))/(NumTarget_z-1):Dist_IntrFc(ifl+1);    % Target coordinate along z-axis
                for zin=1:NumTarget_z
                    for xin=1:NumTarget_x
                        indice=xin+((zin-1)*(NumTarget_x));
                        cell_TrgCoord{ifl}(:,indice)=[TrgCoord_x(xin) Y_PlaneCoord TrgCoord_z(zin)];
                    end
                end
            end
            
            
            
            
            disp ('Fluid Target point is done');
            TotNumTarget=NumTarget_x*NumTarget_z;            % Total number of target points
            
            %% Calculation of the Pressure Wave Field in the Fluid
            
            % Calling of subprogramme CompPressureFluid to calculate the Pressure Matrix for Fluid Domain
            [PR1]=CompPressureFluid(As,A1,TransCoord_Btm{1}(:,:),IntrFcCoord_Top{1}(:,:),cell_TrgCoord{1}(:,:),WaveNum_P,NumTarget_x,NumTarget_z,NumSourcePt_Trans,NumIntfcSourceTot);
            [PR4]=CompPressureFluid(A_2,Ar,IntrFcCoord_Btm{2}(:,:),TransCoord_Top{2}(:,:),cell_TrgCoord{3}(:,:),WaveNum_P,NumTarget_x,NumTarget_z,NumIntfcSourceTot,NumSourcePt_Trans);
            disp ('Pressure Calculation is done');
            
            disp('Calculation for XZ plot is done');
            save ForXZPLot.mat
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        % elseif (WaveFieldMode==2)
        %
        % else
        % TargetPoint=[TrgtX TrgtY TrgtZ];
    end
    
end

%% Final Calculation
% Total Pressure and Displacement for Full Signal
% Tot_PR_S = [T_PR_S fliplr(T_PR_S)];
% Tot_D_S = [T_D_S fliplr(T_D_S)];
% Tot_PR_R = [T_PR_R fliplr(T_PR_R)];
% Tot_D_R = [T_D_R fliplr(T_D_R)];
% 
% Tot_PR_L = [T_PR_L fliplr(T_PR_L)];
% Tot_D_L = [T_D_L fliplr(T_D_L)];

% Total Pressure and Displacement in Time Domain
Tot_PR_S = ifft(T_PR_S');
Tot_D_S = ifft(T_D_S');
Tot_PR_R = ifft(T_PR_R');
Tot_D_R = ifft(T_D_R');

Tot_PR_L = ifft(T_PR_L');
Tot_D_L = ifft(T_D_L');
%%
Tot_PR_S = fliplr(Tot_PR_S);
Tot_D_S = fliplr(Tot_D_S);
Tot_PR_R = fliplr(Tot_PR_R);
Tot_D_R = fliplr(Tot_D_R);

Tot_PR_L = fliplr(Tot_PR_L);
Tot_D_L = fliplr(Tot_D_L);

%%
toc;
time_taken = toc/60;
%Saving the Matlab Workspace
save TI-20mm115x3y5deg1MHzTimeSamp5BothActuation.mat
%% END OF CODE
%% Ploting of results



function [IntrFcCoord_Cent,IntrFcCoord_Top,IntrFcCoord_Btm,Sw_IntrFcCoord_Cent,Sw_IntrFcCoord_Top,Sw_IntrFcCoord_Btm] = DiscretizeInterfaces(Length_IntrFc_x,Length_IntrFc_y,NumSourcePt_IntrFc_x,NumSourcePt_IntrFc_y,NumSolidFluidIntrFc,IntrFcCoord_z,IntrFcShift)


display('i am here in geominterface')
%% Define the cell variable
IntrFcCoord_Cent=cell(1,NumSolidFluidIntrFc);      % Cell contains Coordinate of the target points at the interface :  1) x axis 2) y axis 3) z axis coordinates.
IntrFcCoord_Top=cell(1,NumSolidFluidIntrFc);       % Cell contains Coordinate of the point sources below the interface:  1) x axis 2) y axis 3) z axis coordinates.
IntrFcCoord_Btm=cell(1,NumSolidFluidIntrFc);       % Cell contains Coordinate of the point sources above the interface:  1) x axis 2) y axis 3) z axis coordinates.

%Sweeping for Solid Interface Only
Sw_IntrFcCoord_Cent=cell(1,NumSolidFluidIntrFc);   % Cell contains Coordinate of the target points at the interface :  1) x axis 2) y axis 3) z axis coordinates.
Sw_IntrFcCoord_Top=cell(1,NumSolidFluidIntrFc);    % Cell contains Coordinate of the point sources below the interface:  1) x axis 2) y axis 3) z axis coordinates.
Sw_IntrFcCoord_Btm=cell(1,NumSolidFluidIntrFc);    % Cell contains Coordinate of the point sources above the interface:  1) x axis 2) y axis 3) z axis coordinates.



%%
NumSourceTot=NumSourcePt_IntrFc_x*NumSourcePt_IntrFc_y;     % total number of source points (on either side of the interface)
IntrFcArea=Length_IntrFc_x*Length_IntrFc_y;                 % area of the interface
Source_EqivR=(IntrFcArea/(NumSourceTot*2*pi))^0.5;          % equivalent radius of the sources
DistSource_x=Length_IntrFc_x/NumSourcePt_IntrFc_x;		    % distance betxeen the sources along x axis
DistSource_y=Length_IntrFc_y/NumSourcePt_IntrFc_y;			% distance betxeen the sources along y axis

%% Distribution of sources on and either side of the interface: in=interface no. index.

X = -((NumSourcePt_IntrFc_x-1)/2)*DistSource_x+IntrFcShift:DistSource_x:((NumSourcePt_IntrFc_x-1)/2)*DistSource_x+IntrFcShift;
Y = -((NumSourcePt_IntrFc_y-1)/2)*DistSource_y:DistSource_y:((NumSourcePt_IntrFc_y-1)/2)*DistSource_y;
Z = IntrFcCoord_z(:,3);

for in=1:NumSolidFluidIntrFc
    for indy=1:NumSourcePt_IntrFc_y
        for indx=1:NumSourcePt_IntrFc_x
            indice=indx+(indy-1)*NumSourcePt_IntrFc_x;
            
            IntrFcCoord_Cent{in}(:,indice) = [ X(indx) Y(indy) Z(in)];
            IntrFcCoord_Top{in}(:,indice) = [ X(indx) Y(indy) Z(in)+Source_EqivR];
            IntrFcCoord_Btm{in}(:,indice) = [ X(indx) Y(indy) Z(in)-Source_EqivR];
            
        end
    end
end


%% For sweeping purpose
NumTrgGreenExtnd_x=(2*NumSourcePt_IntrFc_x-1);       % Here we take extended number of point sources for sweeping
NumTrgGreenExtnd_y=(2*NumSourcePt_IntrFc_y-1);
% NumSourceTot=NumTrgGreenExtnd_x*NumSourcePt_IntrFc_y;      % Total number of extended number of point sources


%% Distribution of sources on and either side of the interface: in=interface no. index.

X = -((NumTrgGreenExtnd_x-1)/2)*DistSource_x+IntrFcShift:DistSource_x:((NumTrgGreenExtnd_x-1)/2)*DistSource_x+IntrFcShift;
Y = -((NumTrgGreenExtnd_y-1)/2)*DistSource_y:DistSource_y:((NumTrgGreenExtnd_y-1)/2)*DistSource_y;
Z = IntrFcCoord_z(:,3);

for in=1:NumSolidFluidIntrFc
    for indy=1:NumTrgGreenExtnd_y
        for indx=1:NumTrgGreenExtnd_x
            
            indice=indx+(indy-1)*NumTrgGreenExtnd_x;
            Sw_IntrFcCoord_Cent{in}(:,indice) = [ X(indx) Y(indy) Z(in)];
            Sw_IntrFcCoord_Top{in}(:,indice) = [ X(indx) Y(indy) Z(in)+Source_EqivR];
            Sw_IntrFcCoord_Btm{in}(:,indice) = [ X(indx) Y(indy) Z(in)-Source_EqivR];
            
            
        end
    end
end

display('i am here in geominterface and its done')
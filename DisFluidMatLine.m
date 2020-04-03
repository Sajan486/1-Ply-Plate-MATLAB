function[DF3ls,DF3li]=DisFluidMatLine(Line,TransCoord_Btm,IntrFcCoord_Top,WaveNum_P,freq,Fluid_rho,NumSourceTot,NumSourcePt_Trans,NumLinePt)

DF3ls = zeros(NumLinePt,NumSourcePt_Trans);
DF3li = zeros(NumLinePt,NumSourceTot);

imag=sqrt(-1);
cons=1/(Fluid_rho*((2*pi*freq)^2));
%*********************************************************************************
% DF3ss
%********************************************************************************
for indice=1:NumSourcePt_Trans                            %total number of sources

CoordCentSourcePt_x=TransCoord_Btm(1,indice);
CoordCentSourcePt_y=TransCoord_Btm(2,indice);
CoordCentSourcePt_z=TransCoord_Btm(3,indice);

for i=1:NumLinePt                  % total no. of target point
R=(((Line(1,i)-CoordCentSourcePt_x)^2)+((Line(2,i)-CoordCentSourcePt_y)^2)+((Line(3,i)-CoordCentSourcePt_z)^2))^0.5;

R3=((Line(3,i)-CoordCentSourcePt_z)/R);

DF3ls(i,indice)=cons*((imag*WaveNum_P*R3*(exp(imag*WaveNum_P*R))/R)-(R3*(exp(imag*WaveNum_P*R))/(R^2)));

end

end
display( 'i am here in disfluidmmata DF3ss')

%*********************************************************************************
% DF3si
%********************************************************************************

for indice=1:NumSourceTot                            %total number of sources

CoordCentSourcePt_x=IntrFcCoord_Top(1,indice);
CoordCentSourcePt_y=IntrFcCoord_Top(2,indice);
CoordCentSourcePt_z=IntrFcCoord_Top(3,indice);

for i=1:NumLinePt                  % total no. of target point
R=(((Line(1,i)-CoordCentSourcePt_x)^2)+((Line(2,i)-CoordCentSourcePt_y)^2)+((Line(3,i)-CoordCentSourcePt_z)^2))^0.5;

R3=((Line(3,i)-CoordCentSourcePt_z)/R);

DF3li(i,indice)=cons*((imag*WaveNum_P*R3*(exp(imag*WaveNum_P*R))/R)-(R3*(exp(imag*WaveNum_P*R))/(R^2)));

end

end
display( 'i am here in disfluidmmata DF3si')

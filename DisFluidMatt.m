function[DF3ss,DF3is,DF3si,DF3ii]=DisFluidMatt(TransCoord_Cent,TransCoord_Btm,IntrFcCoord_Cent,IntrFcCoord_Top,WaveNum_P,freq,Fluid_rho,NumSourceTot,NumSourcePt_Trans)

DF3ss = zeros(NumSourcePt_Trans,NumSourcePt_Trans);
DF3is = zeros(NumSourceTot,NumSourcePt_Trans);
DF3si = zeros(NumSourcePt_Trans,NumSourceTot);
DF3ii = zeros(NumSourceTot,NumSourceTot);

imag=sqrt(-1);
cons=1/(Fluid_rho*((2*pi*freq)^2));
%*********************************************************************************
% DF3ss
%********************************************************************************
for indice=1:NumSourcePt_Trans                            %total number of sources

CoordCentSourcePt_x=TransCoord_Btm(1,indice);
CoordCentSourcePt_y=TransCoord_Btm(2,indice);
CoordCentSourcePt_z=TransCoord_Btm(3,indice);

for i=1:NumSourcePt_Trans                  % total no. of target point
R=(((TransCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((TransCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((TransCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;

R3=((TransCoord_Cent(3,i)-CoordCentSourcePt_z)/R);

DF3ss(i,indice)=cons*((imag*WaveNum_P*R3*(exp(imag*WaveNum_P*R))/R)-(R3*(exp(imag*WaveNum_P*R))/(R^2)));

end

end
display( 'i am here in disfluidmmata DF3ss')
%*********************************************************************************
% DF3is
%********************************************************************************

for indice=1:NumSourcePt_Trans                            %total number of sources

CoordCentSourcePt_x=TransCoord_Btm(1,indice);
CoordCentSourcePt_y=TransCoord_Btm(2,indice);
CoordCentSourcePt_z=TransCoord_Btm(3,indice);

for i=1:NumSourceTot                  % total no. of target point
R=(((IntrFcCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((IntrFcCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((IntrFcCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;

R3=((IntrFcCoord_Cent(3,i)-CoordCentSourcePt_z)/R);

DF3is(i,indice)=cons*((imag*WaveNum_P*R3*(exp(imag*WaveNum_P*R))/R)-(R3*(exp(imag*WaveNum_P*R))/(R^2)));

end

end
display( 'i am here in disfluidmmata DF3is')
%*********************************************************************************
% DF3si
%********************************************************************************

for indice=1:NumSourceTot                            %total number of sources

CoordCentSourcePt_x=IntrFcCoord_Top(1,indice);
CoordCentSourcePt_y=IntrFcCoord_Top(2,indice);
CoordCentSourcePt_z=IntrFcCoord_Top(3,indice);

for i=1:NumSourcePt_Trans                  % total no. of target point
R=(((TransCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((TransCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((TransCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;

R3=((TransCoord_Cent(3,i)-CoordCentSourcePt_z)/R);

DF3si(i,indice)=cons*((imag*WaveNum_P*R3*(exp(imag*WaveNum_P*R))/R)-(R3*(exp(imag*WaveNum_P*R))/(R^2)));

end

end
display( 'i am here in disfluidmmata DF3si')
%*********************************************************************************
% DF3ii
%********************************************************************************



for indice=1:NumSourceTot                            %total number of sources

CoordCentSourcePt_x=IntrFcCoord_Top(1,indice);
CoordCentSourcePt_y=IntrFcCoord_Top(2,indice);
CoordCentSourcePt_z=IntrFcCoord_Top(3,indice);

for i=1:NumSourceTot                  % total no. of target point
R=(((IntrFcCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((IntrFcCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((IntrFcCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;

R3=((IntrFcCoord_Cent(3,i)-CoordCentSourcePt_z)/R);

DF3ii(i,indice)=cons*((imag*WaveNum_P*R3*(exp(imag*WaveNum_P*R))/R)-(R3*(exp(imag*WaveNum_P*R))/(R^2)));

end

end
display( 'i am here in disfluidmmata DF3ii')
%*************************************************************************************
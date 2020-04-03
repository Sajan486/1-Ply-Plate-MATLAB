function [Mls,Mli,Qls,Qli]=PressureFluidMatLine(Line,TransCoord_Source,IntfaceCoord_Source,WaveNum_P,Fluid_rho,WaveVel_P,NumSourceTot,NumSourcePt_Trans,RotationTrans,InterfcIndex,NumLinePt)

RotationTrans=abs(RotationTrans)*pi/180;

%%
%********************************************************
%Calculation of Mss
%*********************************************************

imag=sqrt(-1);

R = zeros(NumSourcePt_Trans,NumLinePt);
x = zeros(NumSourcePt_Trans,NumLinePt);
z = zeros(NumSourcePt_Trans,NumLinePt);

for indice=1:NumSourcePt_Trans
    
    CoordCentSourcePt_x=TransCoord_Source(1,indice);
    CoordCentSourcePt_y=TransCoord_Source(2,indice);
    CoordCentSourcePt_z=TransCoord_Source(3,indice);
    
    
    
    for i=1:NumLinePt                  % total no. of sources
        R(indice,i)=(((Line(1,i)-CoordCentSourcePt_x)^2)+((Line(2,i)-CoordCentSourcePt_y)^2)+((Line(3,i)-CoordCentSourcePt_z)^2))^0.5;
        x(indice,i)=((Line(1,i)-CoordCentSourcePt_x));
        z(indice,i)=((Line(3,i)-CoordCentSourcePt_z));
    end
    
end

Mls = zeros(NumLinePt,NumSourcePt_Trans);
Qls = zeros(NumLinePt,NumSourcePt_Trans);

for is=1:NumSourcePt_Trans
    for jt=1:NumLinePt
        MxSS=((x(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        MzSS=((z(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        if InterfcIndex ==1
        Mls(jt,is)=-MxSS*sin(RotationTrans)+MzSS*cos(RotationTrans);
        elseif InterfcIndex ==2 
        Mls(jt,is)=-MxSS*sin(RotationTrans)-MzSS*cos(RotationTrans); 
        end
        Qls(jt,is)=(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
        
    end
end




display( 'i am here in mmat Mls')
%%
%********************************************************
%Calculation of Msi
%*********************************************************
R = zeros(NumSourceTot,NumLinePt);
x = zeros(NumSourceTot,NumLinePt);
z = zeros(NumSourceTot,NumLinePt);

for indice=1:NumSourceTot
    
    CoordCentSourcePt_x=IntfaceCoord_Source(1,indice);
    CoordCentSourcePt_y=IntfaceCoord_Source(2,indice);
    CoordCentSourcePt_z=IntfaceCoord_Source(3,indice);
    
    for i=1:NumLinePt                  % total no. of sources
        R(indice,i)=(((Line(1,i)-CoordCentSourcePt_x)^2)+((Line(2,i)-CoordCentSourcePt_y)^2)+((Line(3,i)-CoordCentSourcePt_z)^2))^0.5;
        x(indice,i)=((Line(1,i)-CoordCentSourcePt_x));
        z(indice,i)=((Line(3,i)-CoordCentSourcePt_z));
    end
    
end

Mli = zeros(NumLinePt,NumSourceTot);
Qli = zeros(NumLinePt,NumSourceTot);

for is=1:NumSourceTot
    for jt=1:NumLinePt
        MxS1=((x(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        MzS1=((z(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        if InterfcIndex ==1
        Mli(jt,is)=-MxS1*sin(RotationTrans)+MzS1*cos(RotationTrans);
        elseif InterfcIndex ==2
        Mli(jt,is)=-MxS1*sin(RotationTrans)-MzS1*cos(RotationTrans);
        end
        
        Qli(jt,is)=(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
        
    end
end

display( 'i am here in mmat Mli')

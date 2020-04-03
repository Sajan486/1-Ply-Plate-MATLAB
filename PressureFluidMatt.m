function [Mss,Mis,Msi,Mii,Qss,Qis,Qsi,Qii]=PressureFluidMatt(TransCoord_Cent,TransCoord_Source,IntfaceCoord_Cent,IntfaceCoord_Source,WaveNum_P,Fluid_rho,WaveVel_P,NumSourceTot,NumSourcePt_Trans,RotationTrans,InterfcIndex)

RotationTrans=abs(RotationTrans)*pi/180;

%%
%********************************************************
%Calculation of Mss
%*********************************************************

imag=sqrt(-1);

R = zeros(NumSourcePt_Trans,NumSourcePt_Trans);
x = zeros(NumSourcePt_Trans,NumSourcePt_Trans);
z = zeros(NumSourcePt_Trans,NumSourcePt_Trans);

for indice=1:NumSourcePt_Trans
    
    CoordCentSourcePt_x=TransCoord_Source(1,indice);
    CoordCentSourcePt_y=TransCoord_Source(2,indice);
    CoordCentSourcePt_z=TransCoord_Source(3,indice);
    
    
    
    for i=1:NumSourcePt_Trans                  % total no. of sources
        R(indice,i)=(((TransCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((TransCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((TransCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;
        x(indice,i)=((TransCoord_Cent(1,i)-CoordCentSourcePt_x));
        z(indice,i)=((TransCoord_Cent(3,i)-CoordCentSourcePt_z));
    end
    
end

Mss = zeros(NumSourcePt_Trans,NumSourcePt_Trans);
Qss = zeros(NumSourcePt_Trans,NumSourcePt_Trans);

for is=1:NumSourcePt_Trans
    for jt=1:NumSourcePt_Trans
        MxSS=((x(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        MzSS=((z(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        if InterfcIndex ==1
        Mss(is,jt)=-MxSS*sin(RotationTrans)+MzSS*cos(RotationTrans);
        elseif InterfcIndex ==2 
        Mss(is,jt)=-MxSS*sin(RotationTrans)-MzSS*cos(RotationTrans); 
        end
        Qss(is,jt)=(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
        
    end
end




display( 'i am here in mmat Mss')
%%
%********************************************************
%Calculation of Mis
%*********************************************************
R = zeros(NumSourcePt_Trans,NumSourceTot);
x = zeros(NumSourcePt_Trans,NumSourceTot);
z = zeros(NumSourcePt_Trans,NumSourceTot);

for indice=1:NumSourcePt_Trans
    
    CoordCentSourcePt_x=TransCoord_Source(1,indice);
    CoordCentSourcePt_y=TransCoord_Source(2,indice);
    CoordCentSourcePt_z=TransCoord_Source(3,indice);
    
    for i=1:NumSourceTot                  % total no. of sources
        R(indice,i)=(((IntfaceCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((IntfaceCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((IntfaceCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;
        x(indice,i)=((IntfaceCoord_Cent(1,i)-CoordCentSourcePt_x));
        z(indice,i)=((IntfaceCoord_Cent(3,i)-CoordCentSourcePt_z));
    end
    
end

Mis = zeros(NumSourceTot,NumSourcePt_Trans);
Qis = zeros(NumSourceTot,NumSourcePt_Trans);

for is=1:NumSourcePt_Trans
    for jt=1:NumSourceTot
        
        %Mx1S=((x(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        %Mz1S=((z(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        %Mis(jt,is)=-Mx1S*sin(Rotation_Trans2)+Mz1S*cos(Rotation_Trans2);
        Qis(jt,is)=(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
    end
end

% Mis is not required as velocity at the interface is not spescribed

display( 'i am here in mmat Mis')
%%
%********************************************************
%Calculation of Msi
%*********************************************************
R = zeros(NumSourceTot,NumSourcePt_Trans);
x = zeros(NumSourceTot,NumSourcePt_Trans);
z = zeros(NumSourceTot,NumSourcePt_Trans);

for indice=1:NumSourceTot
    
    CoordCentSourcePt_x=IntfaceCoord_Source(1,indice);
    CoordCentSourcePt_y=IntfaceCoord_Source(2,indice);
    CoordCentSourcePt_z=IntfaceCoord_Source(3,indice);
    
    for i=1:NumSourcePt_Trans                  % total no. of sources
        R(indice,i)=(((TransCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((TransCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((TransCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;
        x(indice,i)=((TransCoord_Cent(1,i)-CoordCentSourcePt_x));
        z(indice,i)=((TransCoord_Cent(3,i)-CoordCentSourcePt_z));
    end
    
end

Msi = zeros(NumSourcePt_Trans,NumSourceTot);
Qsi = zeros(NumSourcePt_Trans,NumSourceTot);

for is=1:NumSourceTot
    for jt=1:NumSourcePt_Trans
        MxS1=((x(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        MzS1=((z(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        if InterfcIndex ==1
        Msi(jt,is)=-MxS1*sin(RotationTrans)+MzS1*cos(RotationTrans);
        elseif InterfcIndex ==2
        Msi(jt,is)=-MxS1*sin(RotationTrans)-MzS1*cos(RotationTrans);
        end
        
        Qsi(jt,is)=(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
        
    end
end

display( 'i am here in mmat Msi')
%%
%********************************************************
%Calculation of Mii
%*********************************************************
R = zeros(NumSourceTot,NumSourceTot);
x = zeros(NumSourceTot,NumSourceTot);
z = zeros(NumSourceTot,NumSourceTot);

for indice=1:NumSourceTot
    
    CoordCentSourcePt_x=IntfaceCoord_Source(1,indice);
    CoordCentSourcePt_y=IntfaceCoord_Source(2,indice);
    CoordCentSourcePt_z=IntfaceCoord_Source(3,indice);
    
    for i=1:NumSourceTot                  % total no. of sources
        R(indice,i)=(((IntfaceCoord_Cent(1,i)-CoordCentSourcePt_x)^2)+((IntfaceCoord_Cent(2,i)-CoordCentSourcePt_y)^2)+((IntfaceCoord_Cent(3,i)-CoordCentSourcePt_z)^2))^0.5;
        x(indice,i)=((IntfaceCoord_Cent(1,i)-CoordCentSourcePt_x));
        z(indice,i)=((IntfaceCoord_Cent(3,i)-CoordCentSourcePt_z));
    end
    
end

Mii = zeros(NumSourceTot,NumSourceTot);
Qii = zeros(NumSourceTot,NumSourceTot);

for is=1:NumSourceTot
    for jt=1:NumSourceTot
        
        %Mx11=((x(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        %Mz11=((z(is,jt)*exp(imag*WaveNum_P*R(is,jt)))*(imag*WaveNum_P-(1/R(is,jt))))/(imag*WaveNum_P*WaveVel_P*Fluid_rho*((R(is,jt))^2));
        %Mii(is,jt)=Mx11*sin(Rotation_Trans2)+Mz11*cos(Rotation_Trans2);
        Qii(is,jt)=(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
    end
end
% Mii is not required as velocity at the interface is not spescribed
clear R x z;
display( 'i am here in mmat Mii')


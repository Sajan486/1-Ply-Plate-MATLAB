function [PR_FL] = CompPressureFluid(As,A1,TransCoord_Btm,IntrFcCoord_Top,cell_coordt,WaveNum_P,NumSourcePt_Trans1,NumSourcePt_Trans,NumSourceTot)

TotNumTarget=NumSourcePt_Trans1;
imag=sqrt(-1);

R = zeros(NumSourcePt_Trans,TotNumTarget);
P = zeros(TotNumTarget,NumSourcePt_Trans);
P1 = zeros(1,TotNumTarget);
for is=1:NumSourcePt_Trans
    CoordCentSourcePt_x=TransCoord_Btm(1,is);
    CoordCentSourcePt_y=TransCoord_Btm(2,is);
    CoordCentSourcePt_z=TransCoord_Btm(3,is);
    
    for jt=1:TotNumTarget
        R(is,jt)=(((cell_coordt(1,jt)-CoordCentSourcePt_x)^2)+((cell_coordt(2,jt)-CoordCentSourcePt_y)^2)+((cell_coordt(3,jt)-CoordCentSourcePt_z)^2))^0.5;
    end
end

for jt=1:TotNumTarget
    P1(jt)=0.0;
    for is=1:NumSourcePt_Trans
        P(jt,is)=As(is)*(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
        P1(jt)=P1(jt)+P(jt,is);
    end
end

R = zeros(NumSourceTot,TotNumTarget);
P = zeros(TotNumTarget,NumSourceTot);
P2 = zeros(1,TotNumTarget);
for is=1:NumSourceTot
    CoordCentSourcePt_x=IntrFcCoord_Top(1,is);
    CoordCentSourcePt_y=IntrFcCoord_Top(2,is);
    CoordCentSourcePt_z=IntrFcCoord_Top(3,is);
    
    for jt=1:TotNumTarget
        R(is,jt)=(((cell_coordt(1,jt)-CoordCentSourcePt_x)^2)+((cell_coordt(2,jt)-CoordCentSourcePt_y)^2)+((cell_coordt(3,jt)-CoordCentSourcePt_z)^2))^0.5;
    end
end

for jt=1:TotNumTarget
    P2(jt)=0.0;
    for is=1:NumSourceTot
        P(jt,is)=A1(is)*(exp(imag*WaveNum_P*R(is,jt)))/(R(is,jt));
        P2(jt)=P2(jt)+P(jt,is);
    end
end

PR_FL = zeros(1,TotNumTarget);
for jt=1:TotNumTarget
    PR_FL(jt)=P1(jt)+P2(jt);
end

clear R P P1 P2;
%*******************************************************************************************************************************************



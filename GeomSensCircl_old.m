function [TransCoord_Cent,TransCoord_Top,TransCoord_Btm,NumSourcePt_Trans] = GeomSensCircl(InnerR_Trans,OuterR_Trans,NumSourcePt_Trans,NumTrans,TransCoord_z,Rotation_Trans1,Rotation_Trans2,Dist_Intface2)
%%

disp('i am here in geomsens')

S=pi*(OuterR_Trans^2-InnerR_Trans^2);  % total surface area of the transducer
ds=S/NumSourcePt_Trans;            % surface area of each source
Source_EqivR=(S/(NumSourcePt_Trans*2*pi))^0.5;	 % equivalent radius of the sources

if InnerR_Trans==0
    InnerR_Trans=Source_EqivR*2^0.5;
end


Mopt=(OuterR_Trans-InnerR_Trans)/(Source_EqivR*sqrt(2*pi));
M=round(Mopt);      % number of annular ring
deltaR=(OuterR_Trans-InnerR_Trans)/M;              % thickness of each annular ring

R = zeros(1,M);
Dteta = zeros(1,M);
W = zeros(1,M);
NumSourcePt_Trans = 0;
for i=1:M
    
    R(i)=InnerR_Trans+(i-1/2)*deltaR;                 % radius of each annular ring
    Dteta(i)=ds/(deltaR*R(i));
    W(i)=2*pi/Dteta(i);                     % number of sources in each annular ring
    NumSourcePt_Trans = NumSourcePt_Trans + round(W(i));                   %total number of sources
end

W=round(W);                             % round off number of number of sources
Dteta=2*pi./W;                           % new ds/R*deltaR

rcoord = zeros(1,NumSourcePt_Trans);
imcoord = zeros(1,NumSourcePt_Trans);

indice=0;
for an=1:M                              % an  index for number of annular ring
    for s=1:W(an)                       % s index for number of sources
        indice=indice+1;
        
        rcoord(indice)=R(an)*(cos((s-1/2)*Dteta(an)));
        imcoord(indice)=R(an)*(sin((s-1/2)*Dteta(an)));
    end
end

%************************************************************************************************
% coordinates of the points (center of the sources) located on and either side of the interface
%************************************************************************************************
coord = cell(1,NumTrans);
coord1 = cell(1,NumTrans);
coord2 = cell(1,NumTrans);
TransCoord_Cent = cell(1,NumTrans);
TransCoord_Top = cell(1,NumTrans);
TransCoord_Btm = cell(1,NumTrans);

for in=1:NumTrans
    for ind=1:NumSourcePt_Trans
        coord{in}(:,ind)=[rcoord(ind) imcoord(ind) TransCoord_z(in,3)];
        coord1{in}(:,ind)=[rcoord(ind) imcoord(ind) TransCoord_z(in,3)+Source_EqivR];
        coord2{in}(:,ind)=[rcoord(ind) imcoord(ind) TransCoord_z(in,3)-Source_EqivR];
        
    end
end

Rotation_Trans1=Rotation_Trans1*pi/180;
Rotation_Trans2=Rotation_Trans2*pi/180;
c1=cos(Rotation_Trans1);
c2=cos(Rotation_Trans2);
s1=sin(Rotation_Trans1);
s2=sin(Rotation_Trans2);

%*******************************************************************************************
% defining triplet source coordinate
%********************************************************************************************
% NOT REQUIRED FOR FLUID REFER TO SOLID, HOW TO DEFINE TRIPLET.
%*****************************************************************************************************
% Defining coordinate due to rotation
%*****************************************************************************************************
for indice=1:NumSourcePt_Trans
    TransCoord_Cent{1}(:,indice) = [c1*coord{1}(1,indice)-s1*coord{1}(3,indice)  coord{1}(2,indice)  s1*coord{1}(1,indice)+c1*coord{1}(3,indice)];%+OuterR_Trans*s1];
    TransCoord_Top{1}(:,indice) = [c1*coord1{1}(1,indice)-s1*coord1{1}(3,indice)-Source_EqivR*s1 coord1{1}(2,indice) s1*coord1{1}(1,indice)+c1*coord1{1}(3,indice)+Source_EqivR*c1];%+OuterR_Trans*s1];%
    TransCoord_Btm{1}(:,indice) = [c1*coord2{1}(1,indice)-s1*coord2{1}(3,indice)+Source_EqivR*s1 coord2{1}(2,indice) s1*coord2{1}(1,indice)+c1*coord2{1}(3,indice)-Source_EqivR*c1];%+OuterR_Trans*s1];
    
    TransCoord_Cent{2}(:,indice) = [c2*coord{2}(1,indice)-s2*coord{2}(3,indice)+(Dist_Intface2*s2)  coord{2}(2,indice)  s2*coord{2}(1,indice)+c2*coord{2}(3,indice)+(Dist_Intface2*(1-c2))];%+OuterR_Trans*s2];
    TransCoord_Top{2}(:,indice) = [c2*coord1{2}(1,indice)-s2*coord1{2}(3,indice)+((Dist_Intface2-Source_EqivR)*s2) coord1{2}(2,indice) s2*coord1{2}(1,indice)+c2*coord1{2}(3,indice)+(Dist_Intface2*(1-c2))+Source_EqivR*c2];%+OuterR_Trans*s2];
    TransCoord_Btm{2}(:,indice) = [c2*coord2{2}(1,indice)-s2*coord2{2}(3,indice)+((Dist_Intface2+Source_EqivR)*s2) coord2{2}(2,indice) s2*coord2{2}(1,indice)+c2*coord2{2}(3,indice)+(Dist_Intface2*(1-c2))-Source_EqivR*c2];%+OuterR_Trans*s2];
end
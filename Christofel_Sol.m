%*********************************************************************************
%  Solution of Elastic Wave Propagation in Anisotropic Solids 
%  Definition of Materials and Solution of Christoffel's Equation
%**********************************************************************************

function[cv,Fi]=Christofel_Sol(C,NumTarget_x,NumTarget_y,NumTarget_z,Solid_rho,freq)   %


%Christoffel Acoustic Tensor Components
L11=(C(1,1)*(NumTarget_x^2))+(C(6,6)*(NumTarget_y^2))+(C(5,5)*(NumTarget_z^2))+(2*C(5,6)*(NumTarget_y*NumTarget_z))+(2*C(1,5)*(NumTarget_x*NumTarget_z))+(2*C(1,6)*(NumTarget_x*NumTarget_y));
L12=(C(1,6)*(NumTarget_x^2))+(C(2,6)*(NumTarget_y^2))+(C(4,5)*(NumTarget_z^2))+((C(1,2)+C(6,6))*(NumTarget_x*NumTarget_y))+((C(1,4)+C(5,6))*(NumTarget_x*NumTarget_z))+((C(4,6)+C(2,5))*(NumTarget_y*NumTarget_z));
L13=(C(1,5)*(NumTarget_x^2))+(C(4,6)*(NumTarget_y^2))+(C(3,5)*(NumTarget_z^2))+((C(1,4)+C(5,6))*(NumTarget_x*NumTarget_y))+((C(1,3)+C(5,5))*(NumTarget_x*NumTarget_z))+((C(3,6)+C(4,5))*(NumTarget_y*NumTarget_z));
L22=(C(6,6)*(NumTarget_x^2))+(C(2,2)*(NumTarget_y^2))+(C(4,4)*(NumTarget_z^2))+(2*C(2,6)*(NumTarget_x*NumTarget_y))+(2*C(4,6)*(NumTarget_x*NumTarget_z))+(2*C(2,4)*(NumTarget_y*NumTarget_z));
L23=(C(5,6)*(NumTarget_x^2))+(C(2,4)*(NumTarget_y^2))+(C(3,4)*(NumTarget_z^2))+((C(4,6)+C(2,5))*(NumTarget_x*NumTarget_y))+((C(3,6)+C(4,5))*(NumTarget_x*NumTarget_z))+((C(2,3)+C(4,4))*(NumTarget_y*NumTarget_z));
L33=(C(5,5)*(NumTarget_x^2))+(C(4,4)*(NumTarget_y^2))+(C(3,3)*(NumTarget_z^2))+(2*C(4,5)*(NumTarget_x*NumTarget_y))+(2*C(3,5)*(NumTarget_x*NumTarget_z))+(2*C(3,4)*(NumTarget_y*NumTarget_z));
L21=L12;
L31=L13;
L32=L23;

%Put all the L values into the matrix LIJ
LIJ=[L11 L12 L13; L21 L22 L23; L31 L32 L33]; %Christoffel Acoustic Tensor

LIJ=LIJ/(Solid_rho*(freq^2));


[Fi,Kinv]=eig(LIJ);            % eigen values are (1/k^2) , and corresponding eigen vectors 

Fi=(Fi);

cv = [sqrt((freq^2)*Kinv(1,1)) sqrt((freq^2)*Kinv(2,2)) sqrt((freq^2)*Kinv(3,3))];
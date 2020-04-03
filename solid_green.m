function[u1,u2,u3,S11,S22,S33,S23,S13,S12] = solid_green(Tr, C, Theta, Phi, CV, FI, nn, Fluid_rho, w, dTheta, dPhi)

% Initializing
uu1 = zeros(3,3);
uu2 = zeros(3,3);

U_Dw1 = zeros(3,3);
U_Dw2 = zeros(3,3);
U_Dw3 = zeros(3,3);

F_bar = [1 0 0; 0 1 0; 0 0 1]; %forcing function %% Load along three directions
mag = ((Tr(1))^2+(Tr(2))^2+(Tr(3))^2)^0.5;
cons1 = ((sqrt(-1)*w)/(8*(pi^2)*Fluid_rho)) * F_bar;
cons2 = (1/(8*(pi^2)*Fluid_rho*mag)) * F_bar;
TD_cons1 = -(w^2)/(8*(pi^2)*Fluid_rho);
TD_cons2 = (8*(pi^2)*Fluid_rho*mag^3);

for mode=1:3
    for i=1:nn+1
        for j=1:nn+1
            % v=[sin(Theta(i,j))*cos(Phi(i,j)) sin(Theta(i,j))*sin(Phi(i,j)) cos(Theta(i,j))]; %Position coordinate for the sphere domain
            v=[cos(Theta(i,j))*cos(Phi(i,j)) sin(Theta(i,j))*cos(Phi(i,j)) sin(Phi(i,j))];
            Vec = FI{i,j};%FI(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3); %Christoffel Sol: Eig Vector
            Pij=Vec(:,mode)*transpose(Vec(:,mode));    %Projection Matrix
            gam = sqrt(CV{i,j}(mode));% Gamma
            IntfaceArea=Tr(1)*v(1)+Tr(2)*v(2)+Tr(3)*v(3); %dot(Tr,v);
            
            D_cons = ((1/gam^4)*Pij*exp(sqrt(-1)*w*IntfaceArea)/gam)*(dTheta*dPhi*cos(Phi(i,j)));
            
           if IntfaceArea>0                
                soln1 = (gam^-3*Pij*exp(sqrt(-1)*w*IntfaceArea/gam)*(dTheta*dPhi*cos(Phi(i,j)))); % for displacement
                soln2 = 0;
                %First Part Derivative
                DispDw1=D_cons*v(1); % disp. derivative w.R.t. x1
                DispDw2=D_cons*v(2);  % disp. derivative w.R.t. x2
                DispDw3=D_cons*v(3);  % disp. derivative w.R.t. x3
                %Second Part Derivative
                DispD_2 = 0;
            elseif IntfaceArea==0
                %First Part Displacement
                soln1 = (gam^-3*Pij*exp(sqrt(-1)*w*IntfaceArea/gam)*(dTheta*dPhi*cos(Phi(i,j))));
                %Second Part Displacement
                soln2 = (gam^-2*Pij)*dPhi;
                
                %First Part Derivative
                DispDw1=D_cons*v(1); % disp. derivative w.R.t. x1
                DispDw2=D_cons*v(2);  % disp. derivative w.R.t. x2
                DispDw3=D_cons*v(3);  % disp. derivative w.R.t. x3
                %Second Part Derivative
                DispD_2 = (1/gam^2)*Pij*dPhi;
            else
                soln1 = 0; soln2 = 0;
                DispDw1=zeros(3,3);
                DispDw2=zeros(3,3);
                DispDw3=zeros(3,3);
                DispD_2=zeros(3,3);
                
                
            end 
            
            
            % for Green's Displacement function calculation for each mode
            uu1 = uu1 + cons1 * soln1;    %First part of disp. green's function  % Please note it is nor necessary to have F_bar as multyiplying a Matrix with an Identity matrix result the same matrix. Banerjee 07/18/2016
            uu2 = uu2 + cons2 * soln2;          %Second part of disp. green's function  % Please note it is nor necessary to have F_bar as multyiplying a Matrix with an Identity matrix result the same matrix.
            
            % total Disp. derivatives due to all possible directions
            U_Dw1=U_Dw1+(DispDw1*(TD_cons1))-(DispD_2*((Tr(1))/TD_cons2));
            U_Dw2=U_Dw2+(DispDw2*(TD_cons1))-(DispD_2*((Tr(2))/TD_cons2));
            U_Dw3=U_Dw3+(DispDw3*(TD_cons1))-(DispD_2*((Tr(3))/TD_cons2));
           
        end
    end
       
end

% Total Displacement Green's Function Due to all 3 wave modes.
DGF=uu1+uu2;

%For Displacements Green's function Allocation 
u1 = DGF(1,:);
u2 = DGF(2,:);
u3 = DGF(3,:);

% For Stress Green's Function

% Strains 
eps1 = [U_Dw1(1,1); U_Dw2(2,1); U_Dw3(3,1); 0.5*(U_Dw3(2,1)+U_Dw2(3,1)); 0.5*(U_Dw1(3,1)+U_Dw3(1,1)); 0.5*(U_Dw1(2,1)+U_Dw2(1,1))]; 
eps2 = [U_Dw1(1,2); U_Dw2(2,2); U_Dw3(3,2); 0.5*(U_Dw3(2,2)+U_Dw2(3,2)); 0.5*(U_Dw1(3,2)+U_Dw3(1,2)); 0.5*(U_Dw1(2,2)+U_Dw2(1,2))];
eps3 = [U_Dw1(1,3); U_Dw2(2,3); U_Dw3(3,3); 0.5*(U_Dw3(2,3)+U_Dw2(3,3)); 0.5*(U_Dw1(3,3)+U_Dw3(1,3)); 0.5*(U_Dw1(2,3)+U_Dw2(1,3))];
%%

% Stress Using Constitutive Eqn
sig1 = C*eps1;      % FOR force along 1 - direction
sig2 = C*eps2;      % FOR force along 2 - direction
sig3 = C*eps3;      % FOR force along 3 - direction

% Stress Green's function Allocation 
S11 = [sig1(1) sig2(1) sig3(1)];
S22 = [sig1(2) sig2(2) sig3(2)];
S33 = [sig1(3) sig2(3) sig3(3)];
S12 = [sig1(6) sig2(6) sig3(6)];
S13 = [sig1(5) sig2(5) sig3(5)];
S23 = [sig1(4) sig2(4) sig3(4)];

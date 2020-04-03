function [CV,FI,Theta,Phi,R,NumTestPt]=SphereChristofel(C,Solid_rho,w,AngTestPt)

                               
NumTestPt=double(int16(360/AngTestPt));    % convert degrees into number of test points

[sph(:,:,1),sph(:,:,2),sph(:,:,3)]=sphere(NumTestPt);
[Theta,Phi,R] = cart2sph(sph(:,:,1),sph(:,:,2),sph(:,:,3));
        
CV = cell(NumTestPt+1,NumTestPt+1);       
FI = cell(NumTestPt+1,NumTestPt+1);      

%%Christoffel Solution in all direction of sphere
for i=1:NumTestPt+1                        %Goes point by point in the n+1 by n+1 matrix from left to right
    for j=1:NumTestPt+1                    %from top to bottom
        
        X = R(i,j)*cos(Theta(i,j))*cos(Phi(i,j));
        Y = R(i,j)*sin(Theta(i,j))*cos(Phi(i,j));
        Z = R(i,j)*sin(Phi(i,j));
               
        NumTestPt_x =X /(sqrt(X^2+Y^2+Z^2));
        NumTestPt_y =Y /(sqrt(X^2+Y^2+Z^2));
        NumTestPt_z =Z /(sqrt(X^2+Y^2+Z^2));
    
        [cv,Fi]=Christofel_Sol(C,NumTestPt_x ,NumTestPt_y ,NumTestPt_z ,Solid_rho,w);
        
        %Check for imaginary numbers
        for p=1:3
            if (imag(cv(1,p))~=0)
                disp('Warning imaginary number');
            end
            for q=1:3
                if imag(Fi(p,q))~=0
                    disp('Warning imaginary numbers');
                end
                    
            end
        end
        
        CV{i,j}= cv; %Storing all Christoffel Sol: Eig Value
        FI{i,j}= Fi; %Storing all Christoffel Sol: Eig Vector
        
    end
end

display('The Christoffel Solution is Done');
function[u1,u2,u3,S331,S311,S321,S111]=CompWaveFieldSolid(PlotMode,A_1s,A2s,IntrFcCoord_Btm,IntrFcCoord_Top,Sw_TrgCoord,C, Theta, Phi, CV, FI, NumTarget_x, NumTarget_y, NumTarget_z, ExtndNumTarget_x, ExtndNumTarget_y, NumSourcePt_Intface_x, NumSourcePt_Intface_y, NumTestPt, Solid_rho, w, dTheta, dPhi)
%%
if PlotMode==1
    TotNumTarget = NumTarget_x * NumTarget_y;
    % Sw_TotNumTarget = ExtndNumTarget_x* ExtndNumTarget_y;
    
    %% Calculation Starts for t1 (Point Sources at Interface 1 - Bottom)
    u1t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    u2t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    u3t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    S11t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    S22t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    S33t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    S32t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    S31t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    S12t1=zeros(ExtndNumTarget_x,ExtndNumTarget_y,3);
    
    SrcIndex = (NumSourcePt_Intface_x*NumSourcePt_Intface_y+1)/2;
    
    % Coordinate of the central point source just 'Source_EqivR' distance below the interface
    CoordCentSourcePt_x=IntrFcCoord_Btm{1}(1,SrcIndex);
    CoordCentSourcePt_y=IntrFcCoord_Btm{1}(2,SrcIndex);
    CoordCentSourcePt_z=IntrFcCoord_Btm{1}(3,SrcIndex);
    
    for h=1:ExtndNumTarget_y
        for i=1:ExtndNumTarget_x                  % total no. of target points
            
            Sw_TrgLoc = i+(h-1)*ExtndNumTarget_x;
            Tr =[Sw_TrgCoord(Sw_TrgLoc,1)-CoordCentSourcePt_x,Sw_TrgCoord(Sw_TrgLoc,2)-CoordCentSourcePt_y, Sw_TrgCoord(Sw_TrgLoc,3)-CoordCentSourcePt_z];
            
            [u1t1(i,h,:),u2t1(i,h,:),u3t1(i,h,:), S11t1(i,h,:), S22t1(i,h,:), S33t1(i,h,:),S32t1(i,h,:),S31t1(i,h,:),S12t1(i,h,:)]=...
                solid_green(Tr, C, Theta, Phi, CV, FI, NumTestPt, Solid_rho, w, dTheta, dPhi);
            
        end
        
    end
     
    S33_t1=zeros(1,TotNumTarget);
    S31_t1=zeros(1,TotNumTarget);
    S32_t1=zeros(1,TotNumTarget);
    S11_t1=zeros(1,TotNumTarget);
    u1_t1=zeros(1,TotNumTarget);
    u2_t1=zeros(1,TotNumTarget);
    u3_t1=zeros(1,TotNumTarget);
    
    for h = 1:NumTarget_y
        for i=1:NumTarget_x  %for target points in X-Y plane
            
            index3 = i+(h-1)*NumTarget_x;
                    
            for j=1:NumSourcePt_Intface_y
                for k=1:NumSourcePt_Intface_x   %for source points in X-Y plane
                    
                    index = k+(j-1)*NumSourcePt_Intface_x;
                    index1 = NumSourcePt_Intface_x+i-k;
                    index2 = NumSourcePt_Intface_y-h+j;
                    
                    %For Stress
                    S33t=(A_1s(1,index)*S33t1(index1,index2,1))+(A_1s(2,index)*S33t1(index1,index2,2))+(A_1s(3,index)*S33t1(index1,index2,3));
                    S31t=(A_1s(1,index)*S31t1(index1,index2,1))+(A_1s(2,index)*S31t1(index1,index2,2))+(A_1s(3,index)*S31t1(index1,index2,3));
                    S32t=(A_1s(1,index)*S32t1(index1,index2,1))+(A_1s(2,index)*S32t1(index1,index2,2))+(A_1s(3,index)*S32t1(index1,index2,3));
                    S11t=(A_1s(1,index)*S11t1(index1,index2,1))+(A_1s(2,index)*S11t1(index1,index2,2))+(A_1s(3,index)*S11t1(index1,index2,3));
                    
                    S33_t1(index3)=S33_t1(index3)+S33t;
                    S31_t1(index3)=S31_t1(index3)+S31t;
                    S32_t1(index3)=S32_t1(index3)+S32t;
                    S11_t1(index3)=S11_t1(index3)+S11t;
                    
                    %For Displacement
                    u1t=(A_1s(1,index)*u1t1(index1,index2,1))+(A_1s(2,index)*u1t1(index1,index2,2))+(A_1s(3,index)*u1t1(index1,index2,3));
                    u2t=(A_1s(1,index)*u2t1(index1,index2,1))+(A_1s(2,index)*u2t1(index1,index2,2))+(A_1s(3,index)*u2t1(index1,index2,3));
                    u3t=(A_1s(1,index)*u3t1(index1,index2,1))+(A_1s(2,index)*u3t1(index1,index2,2))+(A_1s(3,index)*u3t1(index1,index2,3));
                    
                    u1_t1(index3)=u1_t1(index3)+u1t;
                    u2_t1(index3)=u2_t1(index3)+u2t;
                    u3_t1(index3)=u3_t1(index3)+u3t;
                    
                end
            end
        end
    end
    
    %% Calculation Starts for t2 (Point Sources at Interface 2 - Top)
    %For Force along 1 direction
    u1t2(:,:,1)=rot90(u1t1(:,:,1),2);
    u2t2(:,:,1)=rot90(u2t1(:,:,1),2);
    u3t2(:,:,1)=-rot90(u3t1(:,:,1),2);
    S11t2(:,:,1)=rot90(S11t1(:,:,1),2);
    %S22t2(:,:,1)=rot90(S22t1(:,:,1),2);
    S33t2(:,:,1)=rot90(S33t1(:,:,1),2);
    S32t2(:,:,1)=-rot90(S32t1(:,:,1),2);
    S31t2(:,:,1)=-rot90(S31t1(:,:,1),2);
    %S12t2(:,:,1)=rot90(S12t1(:,:,1),2);
    
    %For Force along 2 direction
    u1t2(:,:,2)=rot90(u1t1(:,:,2),2);
    u2t2(:,:,2)=rot90(u2t1(:,:,2),2);
    u3t2(:,:,2)=-rot90(u3t1(:,:,2),2);
    S11t2(:,:,2)=rot90(S11t1(:,:,2),2);
    %S22t2(:,:,2)=rot90(S22t1(:,:,2),2);
    S33t2(:,:,2)=rot90(S33t1(:,:,2),2);
    S32t2(:,:,2)=-rot90(S32t1(:,:,2),2);
    S31t2(:,:,2)=-rot90(S31t1(:,:,2),2);
    %S12t2(:,:,2)=rot90(S12t1(:,:,2),2);
    
    %For Force along 3 direction
    u1t2(:,:,3)=-rot90(u1t1(:,:,3),2);
    u2t2(:,:,3)=-rot90(u2t1(:,:,3),2);
    u3t2(:,:,3)=rot90(u3t1(:,:,3),2);
    S11t2(:,:,3)=-rot90(S11t1(:,:,3),2);
    %S22t2(:,:,3)=-rot90(S22t1(:,:,3),2);
    S33t2(:,:,3)=-rot90(S33t1(:,:,3),2);
    S32t2(:,:,3)=rot90(S32t1(:,:,3),2);
    S31t2(:,:,3)=rot90(S31t1(:,:,3),2);
    %S12t2(:,:,3)=-rot90(S12t1(:,:,3),2);

%     %For Force along 1 direction
%     u1t2(:,:,1)=-rot90(u1t1(:,:,1),2);
%     u2t2(:,:,1)=-rot90(u2t1(:,:,1),2);
%     u3t2(:,:,1)=rot90(u3t1(:,:,1),2);
%     S11t2(:,:,1)=-rot90(S11t1(:,:,1),2);
%     %S22t2(:,:,1)=-rot90(S22t1(:,:,1),2);
%     S33t2(:,:,1)=-rot90(S33t1(:,:,1),2);
%     S32t2(:,:,1)=rot90(S32t1(:,:,1),2);
%     S31t2(:,:,1)=rot90(S31t1(:,:,1),2);
%     %S12t2(:,:,1)=-rot90(S12t1(:,:,1),2);
%     
%     %For Force along 2 direction
%     u1t2(:,:,2)=-rot90(u1t1(:,:,2),2);
%     u2t2(:,:,2)=-rot90(u2t1(:,:,2),2);
%     u3t2(:,:,2)=rot90(u3t1(:,:,2),2);
%     S11t2(:,:,2)=-rot90(S11t1(:,:,2),2);
%     %S22t2(:,:,2)=-rot90(S22t1(:,:,2),2);
%     S33t2(:,:,2)=-rot90(S33t1(:,:,2),2);
%     S32t2(:,:,2)=rot90(S32t1(:,:,2),2);
%     S31t2(:,:,2)=rot90(S31t1(:,:,2),2);
%     %S12t2(:,:,2)=-rot90(S12t1(:,:,2),2);
%     
%     %For Force along 3 direction
%     u1t2(:,:,3)=rot90(u1t1(:,:,3),2);
%     u2t2(:,:,3)=rot90(u2t1(:,:,3),2);
%     u3t2(:,:,3)=-rot90(u3t1(:,:,3),2);
%     S11t2(:,:,3)=rot90(S11t1(:,:,3),2);
%     %S22t2(:,:,3)=rot90(S22t1(:,:,3),2);
%     S33t2(:,:,3)=rot90(S33t1(:,:,3),2);
%     S32t2(:,:,3)=-rot90(S32t1(:,:,3),2);
%     S31t2(:,:,3)=-rot90(S31t1(:,:,3),2);
%     %S12t2(:,:,3)=rot90(S12t1(:,:,3),2);
    
    S33_t2=zeros(1,TotNumTarget);
    S31_t2=zeros(1,TotNumTarget);
    S32_t2=zeros(1,TotNumTarget);
    S11_t2=zeros(1,TotNumTarget);
    u1_t2=zeros(1,TotNumTarget);
    u2_t2=zeros(1,TotNumTarget);
    u3_t2=zeros(1,TotNumTarget);
    
    for h = 1:NumTarget_y
        for i=1:NumTarget_x  %for target points in X-Y plane
            
            index3 = i+(h-1)*NumTarget_x;
                    
            for j=1:NumSourcePt_Intface_y
                for k=1:NumSourcePt_Intface_x   %for source points in X-Y plane
                    
                    index = k+(j-1)*NumSourcePt_Intface_x;
                    index1 = NumSourcePt_Intface_x+i-k;
                    index2 = NumSourcePt_Intface_y-h+j;
                    
                    %For Stress
                    S33t=(A2s(1,index)*S33t2(index1,index2,1))+(A2s(2,index)*S33t2(index1,index2,2))+(A2s(3,index)*S33t2(index1,index2,3));
                    S31t=(A2s(1,index)*S31t2(index1,index2,1))+(A2s(2,index)*S31t2(index1,index2,2))+(A2s(3,index)*S31t2(index1,index2,3));
                    S32t=(A2s(1,index)*S32t2(index1,index2,1))+(A2s(2,index)*S32t2(index1,index2,2))+(A2s(3,index)*S32t2(index1,index2,3));
                    S11t=(A2s(1,index)*S11t2(index1,index2,1))+(A2s(2,index)*S11t2(index1,index2,2))+(A2s(3,index)*S11t2(index1,index2,3));
                    
                    S33_t2(index3)=S33_t2(index3)+S33t;
                    S31_t2(index3)=S31_t2(index3)+S31t;
                    S32_t2(index3)=S32_t2(index3)+S32t;
                    S11_t2(index3)=S11_t2(index3)+S11t;
                    
                    %For Displacement
                    u1t=(A2s(1,index)*u1t2(index1,index2,1))+(A2s(2,index)*u1t2(index1,index2,2))+(A2s(3,index)*u1t2(index1,index2,3));
                    u2t=(A2s(1,index)*u2t2(index1,index2,1))+(A2s(2,index)*u2t2(index1,index2,2))+(A2s(3,index)*u2t2(index1,index2,3));
                    u3t=(A2s(1,index)*u3t2(index1,index2,1))+(A2s(2,index)*u3t2(index1,index2,2))+(A2s(3,index)*u3t2(index1,index2,3));
                    
                    u1_t2(index3)=u1_t2(index3)+u1t;
                    u2_t2(index3)=u2_t2(index3)+u2t;
                    u3_t2(index3)=u3_t2(index3)+u3t;
                    
                end
            end
        end
    end
    %% Final Calculation for Stress and Displacement
    S331=zeros(1,TotNumTarget);
    S311=zeros(1,TotNumTarget);
    S321=zeros(1,TotNumTarget);
    S111=zeros(1,TotNumTarget);
    u1=zeros(1,TotNumTarget);
    u2=zeros(1,TotNumTarget);
    u3=zeros(1,TotNumTarget);
    
    for i=1:TotNumTarget
        
        S331(i)=S33_t1(i)+S33_t2(i);
        S311(i)=S31_t1(i)+S31_t2(i);
        S321(i)=S32_t1(i)+S32_t2(i);
        S111(i)=S11_t1(i)+S11_t2(i);
        u1(i)=u1_t1(i)+u1_t2(i);
        u2(i)=u2_t1(i)+u2_t2(i);
        u3(i)=u3_t1(i)+u3_t2(i);
        
    end
end
%%
if PlotMode ==2
    TotNumTarget = NumTarget_y * NumTarget_z;
    Sw_TotNumTarget = ExtndNumTarget_y* NumTarget_z;
    %% Calculation Starts for t1
    u1t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    u2t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    u3t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    S11t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    S22t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    S33t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    S32t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    S31t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    S12t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_x,3);
    
    for g = 1:NumSourcePt_Intface_x % For Retrieving Centre source point for each line of source points in x axis.
        
        SrcIndex = (g-1)*NumSourcePt_Intface_y + (NumSourcePt_Intface_y+1)/2;
        
        % Coordinate of the central point source just 'Source_EqivR' distance below the interface
        CoordCentSourcePt_x=IntrFcCoord_Btm{1}(1,SrcIndex);
        CoordCentSourcePt_y=IntrFcCoord_Btm{1}(2,SrcIndex);
        CoordCentSourcePt_z=IntrFcCoord_Btm{1}(3,SrcIndex);
        
       for h=1:NumTarget_z
            for i=1:ExtndNumTarget_y                  % total no. of target points
                
                Sw_TrgLoc = i+(h-1)*ExtndNumTarget_y;
                Tr =[Sw_TrgCoord(Sw_TrgLoc,1)-CoordCentSourcePt_x,Sw_TrgCoord(Sw_TrgLoc,2)-CoordCentSourcePt_y, Sw_TrgCoord(Sw_TrgLoc,3)-CoordCentSourcePt_z];
                
                [u1t1(Sw_TrgLoc,g,:),u2t1(Sw_TrgLoc,g,:),u3t1(Sw_TrgLoc,g,:), S11t1(Sw_TrgLoc,g,:), S22t1(Sw_TrgLoc,g,:),S33t1(Sw_TrgLoc,g,:),S32t1(Sw_TrgLoc,g,:)...
                    ,S31t1(Sw_TrgLoc,g,:),S12t1(Sw_TrgLoc,g,:)] = solid_green(Tr, C, Theta, Phi, CV, FI, NumTestPt, Solid_rho, w, dTheta, dPhi);
                
            end
            
        end
        
    end
    
    S33_t1=zeros(1,TotNumTarget);
    S31_t1=zeros(1,TotNumTarget);
    S32_t1=zeros(1,TotNumTarget);
    S11_t1=zeros(1,TotNumTarget);
    u1_t1=zeros(1,TotNumTarget);
    u2_t1=zeros(1,TotNumTarget);
    u3_t1=zeros(1,TotNumTarget);
    
    for h = 1:NumTarget_z
        for i=1:NumTarget_y  %for target points in Y-Z plane
            
            for j=1:NumSourcePt_Intface_x
                for k=1:NumSourcePt_Intface_y   %for source points in X-Y plane
                    
                    index = k+(j-1)*NumSourcePt_Intface_y;
                    index1 = i+(h-1)*ExtndNumTarget_y+NumSourcePt_Intface_y-k;
                    index2 = i+(h-1)*NumTarget_y;
                    
                    %For Stress
                    S33t=(A_1s(1,index)*S33t1(index1,j,1))+(A_1s(2,index)*S33t1(index1,j,2))+(A_1s(3,index)*S33t1(index1,j,3));
                    S31t=(A_1s(1,index)*S31t1(index1,j,1))+(A_1s(2,index)*S31t1(index1,j,2))+(A_1s(3,index)*S31t1(index1,j,3));
                    S32t=(A_1s(1,index)*S32t1(index1,j,1))+(A_1s(2,index)*S32t1(index1,j,2))+(A_1s(3,index)*S32t1(index1,j,3));
                    S11t=(A_1s(1,index)*S11t1(index1,j,1))+(A_1s(2,index)*S11t1(index1,j,2))+(A_1s(3,index)*S11t1(index1,j,3));
                    
                    S33_t1(index2)=S33_t1(index2)+S33t;
                    S31_t1(index2)=S31_t1(index2)+S31t;
                    S32_t1(index2)=S32_t1(index2)+S32t;
                    S11_t1(index2)=S11_t1(index2)+S11t;
                    
                    %For Displacement
                    u1t=(A_1s(1,index)*u1t1(index1,j,1))+(A_1s(2,index)*u1t1(index1,j,2))+(A_1s(3,index)*u1t1(index1,j,3));
                    u2t=(A_1s(1,index)*u2t1(index1,j,1))+(A_1s(2,index)*u2t1(index1,j,2))+(A_1s(3,index)*u2t1(index1,j,3));
                    u3t=(A_1s(1,index)*u3t1(index1,j,1))+(A_1s(2,index)*u3t1(index1,j,2))+(A_1s(3,index)*u3t1(index1,j,3));
                    
                    u1_t1(index2)=u1_t1(index2)+u1t;
                    u2_t1(index2)=u2_t1(index2)+u2t;
                    u3_t1(index2)=u3_t1(index2)+u3t;
                    
                end
            end
        end
    end
    
    %% Calculation Starts for t2
    %For Force along 1 direction
    u1t2(:,:,1)=rot90(u1t1(:,:,1),2);
    u2t2(:,:,1)=rot90(u2t1(:,:,1),2);
    u3t2(:,:,1)=-rot90(u3t1(:,:,1),2);
    S11t2(:,:,1)=rot90(S11t1(:,:,1),2);
    %S22t2(:,:,1)=rot90(S22t1(:,:,1),2);
    S33t2(:,:,1)=rot90(S33t1(:,:,1),2);
    S32t2(:,:,1)=-rot90(S32t1(:,:,1),2);
    S31t2(:,:,1)=-rot90(S31t1(:,:,1),2);
    %S12t2(:,:,1)=rot90(S12t1(:,:,1),2);
    
    %For Force along 2 direction
    u1t2(:,:,2)=rot90(u1t1(:,:,2),2);
    u2t2(:,:,2)=rot90(u2t1(:,:,2),2);
    u3t2(:,:,2)=-rot90(u3t1(:,:,2),2);
    S11t2(:,:,2)=rot90(S11t1(:,:,2),2);
    %S22t2(:,:,2)=rot90(S22t1(:,:,2),2);
    S33t2(:,:,2)=rot90(S33t1(:,:,2),2);
    S32t2(:,:,2)=-rot90(S32t1(:,:,2),2);
    S31t2(:,:,2)=-rot90(S31t1(:,:,2),2);
    %S12t2(:,:,2)=rot90(S12t1(:,:,2),2);
    
    %For Force along 3 direction
    u1t2(:,:,3)=-rot90(u1t1(:,:,3),2);
    u2t2(:,:,3)=-rot90(u2t1(:,:,3),2);
    u3t2(:,:,3)=rot90(u3t1(:,:,3),2);
    S11t2(:,:,3)=-rot90(S11t1(:,:,3),2);
    %S22t2(:,:,3)=-rot90(S22t1(:,:,3),2);
    S33t2(:,:,3)=-rot90(S33t1(:,:,3),2);
    S32t2(:,:,3)=rot90(S32t1(:,:,3),2);
    S31t2(:,:,3)=rot90(S31t1(:,:,3),2);
    %S12t2(:,:,3)=-rot90(S12t1(:,:,3),2);

    S33_t2=zeros(1,TotNumTarget);
    S31_t2=zeros(1,TotNumTarget);
    S32_t2=zeros(1,TotNumTarget);
    S11_t2=zeros(1,TotNumTarget);
    u1_t2=zeros(1,TotNumTarget);
    u2_t2=zeros(1,TotNumTarget);
    u3_t2=zeros(1,TotNumTarget);
    
    for h = 1:NumTarget_z
        for i=1:NumTarget_y  %for target points in Y-Z plane
            
            index2 = i+(h-1)*NumTarget_y;
                    
            for j=1:NumSourcePt_Intface_x
                for k=1:NumSourcePt_Intface_y   %for source points in X-Y plane
                    
                    index = k+(j-1)*NumSourcePt_Intface_y;
                    index1 = i+(h-1)*ExtndNumTarget_y+NumSourcePt_Intface_y-k;
                    
                    %For Stress
                    S33t=(A2s(1,index)*S33t2(index1,j,1))+(A2s(2,index)*S33t2(index1,j,2))+(A2s(3,index)*S33t2(index1,j,3));
                    S31t=(A2s(1,index)*S31t2(index1,j,1))+(A2s(2,index)*S31t2(index1,j,2))+(A2s(3,index)*S31t2(index1,j,3));
                    S32t=(A2s(1,index)*S32t2(index1,j,1))+(A2s(2,index)*S32t2(index1,j,2))+(A2s(3,index)*S32t2(index1,j,3));
                    S11t=(A2s(1,index)*S11t2(index1,j,1))+(A2s(2,index)*S11t2(index1,j,2))+(A2s(3,index)*S11t2(index1,j,3));
                    
                    S33_t2(index2)=S33_t2(index2)+S33t;
                    S31_t2(index2)=S31_t2(index2)+S31t;
                    S32_t2(index2)=S32_t2(index2)+S32t;
                    S11_t2(index2)=S11_t2(index2)+S11t;
                    
                    %For Displacement
                    u1t=(A2s(1,index)*u1t2(index1,j,1))+(A2s(2,index)*u1t2(index1,j,2))+(A2s(3,index)*u1t2(index1,j,3));
                    u2t=(A2s(1,index)*u2t2(index1,j,1))+(A2s(2,index)*u2t2(index1,j,2))+(A2s(3,index)*u2t2(index1,j,3));
                    u3t=(A2s(1,index)*u3t2(index1,j,1))+(A2s(2,index)*u3t2(index1,j,2))+(A2s(3,index)*u3t2(index1,j,3));
                    
                    u1_t2(index2)=u1_t2(index2)+u1t;
                    u2_t2(index2)=u2_t2(index2)+u2t;
                    u3_t2(index2)=u3_t2(index2)+u3t;
                    
                end
            end
        end
    end
    %% Final Calculation for Stress and Displacement
    S331=zeros(1,TotNumTarget);
    S311=zeros(1,TotNumTarget);
    S321=zeros(1,TotNumTarget);
    S111=zeros(1,TotNumTarget);
    u1=zeros(1,TotNumTarget);
    u2=zeros(1,TotNumTarget);
    u3=zeros(1,TotNumTarget);
    
    for i=1:TotNumTarget
        
        S331(i)=S33_t1(i)+S33_t2(i);
        S311(i)=S31_t1(i)+S31_t2(i);
        S321(i)=S32_t1(i)+S32_t2(i);
        S111(i)=S11_t1(i)+S11_t2(i);
        u1(i)=u1_t1(i)+u1_t2(i);
        u2(i)=u2_t1(i)+u2_t2(i);
        u3(i)=u3_t1(i)+u3_t2(i);
        
    end
end
%%
if PlotMode ==3
    TotNumTarget = NumTarget_x * NumTarget_z;
    Sw_TotNumTarget = ExtndNumTarget_x* NumTarget_z;
    %% Calculation Starts for t1
    u1t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    u2t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    u3t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    S11t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    S22t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    S33t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    S32t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    S31t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    S12t1=zeros(Sw_TotNumTarget,NumSourcePt_Intface_y,3);
    
    for g = 1:NumSourcePt_Intface_y % For Retrieving Centre source point for each line of source points in y axis.
        
        SrcIndex = (g-1)*NumSourcePt_Intface_x + (NumSourcePt_Intface_x+1)/2;
        
        % Coordinate of the central point source just 'Source_EqivR' distance below the interface
        CoordCentSourcePt_x=IntrFcCoord_Btm{1}(1,SrcIndex);
        CoordCentSourcePt_y=IntrFcCoord_Btm{1}(2,SrcIndex);
        CoordCentSourcePt_z=IntrFcCoord_Btm{1}(3,SrcIndex);
        
        for h=1:NumTarget_z
            for i=1:ExtndNumTarget_x                  % total no. of target points
                
                Sw_TrgLoc = i+(h-1)*ExtndNumTarget_x;
                Tr =[Sw_TrgCoord(Sw_TrgLoc,1)-CoordCentSourcePt_x,Sw_TrgCoord(Sw_TrgLoc,2)-CoordCentSourcePt_y, Sw_TrgCoord(Sw_TrgLoc,3)-CoordCentSourcePt_z];
                
                [u1t1(Sw_TrgLoc,g,:),u2t1(Sw_TrgLoc,g,:),u3t1(Sw_TrgLoc,g,:), S11t1(Sw_TrgLoc,g,:), S22t1(Sw_TrgLoc,g,:), ...
                    S33t1(Sw_TrgLoc,g,:),S32t1(Sw_TrgLoc,g,:),S31t1(Sw_TrgLoc,g,:),S12t1(Sw_TrgLoc,g,:)]=...
                    solid_green(Tr, C, Theta, Phi, CV, FI, NumTestPt, Solid_rho, w, dTheta, dPhi);
                
            end
            
        end
        
    end
    
    S33_t1=zeros(1,TotNumTarget);
    S31_t1=zeros(1,TotNumTarget);
    S32_t1=zeros(1,TotNumTarget);
    S11_t1=zeros(1,TotNumTarget);
    u1_t1=zeros(1,TotNumTarget);
    u2_t1=zeros(1,TotNumTarget);
    u3_t1=zeros(1,TotNumTarget);
    
    for h = 1:NumTarget_z
        for i=1:NumTarget_x  %for target points in X-Z plane
            
            index2 = i+(h-1)*NumTarget_x;
                    
            for j=1:NumSourcePt_Intface_y
                for k=1:NumSourcePt_Intface_x   %for source points in X-Y plane
                    
                    index = k+(j-1)*NumSourcePt_Intface_x;
                    index1 = i+(h-1)*ExtndNumTarget_x+NumSourcePt_Intface_x-k;
                    
                    %For Stress
                    S33t=(A_1s(1,index)*S33t1(index1,j,1))+(A_1s(2,index)*S33t1(index1,j,2))+(A_1s(3,index)*S33t1(index1,j,3));
                    S31t=(A_1s(1,index)*S31t1(index1,j,1))+(A_1s(2,index)*S31t1(index1,j,2))+(A_1s(3,index)*S31t1(index1,j,3));
                    S32t=(A_1s(1,index)*S32t1(index1,j,1))+(A_1s(2,index)*S32t1(index1,j,2))+(A_1s(3,index)*S32t1(index1,j,3));
                    S11t=(A_1s(1,index)*S11t1(index1,j,1))+(A_1s(2,index)*S11t1(index1,j,2))+(A_1s(3,index)*S11t1(index1,j,3));
                    
                    S33_t1(index2)=S33_t1(index2)+S33t;
                    S31_t1(index2)=S31_t1(index2)+S31t;
                    S32_t1(index2)=S32_t1(index2)+S32t;
                    S11_t1(index2)=S11_t1(index2)+S11t;
                    
                    %For Displacement
                    u1t=(A_1s(1,index)*u1t1(index1,j,1))+(A_1s(2,index)*u1t1(index1,j,2))+(A_1s(3,index)*u1t1(index1,j,3));
                    u2t=(A_1s(1,index)*u2t1(index1,j,1))+(A_1s(2,index)*u2t1(index1,j,2))+(A_1s(3,index)*u2t1(index1,j,3));
                    u3t=(A_1s(1,index)*u3t1(index1,j,1))+(A_1s(2,index)*u3t1(index1,j,2))+(A_1s(3,index)*u3t1(index1,j,3));
                    
                    u1_t1(index2)=u1_t1(index2)+u1t;
                    u2_t1(index2)=u2_t1(index2)+u2t;
                    u3_t1(index2)=u3_t1(index2)+u3t;
                    
                end
            end
        end
    end
    %% Calculation Starts for t2
    %For Force along 1 direction
    u1t2(:,:,1)=rot90(u1t1(:,:,1),2);
    u2t2(:,:,1)=rot90(u2t1(:,:,1),2);
    u3t2(:,:,1)=-rot90(u3t1(:,:,1),2);
    S11t2(:,:,1)=rot90(S11t1(:,:,1),2);
    %S22t2(:,:,1)=rot90(S22t1(:,:,1),2);
    S33t2(:,:,1)=rot90(S33t1(:,:,1),2);
    S32t2(:,:,1)=-rot90(S32t1(:,:,1),2);
    S31t2(:,:,1)=-rot90(S31t1(:,:,1),2);
    %S12t2(:,:,1)=rot90(S12t1(:,:,1),2);
    
    %For Force along 2 direction
    u1t2(:,:,2)=rot90(u1t1(:,:,2),2);
    u2t2(:,:,2)=rot90(u2t1(:,:,2),2);
    u3t2(:,:,2)=-rot90(u3t1(:,:,2),2);
    S11t2(:,:,2)=rot90(S11t1(:,:,2),2);
    %S22t2(:,:,2)=rot90(S22t1(:,:,2),2);
    S33t2(:,:,2)=rot90(S33t1(:,:,2),2);
    S32t2(:,:,2)=-rot90(S32t1(:,:,2),2);
    S31t2(:,:,2)=-rot90(S31t1(:,:,2),2);
    %S12t2(:,:,2)=rot90(S12t1(:,:,2),2);
    
    %For Force along 3 direction
    u1t2(:,:,3)=-rot90(u1t1(:,:,3),2);
    u2t2(:,:,3)=-rot90(u2t1(:,:,3),2);
    u3t2(:,:,3)=rot90(u3t1(:,:,3),2);
    S11t2(:,:,3)=-rot90(S11t1(:,:,3),2);
    %S22t2(:,:,3)=-rot90(S22t1(:,:,3),2);
    S33t2(:,:,3)=-rot90(S33t1(:,:,3),2);
    S32t2(:,:,3)=rot90(S32t1(:,:,3),2);
    S31t2(:,:,3)=rot90(S31t1(:,:,3),2);
    %S12t2(:,:,3)=-rot90(S12t1(:,:,3),2);

    S33_t2=zeros(1,TotNumTarget);
    S31_t2=zeros(1,TotNumTarget);
    S32_t2=zeros(1,TotNumTarget);
    S11_t2=zeros(1,TotNumTarget);
    u1_t2=zeros(1,TotNumTarget);
    u2_t2=zeros(1,TotNumTarget);
    u3_t2=zeros(1,TotNumTarget);
    
    for h = 1:NumTarget_z
        for i=1:NumTarget_x  %for target points in X-Z plane
            
            index2 = i+(h-1)*NumTarget_x;
                    
            for j=1:NumSourcePt_Intface_y
                for k=1:NumSourcePt_Intface_x   %for source points in X-Y plane
                    
                    index = k+(j-1)*NumSourcePt_Intface_x;
                    index1 = i+(h-1)*ExtndNumTarget_x+NumSourcePt_Intface_x-k;
                    
                    %For Stress
                    S33t=(A2s(1,index)*S33t2(index1,j,1))+(A2s(2,index)*S33t2(index1,j,2))+(A2s(3,index)*S33t2(index1,j,3));
                    S31t=(A2s(1,index)*S31t2(index1,j,1))+(A2s(2,index)*S31t2(index1,j,2))+(A2s(3,index)*S31t2(index1,j,3));
                    S32t=(A2s(1,index)*S32t2(index1,j,1))+(A2s(2,index)*S32t2(index1,j,2))+(A2s(3,index)*S32t2(index1,j,3));
                    S11t=(A2s(1,index)*S11t2(index1,j,1))+(A2s(2,index)*S11t2(index1,j,2))+(A2s(3,index)*S11t2(index1,j,3));
                    
                    S33_t2(index2)=S33_t2(index2)+S33t;
                    S31_t2(index2)=S31_t2(index2)+S31t;
                    S32_t2(index2)=S32_t2(index2)+S32t;
                    S11_t2(index2)=S11_t2(index2)+S11t;
                    
                    %For Displacement
                    u1t=(A2s(1,index)*u1t2(index1,j,1))+(A2s(2,index)*u1t2(index1,j,2))+(A2s(3,index)*u1t2(index1,j,3));
                    u2t=(A2s(1,index)*u2t2(index1,j,1))+(A2s(2,index)*u2t2(index1,j,2))+(A2s(3,index)*u2t2(index1,j,3));
                    u3t=(A2s(1,index)*u3t2(index1,j,1))+(A2s(2,index)*u3t2(index1,j,2))+(A2s(3,index)*u3t2(index1,j,3));
                    
                    u1_t2(index2)=u1_t2(index2)+u1t;
                    u2_t2(index2)=u2_t2(index2)+u2t;
                    u3_t2(index2)=u3_t2(index2)+u3t;
                    
                end
            end
        end
    end
    %% Final Calculation for Stress and Displacement
    S331=zeros(1,TotNumTarget);
    S311=zeros(1,TotNumTarget);
    S321=zeros(1,TotNumTarget);
    S111=zeros(1,TotNumTarget);
    u1=zeros(1,TotNumTarget);
    u2=zeros(1,TotNumTarget);
    u3=zeros(1,TotNumTarget);
    
    for i=1:TotNumTarget
        
        S331(i)=S33_t1(i)+S33_t2(i);
        S311(i)=S31_t1(i)+S31_t2(i);
        S321(i)=S32_t1(i)+S32_t2(i);
        S111(i)=S11_t1(i)+S11_t2(i);
        u1(i)=u1_t1(i)+u1_t2(i);
        u2(i)=u2_t1(i)+u2_t2(i);
        u3(i)=u3_t1(i)+u3_t2(i);
        
    end
    
end

%*******************************************************************************************************************

%******************************************************************************************************************
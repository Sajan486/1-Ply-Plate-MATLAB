% function[u1,u2,u3]=solid_disp(cons,WaveNum_P,WaveNum_S,InnerR_Trans,OuterR_Trans,R3,R,lam,miu);
function[DS1ii,DS2ii,DS3ii, S11ii, S22ii, S33ii, S12ii, S31ii, S32ii]=Aniso_green(Sw_IntrFcCoord_Cent, CoordCentSourcePt, C, Theta, Phi, CV, FI, NumTestPt, NumSourcePt_Intface_x, NumSourcePt_Intface_y,  Solid_rho, w, dTheta, dPhi)

%tcoordi = total target point
%%Displacement Green's Function

%% Initialization
NumSourceTot=NumSourcePt_Intface_x*NumSourcePt_Intface_y;

% For Final allocation of Displacements and Stresses
DS1ii = zeros(NumSourceTot,NumSourceTot,3);
DS2ii = zeros(NumSourceTot,NumSourceTot,3);
DS3ii = zeros(NumSourceTot,NumSourceTot,3);

S11ii = zeros(NumSourceTot,NumSourceTot,3);
S22ii = zeros(NumSourceTot,NumSourceTot,3);
S33ii = zeros(NumSourceTot,NumSourceTot,3);
S12ii = zeros(NumSourceTot,NumSourceTot,3);
S31ii = zeros(NumSourceTot,NumSourceTot,3);
S32ii = zeros(NumSourceTot,NumSourceTot,3);


% ntt = (2*NumSourcePt_Intface_x-1)*NumSourcePt_Intface_y;
NumTrgGreenExtnd_x=(2*NumSourcePt_Intface_x-1);
NumTrgGreenExtnd_y=(2*NumSourcePt_Intface_y-1);

% For Temporary allocation of Displacements and Stresses
tu1_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tu2_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tu3_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);

tS11_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tS22_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tS33_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tS12_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tS31_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);
tS32_11 = zeros(NumTrgGreenExtnd_x,NumTrgGreenExtnd_y,3);

%% Calculation Starts

CoordCentSourcePt_x=CoordCentSourcePt(1);
CoordCentSourcePt_y=CoordCentSourcePt(2);
CoordCentSourcePt_z=CoordCentSourcePt(3);

for h=1:NumTrgGreenExtnd_y
    for i=1:NumTrgGreenExtnd_x                  % total no. of target point
        
        Sw_TrgLoc = i+(h-1)*NumTrgGreenExtnd_x;
        Tr = [(Sw_IntrFcCoord_Cent(1,Sw_TrgLoc)-CoordCentSourcePt_x), (Sw_IntrFcCoord_Cent(2,Sw_TrgLoc)-CoordCentSourcePt_y), (Sw_IntrFcCoord_Cent(3,Sw_TrgLoc)-CoordCentSourcePt_z)];
        
        [u1,u2,u3,S11,S22,S33,S23,S13,S12] = solid_green(Tr, C, Theta, Phi, CV, FI, NumTestPt,  Solid_rho, w, dTheta, dPhi);
        
        tu1_11(i,h,:) = u1;
        tu2_11(i,h,:) = u2;
        tu3_11(i,h,:) = u3;
        tS11_11(i,h,:) = S11;
        tS22_11(i,h,:) = S22;
        tS33_11(i,h,:) = S33;
        tS12_11(i,h,:) = S12;
        tS31_11(i,h,:) = S13;
        tS32_11(i,h,:) = S23;
        
    end
end

for indicey=1:NumSourcePt_Intface_y
    for indicex=1:NumSourcePt_Intface_x                %total number of sources
        
        src_pt = indicex+(indicey-1)*NumSourcePt_Intface_x;
        trg_pts = (NumSourcePt_Intface_x+1-indicex):(NumTrgGreenExtnd_x+1-indicex);     
                
        for dir=1:3                  % total no. of target points
            for indiceyp = 1:NumSourcePt_Intface_y
                
                sw_index_y = NumSourcePt_Intface_y-indicey+indiceyp;
                trg_val = (indiceyp-1)*NumSourcePt_Intface_x+1:(indiceyp)*NumSourcePt_Intface_x;
                
                DS1ii(trg_val,src_pt,dir) = tu1_11(trg_pts,sw_index_y,dir);    % Displacement in mm
                DS2ii(trg_val,src_pt,dir) = tu2_11(trg_pts,sw_index_y,dir);    % Displacement in mm
                DS3ii(trg_val,src_pt,dir) = tu3_11(trg_pts,sw_index_y,dir);    % Displacement in mm
                
                S11ii(trg_val,src_pt,dir) = tS11_11(trg_pts,sw_index_y,dir);  % Stress in GPa
                S22ii(trg_val,src_pt,dir) = tS22_11(trg_pts,sw_index_y,dir);  % Stress in GPa
                S33ii(trg_val,src_pt,dir) = tS33_11(trg_pts,sw_index_y,dir);  % Stress in GPa
                S12ii(trg_val,src_pt,dir) = tS12_11(trg_pts,sw_index_y,dir);  % Stress in GPa
                S31ii(trg_val,src_pt,dir) = tS31_11(trg_pts,sw_index_y,dir);  % Stress in GPa
                S32ii(trg_val,src_pt,dir) = tS32_11(trg_pts,sw_index_y,dir);  % Stress in GPa
                
                
            end
        end
    end
end
disp('Green displacement and stress calculated')
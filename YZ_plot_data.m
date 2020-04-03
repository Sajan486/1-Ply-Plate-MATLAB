%*****************************************************************************************************************
% Ploting of results
%*****************************************************************************************************************
%%
figure
surfc(yp(:,52:102),zp(:,52:102),abs(S33(:,52:102)),'edgecolor','none'), rotate3d on;
xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Stress33 distribution in Isotropic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',18);
axis equal;
colorbar;
view([0 90]);
print('S_S33','-djpeg')
%%
figure
surfc(yp(:,1:51),zp(:,1:51),abs(S33(:,1:51)),'edgecolor','none'), rotate3d on;
xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Pressure distribution in Isotropic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',18);
axis equal;
colorbar;
view([0 90]);
print('S_p','-djpeg')
%%

% figure(2)
% contour(yp,zp,abs(S33),100)
% xlabel('Y axis in mm','FontSize',15,'FontWeight','bold')
% ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Contour plot of Stress33 distribution','FontSize',20)
% set(gca,'FontWeight','bold','FontSize',18);
% print('C_S33','-djpeg')

% surfc(yp(:,1:n2),zp(,abs(P(:,1:n2)), rotate3d on;
% xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
% ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Stress31 distribution','FontSize',20)

% figure(4)
% contour(yp,zp,abs(P(:,1:n2)),60)
% xlabel('Y axis in mm','FontSize',15,'FontWeight','bold')
% ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Contour plot of Stress31 distribution','FontSize',20)
%%
figure
surfc(yp(:,52:102),zp(:,52:102),abs(S11(:,52:102)),'edgecolor','none'), rotate3d on;
xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Stress11 distribution in Isotropic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',18);
axis equal;
colorbar;
view([0 90]);
print('S_S11','-djpeg')


%%
% figure(6)
% contour(yp,zp,abs(S11),100)
% xlabel('Y axis in mm','FontSize',15,'FontWeight','bold')
% ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Contour plot of Stress11 distribution','FontSize',20)
% set(gca,'FontWeight','bold','FontSize',18);
% print('C_S11','-djpeg')

% figure(7)
% surfc(yp(:,1:51),zp(:,1:51),abs(S11(:,1:51))), rotate3d on;
% xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
% ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
% title('pressure distribution','FontSize',20)
% figure(8)
% contour(yp(:,1:51),zp(:,1:51),abs(S11(:,1:51)),200)
% xlabel('Y axis in mm','FontSize',15,'FontWeight','bold')
% ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Contour plot of pressure distribution','FontSize',20)

figure
surfc(yp(:,52:102),zp(:,52:102),abs(u_3(:,52:102)),'edgecolor','none'), rotate3d on;
xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   u3 distribution in Isotropic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',18);
axis equal;
colorbar;
view([0 90]);
print('S_u3','-djpeg')

% figure(8)
% contour(yp,zp,abs(u_3),60)
% xlabel('Y axis in mm','FontSize',15,'FontWeight','bold')
% ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Contour plot of u3 distribution','FontSize',20)
% set(gca,'FontWeight','bold','FontSize',18);
% print('C_u3','-djpeg')
%%
figure
surfc(yp(:,52:102),zp(:,52:102),abs(u_1(:,52:102)),'edgecolor','none'), rotate3d on;
xlabel('Y axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   u1 distribution in Isotropic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',18);
axis equal;
colorbar;
view([0 90]);
print('S_u1','-djpeg')

% figure(10)
% contour(yp,zp,abs(u_1),60)
% xlabel('Y axis in mm','FontSize',15,'FontWeight','bold')
% ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% title('Contour plot of u1 distribution','FontSize',20)
% print('C_u1','-djpeg')

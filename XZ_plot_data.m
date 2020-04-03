%*****************************************************************************************************************
% Ploting of results
%*****************************************************************************************************************
%%
figure
surfc(xp(:,52:102),zp(:,52:102),abs(S33(:,52:102)),'edgecolor','none');% , rotate3d on;
xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Stress33 distribution in Monoclinic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',20);
axis equal;
colormap jet;
colorbar;
view([0 90]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 28 14])
print('S_S33','-djpeg')
%%
figure
surfc(xp(:,1:51),zp(:,1:51),abs(S33(:,1:51)),'edgecolor','none');% , rotate3d on;
xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Pressure (P1) distribution in Monoclinic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',20);
axis equal;
colormap jet;
colorbar;
view([0 90]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 28 14])
print('S_p1','-djpeg')
%%
figure
surfc(xp(:,103:153),zp(:,103:153),abs(S33(:,103:153)),'edgecolor','none');% , rotate3d on;
xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Pressure (P4) distribution in Monoclinic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',20);
axis equal;
colormap jet;
colorbar;
view([0 90]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 28 14])
print('S_p4','-djpeg')
%%

% figure(2)
% contour(xp,zp,abs(S33),100)
% %xlabel('X axis in mm','FontSize',15,'FontWeight','bold')
% %ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Contour plot of Stress33 distribution','FontSize',20)
% set(gca,'FontWeight','bold','FontSize',20);
% print('C_S33','-djpeg')

% surfc(xp(:,1:n2),zp(,abs(P(:,1:n2));% , rotate3d on;
% %xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
% %ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Stress31 distribution','FontSize',20)

% figure(4)
% contour(xp,zp,abs(P(:,1:n2)),60)
% %xlabel('X axis in mm','FontSize',15,'FontWeight','bold')
% %ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Contour plot of Stress31 distribution','FontSize',20)
%%
figure
surfc(xp(:,52:102),zp(:,52:102),abs(S11(:,52:102)),'edgecolor','none');% , rotate3d on;
xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   Stress11 distribution in Monoclinic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',20);
axis equal;
colormap jet;
colorbar;
view([0 90]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 28 14])
print('S_S11','-djpeg')


%%
% figure(6)
% contour(xp,zp,abs(S11),100)
% %xlabel('X axis in mm','FontSize',15,'FontWeight','bold')
% %ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Contour plot of Stress11 distribution','FontSize',20)
% set(gca,'FontWeight','bold','FontSize',20);
% print('C_S11','-djpeg')

% figure(7)
% surfc(xp(:,1:51),zp(:,1:51),abs(S11(:,1:51)));% , rotate3d on;
% %xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
% %ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('pressure distribution','FontSize',20)
% figure(8)
% contour(xp(:,1:51),zp(:,1:51),abs(S11(:,1:51)),200)
% %xlabel('X axis in mm','FontSize',15,'FontWeight','bold')
% %ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Contour plot of pressure distribution','FontSize',20)

figure
surfc(xp(:,52:102),zp(:,52:102),abs(u_3(:,52:102)),'edgecolor','none');% , rotate3d on;
xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   u3 distribution in Monoclinic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',20);
axis equal;
colormap jet;
colorbar;
view([0 90]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 28 14])
print('S_u3','-djpeg')

% figure(8)
% contour(xp,zp,abs(u_3),60)
% %xlabel('X axis in mm','FontSize',15,'FontWeight','bold')
% %ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Contour plot of u3 distribution','FontSize',20)
% set(gca,'FontWeight','bold','FontSize',20);
% print('C_u3','-djpeg')
%%
figure
% figure('units','normalized','position',[.15 .15 .45 .45])
surfc(xp(:,52:102),zp(:,52:102),abs(u_1(:,52:102)),'edgecolor','none');% , rotate3d on;
xlabel('X axis  in mm','FontSize',15,'FontWeight','bold')
ylabel(' Z axis in mm','FontSize',15,'FontWeight','bold')
title({'   u1 distribution in Monoclinic Halfspace',' '},'FontSize',19,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',20);
axis equal;
colormap jet;
colorbar;
view([0 90]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 28 14])
print('S_u1','-djpeg')

% figure(10)
% contour(xp,zp,abs(u_1),60)
% %xlabel('X axis in mm','FontSize',15,'FontWeight','bold')
% %ylabel('Z axis in mm','FontSize',15,'FontWeight','bold')
% %title('Contour plot of u1 distribution','FontSize',20)
% print('C_u1','-djpeg')

% close all
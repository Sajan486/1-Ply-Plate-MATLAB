figure
plot(IntrFcCoord_Cent{1}(1,:),abs(S33(:,51)),'Linewidth',3)
hold on
plot(IntrFcCoord_Cent{1}(1,:),abs(S33(:,52)),'ro','Linewidth',2)
hold off
xlabel('X axis  in mm','FontSize',18,'FontWeight','bold')
ylabel(' Stress/Pressure ','FontSize',18,'FontWeight','bold')
title('Stress-Pressure at Fluid 2-Solid Boundary','FontSize',25)
set(gca,'FontWeight','bold','FontSize',18);
% axis equal
legend('Pressure','Stress')
print('match2','-djpeg')
%%
figure
plot(IntrFcCoord_Cent{1}(1,:),abs(S33(:,102)),'Linewidth',3)
hold on
plot(IntrFcCoord_Cent{1}(1,:),abs(S33(:,103)),'ro','Linewidth',2)
hold off
xlabel('X axis  in mm','FontSize',18,'FontWeight','bold')
ylabel(' Stress/Pressure ','FontSize',18,'FontWeight','bold')
title('Stress-Pressure at Fluid 1-Solid Boundary ','FontSize',25)
set(gca,'FontWeight','bold','FontSize',18);
% axis equal
legend('Pressure','Stress')
print('match1','-djpeg')
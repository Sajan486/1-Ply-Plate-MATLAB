plot(TransCoord_Cent{2}(1,:),TransCoord_Cent{2}(3,:))
hold on
 plot(TransCoord_Cent{1}(1,:),TransCoord_Cent{1}(3,:))
axis equal
plot(IntrFcCoord_Cent{1}(1,:),IntrFcCoord_Cent{1}(3,:))
plot(IntrFcCoord_Cent{2}(1,:),IntrFcCoord_Cent{2}(3,:))

plot(TransCoord_Top{2}(1,:),TransCoord_Top{2}(3,:),'ko')
plot(TransCoord_Btm{2}(1,:),TransCoord_Btm{2}(3,:),'ro')

 plot(TransCoord_Top{1}(1,:),TransCoord_Top{1}(3,:),'ko')
 plot(TransCoord_Btm{1}(1,:),TransCoord_Btm{1}(3,:),'ro')

plot(IntrFcCoord_Top{2}(1,:),IntrFcCoord_Top{2}(3,:),'ko')
plot(IntrFcCoord_Btm{2}(1,:),IntrFcCoord_Btm{2}(3,:),'ro')

plot(IntrFcCoord_Top{1}(1,:),IntrFcCoord_Top{1}(3,:),'ko')
plot(IntrFcCoord_Btm{1}(1,:),IntrFcCoord_Btm{1}(3,:),'ro')
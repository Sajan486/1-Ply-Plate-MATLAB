% Plot Comparison
% plot(timeSignal,FtSignal);
% figure;plot(timeSignal,(Tot_D_R(1,:))');

%% Plot at Transducer
% figure;plot(timeSignal,real(mean(Tot_D_R)/(4e-5)), 'r');
% figure;plot(timeSignal,real(mean(Tot_PR_R)/(2.5e-4)), 'r');
% hold on
% plot(TimeStamp,ForceTimeSignal);
% ylim([-2,2])

% yy1 = smooth(TimeStamp,mean(Tot_PR_R,2),0.000001,'rloess');
% plot(TimeStamp,real(yy1), 'r');
% mean(Tot_PR_R,2)

plot(TimeStamp,ForceTimeSignal);
hold on
plot(TimeStamp(1:128),real(mean(Tot_PR_R,2)/4e-3), 'r');
% plot(TimeStamp,real(mean(Tot_D_R,2)/(1e-5)), 'r');
% plot(TimeStamp(1:128),real(Tot_D_L(1:128,41)/(4e-4)), 'r');
ylim([-2 2])
hold off
%% Plot Loop
% For Displacement along line
for i = 1:41
    plot(TimeStamp(1:128),real(Tot_D_L(:,i)), 'r');
    %ylim([-1 1])
    ylim([-0.002 0.002])
    %ylim([-0.0005 0.0005])
<<<<<<< HEAD
    pause(0.4);
end
%%
% For Pressure along line
% for i = 1:41
%     plot(TimeStamp,real(Tot_PR_L(:,i)), 'r');
%     %ylim([-1 1])
%     ylim([-0.002 0.002])
%     %ylim([-0.0005 0.0005])
%     pause(0.4);
% end
=======
    pause(0.3);
end
%%
% For Pressure along line
for i = 1:41
    plot(TimeStamp(1:128),real(Tot_PR_L(:,i)), 'r');
    %ylim([-0.01 0.01])
    ylim([-0.007 0.007])
    %ylim([-0.0005 0.0005])
    pause(0.3);
end
>>>>>>> 39ccb1c748f7a2dacdd05811f0632dbb36e96b6f


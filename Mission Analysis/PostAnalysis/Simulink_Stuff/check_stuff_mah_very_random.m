for i=2:length(MA_time_correct)
    aaaaaa(i,:) = TrajSC(i,:)-TrajSC(i-1,:);
end
figure()
plot(aaaaaa(:,1))
legend('x')
figure()
plot(aaaaaa(:,2))
legend('y')
figure()
plot(aaaaaa(:,3))
legend('z')
figure()
plot(aaaaaa(:,4))
legend('t')
figure()
comet3(TrajSC(:,1),TrajSC(:,2),TrajSC(:,3))
hold on
comet3(MA.SC1.uniform.R(:,1),MA.SC1.uniform.R(:,2),MA.SC1.uniform.R(:,3))
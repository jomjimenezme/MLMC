close all;
figure 
rectangle('Position',[579.9354, 0, 2*199.1502,1400],'FaceColor',[1, 0, 0, 0.2])
hold on
%axis([650 1000 0 1400])
plot(779.086,0,".r", 'MarkerSize',15)
xline(779.086,".r")

rectangle('Position',[584.1003, 0, 2*199.1502,1400],'FaceColor',[0, 0, 1, 0.2])
plot(783.251 ,0,".b", 'MarkerSize',15)
xline(783.251, ".b")

title("TOL 1E-1")
saveas(gcf,'1E1.png')
hold off



%%
close all
figure
rectangle('Position',[759.1706, 0, 2*19.91502,1400],'FaceColor',[1, 0, 0, 0.2])
hold on
%axis([650 1000 0 1400])
plot(779.086,0,".r", 'MarkerSize',15)
xline(779.086,".r")

rectangle('Position',[763.2647, 0, 2*19.91502,1400],'FaceColor',[0, 0, 1, 0.2])
plot(783.251 ,0,".b", 'MarkerSize',15)
xline(783.251, ".b")

title("TOL 1E-2")
saveas(gcf,'1E2.png')
hold off



%%
%1E-3
close all
figure
rectangle('Position',[775.9851, 0, 2*1.991502,1400],'FaceColor',[1, 0, 0, 0.2])
hold on
%axis([650 1000 0 1400])
plot(777.977 ,0,".r", 'MarkerSize',15)
xline(777.977,".r")

rectangle('Position',[775.5248, 0, 2*1.991502,1400],'FaceColor',[0, 0, 1, 0.2])
plot(777.516 ,0,".b", 'MarkerSize',15)
xline(777.516, ".b")

title("TOL 1E-3")
saveas(gcf,'1E3.png')
hold off



%%
%1E-4
close all
figure
rectangle('Position',[778.2109, 0, 2*0.1991502,1400],'FaceColor',[1, 0, 0, 0.2])
hold on
%axis([650 1000 0 1400])
plot(778.410,0,".r", 'MarkerSize',15)
xline(778.410,".r")

rectangle('Position',[778.1140, 0, 2*0.1991502,1400],'FaceColor',[0, 0, 1, 0.2])
plot(778.313 ,0,".b", 'MarkerSize',15)
xline(778.313, ".b")
title("TOL 1E-4")
saveas(gcf,'1E4.png')
hold off




%%
close all
figure
rectangle('Position',[759.1706, 0, 2*19.91502,1400],'FaceColor',[1, 0, 0, 0.2])
hold on
%axis([650 1000 0 1400])



rectangle('Position',[763.2647, 0, 2*19.91502,1400],'FaceColor',[0, 0, 1, 0.2])





rectangle('Position',[775.9851, 0, 2*1.991502,1400],'FaceColor',[1, 0, 0, 0.2])

plot(777.977 ,0,".r", 'MarkerSize',15)
xline(777.977,".r")

rectangle('Position',[775.5248, 0, 2*1.991502,1400],'FaceColor',[0, 0, 1, 0.2])
plot(777.516 ,0,".b", 'MarkerSize',15)
xline(777.516, ".b")

title("TOL 1E-2 and 1E-3")


saveas(gcf,'Zcompare_2-3.png')
hold off





%%



close all;
%1E-3
figure
rectangle('Position',[775.9851, 0, 2*1.991502,1400],'FaceColor',[1, 0, 0, 0.1])
hold on
%axis([650 1000 0 1400])


rectangle('Position',[775.5248, 0, 2*1.991502,1400],'FaceColor',[0, 0, 1, 0.1])







%1E-4
rectangle('Position',[778.2109, 0, 2*0.1991502,1400],'FaceColor',[1, 0, 0, 0.5])

%axis([650 1000 0 1400])
plot(778.410,0,".r", 'MarkerSize',15)
xline(778.410,".r")

rectangle('Position',[778.1140, 0, 2*0.1991502,1400],'FaceColor',[0, 0, 1, 0.5])
plot(778.313 ,0,".b", 'MarkerSize',15)
xline(778.313, ".b")


title("TOL 1E-3 and 1E-4")
saveas(gcf,'Zcompare_3-4.png')
hold off





%%







%1E-4
close all
load("plain_values.mat")
figure
%hist(a,'FaceColor',[1, 0, 0])
hold on
rectangle('Position',[778.2109, 0, 2*0.1991502,10],'FaceColor',[1, 0, 0, 0.2])

%axis([650 1000 0 1400])
plot(778.410,0,".r", 'MarkerSize',15)
xline(778.410,".r")
plot(a,1, "sr")


load("mlmc_values.mat")
rectangle('Position',[778.1140, 0, 2*0.1991502,10],'FaceColor',[0, 0, 1, 0.2])
plot(778.313 ,0,".b", 'MarkerSize',15)
xline(778.313, ".b")
plot(b,2, "sb")

title("TOL 1E-4")
saveas(gcf,'ZSeveral1E-4.png')
hold off




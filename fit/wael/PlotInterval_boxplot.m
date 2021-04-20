function []= PlotInterval_boxplot(xAxisData,MeanData, StandDev, NumOfDev,color)

%% Description:
%% Input: 
%xAxisData: the x-axis data
 %MeanData: the vector containing the means 
 %StandDev: the vector containing the standard dev for each mean
 %NumOfDev: How many stds is the confidence interval, defult is 1
 
 %% Output: an addition to the plot containing the specified interval
 
 %% Function starts here
 %if last input is missing
 if nargin==3
     NumOfDev=1;
     color='r';
 end
 %% plots the interval
 hold on
 ebs=0.06;
 for i=1:length(xAxisData)
     up_lim=MeanData(i)+StandDev(i)*NumOfDev;
     down_lim=MeanData(i)-StandDev(i)*NumOfDev;
     
     plot([xAxisData(i) xAxisData(i)],[down_lim up_lim],color)
     plot([xAxisData(i)-ebs xAxisData(i)+ebs],[up_lim up_lim],color)
     plot([xAxisData(i)-ebs xAxisData(i)+ebs],[down_lim down_lim],color)
 end
 plot(xAxisData,MeanData,'s','MarkerSize',5,'MarkerFaceColor',color,'MarkerEdgeColor',color,'LineWidth',1.2)
 plot(xAxisData,MeanData,color,'LineWidth',1.2)
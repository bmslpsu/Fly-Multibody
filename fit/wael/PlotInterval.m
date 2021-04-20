function []= PlotInterval(xAxisData,MeanData, StandDev, NumOfDev,color)

%% Description: Plots the interval with a shaded region to indicate the standard deviation
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
 plot(xAxisData,MeanData,color,'LineWidth',2)
 curve1=MeanData+StandDev*NumOfDev;
 curve2=MeanData-StandDev*NumOfDev;
 x2=[xAxisData, fliplr(xAxisData)];
 inBetween= [curve1, fliplr(curve2)];
 patch(x2, inBetween, color,'FaceAlpha',.2,'EdgeColor',color);
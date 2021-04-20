function [p_values]=GenerateBodePlotByArea(root_area,NumberOfBins,U_D,Frf_D)
figure
AreaBins=linspace(0.6,0.9,NumberOfBins+1);
[fileName,~] = uigetfile(root_area,'MultiSelect', 'on','*.mat');
fly_array=U_D.fly{1};
for i=1:NumberOfBins
    Gain=[];
    Phase=[];
    for k=1:length(fileName)
        
        for j=1:length(fly_array)
            temp=split(fileName{k},[".","_"]);
            fly_num=str2double(temp{2});
            if fly_num==fly_array(j)
                load([root_area fileName{k}]);
                if save_area>=AreaBins(i) && save_area<AreaBins(i+1)
                    FlyData(i).WingAreaRange=AreaBins(i:i+1);
                    Gain=[Frf_D.ref2body.fly.gain(:,j) Gain];
                    Phase=[Frf_D.ref2body.fly.phase(:,j) Phase];
                end
            end
        end
    end
    colors=['b','g','r','k','c','m'];
    FlyData(i).Gain=Gain;
    FlyData(i).Phase=Phase;
    FlyData(i).AreaRange=AreaBins(i:i+1);
    subplot(2,1,1)
    PlotInterval_boxplot(Frf_D.IOFv{1},mean(Gain,2),std(Gain')',1,colors(i))
    subplot(2,1,2)
    PlotInterval_boxplot(Frf_D.IOFv{1},mean(Phase,2),std(Phase')',1,colors(i))
    hold on
end

for i=1:9
    [~,p_val1(i)]=ttest2(FlyData(1).Gain(i,:),FlyData(2).Gain(i,:));
end
for i=1:9
    [~,p_val2(i)]=ttest2(FlyData(1).Gain(i,:),FlyData(3).Gain(i,:));
end
for i=1:9
    [~,p_val3(i)]=ttest2(FlyData(2).Gain(i,:),FlyData(3).Gain(i,:));
end
p_values=[p_val1;p_val2;p_val3];
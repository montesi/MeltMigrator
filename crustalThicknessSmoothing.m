function FiltedCrust=crustalThicknessSmoothing(Geometry,Crust,SmoothingWidth,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crustalThicknessSmoothing.m
% Smooth and plot crustal thickness profile along plate boundary segments 
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

figure(indFigure); clf; hold on;

% % generate plot frames
yMax=[20,30];
Label(1).String='Thickness [km]';
Label(2).String='Average and Max. F [%]';
for iSubplot=1:numel(Label)
    subplot(numel(Label),1,iSubplot); hold on;
    for iSegment=1:numel(Geometry.PlateBoundaryLength);
        plot([1,1]*sum(Geometry.PlateBoundaryLength(1:iSegment)),[0,1]*yMax(iSubplot),'k','linewidth',1);
        ylim([0,yMax(iSubplot)]);
        xlim([0,sum(Geometry.PlateBoundaryLength)]);
        title(Label(iSubplot).String);
        box on;
    end
end

if SmoothingWidth==0;
    subplot(numel(Label),1,1)
    plot(Crust.AlongPlateBoundaryCoordinate,abs(Crust.Thickness),'color','r','linewidth',1);
    subplot(numel(Label),1,2)
    plot(Crust.AlongPlateBoundaryCoordinate,Crust.AverageMeltFraction*100,'color','r','linewidth',1);
    subplot(numel(Label),1,2)
    plot(Crust.AlongPlateBoundaryCoordinate,Crust.MaxMeltFraction*100,'color','b','linewidth',1);
else
    EffectiveSmoothingWidth=max(max(find((Crust.AlongPlateBoundaryCoordinate-Crust.AlongPlateBoundaryCoordinate(1))<SmoothingWidth)),1);
    Thickness=[];
    AverageMeltFraction=[];
    MaxMeltFraction=[];
    MeanThickness=[];
    
    % % smooth crustal thickness for each segment
    iStart=1;
    for iSegment=1:numel(Geometry.PlateBoundaryLength);
        [minDis,iEnd]=min(abs(Crust.AlongPlateBoundaryCoordinate-sum(Geometry.PlateBoundaryLength(1:iSegment))));
        FiltedThickness=fastSmooth(Crust.Thickness(iStart:iEnd),EffectiveSmoothingWidth,3,1);
        FiltedAverageMeltFraction=fastSmooth(Crust.AverageMeltFraction(iStart:iEnd),EffectiveSmoothingWidth,3,1);
        FiltedMaxMeltFraction=fastSmooth(Crust.MaxMeltFraction(iStart:iEnd),EffectiveSmoothingWidth,3,1);
        MeanThickness=[MeanThickness,zeros(size(FiltedThickness))+mean(Crust.Thickness(iStart:iEnd))];
        Thickness=[Thickness,FiltedThickness];
        AverageMeltFraction=[AverageMeltFraction,FiltedAverageMeltFraction];
        MaxMeltFraction=[MaxMeltFraction,FiltedMaxMeltFraction];
        iStart=1+iEnd;
    end
    
    subplot(numel(Label),1,1);
    plot(Crust.AlongPlateBoundaryCoordinate,Thickness,'color','r','linewidth',1.5);
    plot(Crust.AlongPlateBoundaryCoordinate,MeanThickness,'color','b','linewidth',1);
    subplot(numel(Label),1,2)
    plot(Crust.AlongPlateBoundaryCoordinate,AverageMeltFraction*100,'color','r','linewidth',1);
    subplot(numel(Label),1,2)
    plot(Crust.AlongPlateBoundaryCoordinate,MaxMeltFraction*100,'color','b','linewidth',1);
end

FiltedCrust.Thickness=Thickness;
FiltedCrust.MeanThickness=MeanThickness;
FiltedCrust.AverageMeltFraction=AverageMeltFraction;
FiltedCrust.MaxMeltFraction=MaxMeltFraction;

return

function plotCrustalAccretionProfile(Geometry,Shoulder,Crust,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotCrustalAccretionProfile.m
% Plot crustal accretion profile (including crustal thickness, maximum & average melt fraction and melting depth 
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% June 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% % initialize figure
figure(indFigure); clf; hold on;
Color=[1,0,0;...
    0,1,0;...
    0.5,0.5,0.5;... 
    1,1,1;...0,0,0;...    
    1,1,0;...
    1,0,1;...
    0,0,1];
ColorSize=size(Color,1);

% % generate plot frames
yMax=[20,30,100];
Label(1).String='Thickness [km]';
Label(2).String='Average and Maximum Melt Fraction [%]';
Label(3).String='Average Melting Depth [km]';
for iSubplot=1:3
    subplot(3,1,iSubplot); hold on;
    for iSegment=1:numel(Geometry.PlateBoundaryLength);
        plot([1 ,1]*sum(Geometry.PlateBoundaryLength(1:iSegment)),[0,1]*yMax(iSubplot),'k','linewidth',1);
        title(Label(iSubplot).String);
        ylim([0,yMax(iSubplot)]);
        xlim([0,sum(Geometry.PlateBoundaryLength)]);
        box on;
    end
end

% % plot contribution from each shoulder
for iShoulder=1:numel(Shoulder);
    if numel(Shoulder(iShoulder).Distribution.Thickness)>0;
        subplot(3,1,1);
        plot(Shoulder(iShoulder).Distribution.Thickness,'color',Color(mod(iShoulder-1,ColorSize)+1,:));
        subplot(3,1,2);
        plot(Shoulder(iShoulder).Distribution.AverageMeltFraction*100,'color',Color(mod(iShoulder-1,ColorSize)+1,:));
        subplot(3,1,2);
        plot(Shoulder(iShoulder).Distribution.MaxMeltFraction*100,'color',Color(mod(iShoulder-1,ColorSize)+1,:));
        subplot(3,1,3);
        plot(Shoulder(iShoulder).Distribution.MeltingDepth,'color',Color(mod(iShoulder-1,ColorSize)+1,:));
    end
end

% % plot integrated values
subplot(3,1,1);
plot(Crust.Thickness,'color','k','linewidth',2);
subplot(3,1,2);
plot(Crust.AverageMeltFraction*100,'color','k','linewidth',2);
subplot(3,1,2);
plot(Crust.MaxMeltFraction*100,'color','k','linewidth',2);
subplot(3,1,3);
plot(Crust.MeltingDepth,'color','k','linewidth',2);

clear
close all

addpath ./cbrewer
load quantiles.mat % quantiles of fake data estimated using gamlss_train_brain.R
age=data.age;
y=data.phenotype;
sex=data.sex;

% color
cb=cbrewer('qual', 'Set1', 8, 'pchip');
cb_red=cb(1,:);
cb_blue=cb(2,:);
cb2=cbrewer('seq', 'Greys', 256, 'pchip');
cb2=cb2(30:end,:);

hf=figure; hf.Color='w';hf.Position=[100 100 500 250];
for i=1:2
    subplot(1,2,i)

    hs=dscatter(age(sex==i-1),y(sex==i-1));
    colormap(cb2);

    hold on
    if i==1
        hp=plot(xF,yF,'Color',cb_red,'LineWidth',1,'LineStyle','--');
        hp(3).LineStyle='-';
        hp(3).LineWidth=2;

    elseif i==2
        hp=plot(xM,yM,'Color',cb_blue,'LineWidth',1,'LineStyle','--');
        hp(3).LineStyle='-';
        hp(3).LineWidth=2;
    end
    ha=gca;
    ha.FontSize=14;
    ha.XLim=[18,100];
    ha.YLim=[1.5,3];
    ha.XTick=[20:10:100];
    ha.YLabel.String='Phenotype'; 
    ha.XLabel.String='Age';
    ha.YGrid='on';
    ha.TickLength=[0.02,0.02];
end
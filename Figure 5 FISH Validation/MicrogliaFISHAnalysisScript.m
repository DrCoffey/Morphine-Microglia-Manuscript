%% Script to batch process IHC FISH Overlap Images for Morphine Microglia Project 2016-2021........ 
close all
clear all

%% PDE10A Options
Options.fishThreshold = .175;
Options.fishFilter = 20;
Options.ihcThreshold = .15;
Options.ihcFilter = 10000;
Options.makeFig3D = 0;
Options.figFold = 'C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Analysed Images\PDE10A';

d=dir('C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\To Analyse\PDE10A');
d=d(3:end,:);

for i=1:height(d);
    fname=fullfile(d(i).folder,d(i).name);
    [colocTable] = microgliaColc3D(fname, Options);
    if i==1;
        mT=colocTable;
    else
        mT=[mT; colocTable];
    end
end
save('PDE10A_Output.mat','Options','mT');

%% Start Here After Image Processing
load('PDE10A_Output.mat');

% Normalized Colocalization
func = @(x,y) 100*(y/x);
B = rowfun(func,mT,'InputVariables',{'cellSize','redColoc'},...
    'OutputVariableName','normColoc');
masterT=[mT B];

% Determine % of microglia with overlap
expThreshold=450;
IDs=unique(masterT.ID);
for i=1:length(IDs)
percentExpress(i,1)=sum(masterT.redColoc>expThreshold & masterT.ID==IDs(i))/sum(masterT.ID==IDs(i))*100;
cellsAnalysed(i,1)=sum(masterT.ID==IDs(i));
end

% Group Stats
avgT=grpstats(masterT,'ID','mean','DataVars',{'cellSize','redColoc','normColoc'});
key=readtable("C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Key.xlsx");
key.Group=categorical(key.Group);
key.Sex=categorical(key.Sex);
key.ID=categorical(key.ID);
avgT=join(avgT,key);
masterT=join(masterT,key);
avgT=[avgT table(percentExpress,cellsAnalysed)];

g=gramm('x',avgT.Group,'y',avgT.percentExpress,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 100]);
g.set_names('x','','y','% Microglia Expressing PDE10A','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','PDE10A Percent Express.png'),'-m3');

g=gramm('x',avgT.Group,'y',avgT.mean_normColoc,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 10]);
g.set_names('x','','y','Normalized PDE10A Expression','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','PDE10A Expression.png'),'-m3');

[p,percentExpress_Anova,stats] = anova1([avgT.percentExpress(avgT.Group=='SS') avgT.percentExpress(avgT.Group=='MS') avgT.percentExpress(avgT.Group=='MN')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,normColoc_Anova,stats] = anova1([avgT.mean_normColoc(avgT.Group=='SS') avgT.mean_normColoc(avgT.Group=='MS') avgT.mean_normColoc(avgT.Group=='MN')]);
figure;
[normColoc_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

save('PDE10A_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');

%% ARPP21
% Script to batch process IHC FISH Overlap Images for Morphine Microglia Project 2016-2021........ 
close all
clear all

%% ARPP21 Options
Options.fishThreshold = .175;
Options.fishFilter = 20;
Options.ihcThreshold = .15;
Options.ihcFilter = 10000;
Options.makeFig3D = 0;
Options.figFold = 'C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Analysed Images\ARPP21';

d=dir('C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\To Analyse\ARPP21');
d=d(3:end,:);

for i=1:height(d);
    fname=fullfile(d(i).folder,d(i).name);
    [colocTable] = microgliaColc3D(fname, Options);
    if i==1;
        mT=colocTable;
    else
        mT=[mT; colocTable];
    end
end
save('ARPP21_Output.mat','Options','mT');

%% Start Here After Image Processing
load('ARPP21_Output.mat');

% Normalized Colocalization
func = @(x,y) 100*(y/x);
B = rowfun(func,mT,'InputVariables',{'cellSize','redColoc'},...
    'OutputVariableName','normColoc');
masterT=[mT B];

% Determine % of microglia with overlap
expThreshold=450;
IDs=unique(masterT.ID);
for i=1:length(IDs)
percentExpress(i,1)=sum(masterT.redColoc>expThreshold & masterT.ID==IDs(i))/sum(masterT.ID==IDs(i))*100;
cellsAnalysed(i,1)=sum(masterT.ID==IDs(i));
end

% Group Stats
avgT=grpstats(masterT,'ID','mean','DataVars',{'cellSize','redColoc','normColoc'});
key=readtable("C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Key.xlsx");
key.Group=categorical(key.Group);
key.Sex=categorical(key.Sex);
key.ID=categorical(key.ID);
avgT=join(avgT,key);
masterT=join(masterT,key);
avgT=[avgT table(percentExpress,cellsAnalysed)];

g=gramm('x',avgT.Group,'y',avgT.percentExpress,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 100]);
g.set_names('x','','y','% Microglia Expressing ARPP21','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','ARPP21 Percent Express.png'),'-m3');

g=gramm('x',avgT.Group,'y',avgT.mean_normColoc,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 10]);
g.set_names('x','','y','Normalized ARPP21 Expression','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','ARPP21 Expression.png'),'-m3');

[p,percentExpress_Anova,stats] = anova1([avgT.percentExpress(avgT.Group=='SS') avgT.percentExpress(avgT.Group=='MS') avgT.percentExpress(avgT.Group=='MN')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,normColoc_Anova,stats] = anova1([avgT.mean_normColoc(avgT.Group=='SS') avgT.mean_normColoc(avgT.Group=='MS') avgT.mean_normColoc(avgT.Group=='MN')]);
figure;
[normColoc_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

save('ARPP21_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');

%% ITPKA
% Script to batch process IHC FISH Overlap Images for Morphine Microglia Project 2016-2021........ 
close all
clear all

%% ITPKA Options
Options.fishThreshold = .175;
Options.fishFilter = 20;
Options.ihcThreshold = .15;
Options.ihcFilter = 10000;
Options.makeFig3D = 0;
Options.figFold = 'C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Analysed Images\ITPKA';

d=dir('C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\To Analyse\ITPKA');
d=d(3:end,:);

for i=1:height(d);
    fname=fullfile(d(i).folder,d(i).name);
    [colocTable] = microgliaColc3D(fname, Options);
    if i==1;
        mT=colocTable;
    else
        mT=[mT; colocTable];
    end
end
save('ITPKA_Output.mat','Options','mT');
%% Start Here After Image Processing
load('ITPKA_Output.mat');

% Normalized Colocalization
func = @(x,y) 100*(y/x);
B = rowfun(func,mT,'InputVariables',{'cellSize','redColoc'},...
    'OutputVariableName','normColoc');
masterT=[mT B];

% Determine % of microglia with overlap
expThreshold=450;
IDs=unique(masterT.ID);
for i=1:length(IDs)
percentExpress(i,1)=sum(masterT.redColoc>expThreshold & masterT.ID==IDs(i))/sum(masterT.ID==IDs(i))*100;
cellsAnalysed(i,1)=sum(masterT.ID==IDs(i));
end

% Group Stats
avgT=grpstats(masterT,'ID','mean','DataVars',{'cellSize','redColoc','normColoc'});
key=readtable("C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Key.xlsx");
key.Group=categorical(key.Group);
key.Sex=categorical(key.Sex);
key.ID=categorical(key.ID);
avgT=join(avgT,key);
masterT=join(masterT,key);
avgT=[avgT table(percentExpress,cellsAnalysed)];

g=gramm('x',avgT.Group,'y',avgT.percentExpress,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 100]);
g.set_names('x','','y','% Microglia Expressing ITPKA','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','ITPKA Percent Express.png'),'-m3');

g=gramm('x',avgT.Group,'y',avgT.mean_normColoc,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 10]);
g.set_names('x','','y','Normalized ITPKA Expression','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','ITPKA Expression.png'),'-m3');

[p,percentExpress_Anova,stats] = anova1([avgT.percentExpress(avgT.Group=='SS') avgT.percentExpress(avgT.Group=='MS') avgT.percentExpress(avgT.Group=='MN')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,normColoc_Anova,stats] = anova1([avgT.mean_normColoc(avgT.Group=='SS') avgT.mean_normColoc(avgT.Group=='MS') avgT.mean_normColoc(avgT.Group=='MN')]);
figure;
[normColoc_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

save('ITPKA_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');

%% DRD1 Options
Options.fishThreshold = .175;
Options.fishFilter = 20;
Options.ihcThreshold = .15;
Options.ihcFilter = 10000;
Options.makeFig3D = 0;
Options.figFold = 'C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Analysed Images\DRD1';

d=dir('C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\To Analyse\DRD1');
d=d(3:end,:);

for i=1:height(d);
    fname=fullfile(d(i).folder,d(i).name);
    [colocTable] = microgliaColc3D(fname, Options);
    if i==1;
        mT=colocTable;
    else
        mT=[mT; colocTable];
    end
end
save('DRD1_Output.mat','Options','mT');

%% Start Here After Image Processing
load('DRD1_Output.mat');

% Normalized Colocalization
func = @(x,y) 100*(y/x);
B = rowfun(func,mT,'InputVariables',{'cellSize','redColoc'},...
    'OutputVariableName','normColoc');
masterT=[mT B];

% Determine % of microglia with overlap
expThreshold=450;
IDs=unique(masterT.ID);
for i=1:length(IDs)
percentExpress(i,1)=sum(masterT.redColoc>expThreshold & masterT.ID==IDs(i))/sum(masterT.ID==IDs(i))*100;
cellsAnalysed(i,1)=sum(masterT.ID==IDs(i));
end

% Group Stats
avgT=grpstats(masterT,'ID','mean','DataVars',{'cellSize','redColoc','normColoc'});
key=readtable("C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Key.xlsx");
key.Group=categorical(key.Group);
key.Sex=categorical(key.Sex);
key.ID=categorical(key.ID);
avgT=join(avgT,key);
masterT=join(masterT,key);
avgT=[avgT table(percentExpress,cellsAnalysed)];

g=gramm('x',avgT.Group,'y',avgT.percentExpress,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 100]);
g.set_names('x','','y','% Microglia Expressing DRD1','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','DRD1 Percent Express.png'),'-m3');

g=gramm('x',avgT.Group,'y',avgT.mean_normColoc,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 10]);
g.set_names('x','','y','Normalized DRD1 Expression','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','DRD1 Expression.png'),'-m3');

[p,percentExpress_Anova,stats] = anova1([avgT.percentExpress(avgT.Group=='SS') avgT.percentExpress(avgT.Group=='MS') avgT.percentExpress(avgT.Group=='MN')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,normColoc_Anova,stats] = anova1([avgT.mean_normColoc(avgT.Group=='SS') avgT.mean_normColoc(avgT.Group=='MS') avgT.mean_normColoc(avgT.Group=='MN')]);
figure;
[normColoc_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

save('DRD1_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');
%% Control Options
Options.fishThreshold = .175;
Options.fishFilter = 20;
Options.ihcThreshold = .15;
Options.ihcFilter = 10000;
Options.makeFig3D = 0;
Options.figFold = 'C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Analysed Images\Control';

d=dir('C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\To Analyse\Control');
d=d(3:end,:);

for i=1:height(d);
    fname=fullfile(d(i).folder,d(i).name);
    [colocTable] = microgliaColc3D(fname, Options);
    if i==1;
        mT=colocTable;
    else
        mT=[mT; colocTable];
    end
end
save('Control_Output.mat','Options','mT');

%% Start Here After Image Processing
load('Control_Output.mat');

% Normalized Colocalization
func = @(x,y) 100*(y/x);
B = rowfun(func,mT,'InputVariables',{'cellSize','redColoc'},...
    'OutputVariableName','normColoc');
masterT=[mT B];

% Determine % of microglia with overlap
expThreshold=450;
IDs=unique(masterT.ID);
for i=1:length(IDs)
percentExpress(i,1)=sum(masterT.redColoc>expThreshold & masterT.ID==IDs(i))/sum(masterT.ID==IDs(i))*100;
cellsAnalysed(i,1)=sum(masterT.ID==IDs(i));
end

% Group Stats
avgT=grpstats(masterT,'ID','mean','DataVars',{'cellSize','redColoc','normColoc'});
key=readtable("C:\Users\DrCoffey\Documents\Neumaier Lab\RNA Scope Experiment\IHC and FISH\Morphine Microglia FISH IHC\Key.xlsx");
key.Group=categorical(key.Group);
key.Sex=categorical(key.Sex);
key.ID=categorical(key.ID);
avgT=join(avgT,key);
masterT=join(masterT,key);
avgT=[avgT table(percentExpress,cellsAnalysed)];

g=gramm('x',avgT.Group,'y',avgT.percentExpress,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 100]);
g.set_names('x','','y','% Microglia Expressing Control','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','Control Percent Express.png'),'-m3');

g=gramm('x',avgT.Group,'y',avgT.mean_normColoc,'color',avgT.Group);
g.stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',0);
g.axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 10]);
g.set_names('x','','y','Normalized Control Expression','color','Groups');
g.set_order_options('x',{'SS','MS','MN'},'color',{'SS','MS','MN'});
g.set_layout_options('legend',0);
figure('Position',[100 100 300 300]);
g.draw();
export_fig(fullfile('.\Figures','Control Expression.png'),'-m3');

[p,percentExpress_Anova,stats] = anova1([avgT.percentExpress(avgT.Group=='SS') avgT.percentExpress(avgT.Group=='MS') avgT.percentExpress(avgT.Group=='MN')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,normColoc_Anova,stats] = anova1([avgT.mean_normColoc(avgT.Group=='SS') avgT.mean_normColoc(avgT.Group=='MS') avgT.mean_normColoc(avgT.Group=='MN')]);
figure;
[normColoc_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

save('Control_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');
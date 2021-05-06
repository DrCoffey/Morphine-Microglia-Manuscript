%% Morphine Naloxone Behavior Figures with Gramm(You must have Gramm installed!)
% Navigate to the folder containing the code!
close all 
clear all

% Determine Image Depth (Z Levels)
files=dir('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\Microglia Morphology Analysis\Analysis Figures');
files=files(3:18);
for i=1:length(files)
    info = imfinfo(fullfile(files(i).folder,files(i).name));
    Slice(i)=size(info,1)/2;
end

% Determine Image Depth (Z Levels)
files=dir('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Scope Experiment\Morphine Microglia FISH IHC\To Analyse\DRD1');
files=files(3:end);
for i=1:length(files)
    info = imfinfo(fullfile(files(i).folder,files(i).name));
    Slice(end+1)=size(info,1);
end

files=dir('C:\Users\Kevin\Documents\GitHub\Manuscripts\Morphine Microglia (unpublished)\Figure 10 Microglia Morphology\Morphology Data');
files=files(3:end);
for i=1:length(files)
t=readtable(fullfile(files(i).folder,files(i).name));
parts=strsplit(files(i).name,'_');
Subject=categorical(repmat(parts(2),[height(t) 1]));
tmp=strtok(parts(4),'.');
A1=categorical(repmat({tmp{1,1}(1)},[height(t) 1]));
A2=categorical(repmat({tmp{1,1}(2)},[height(t) 1]));
Group=categorical(repmat(strtok(parts(4),'.'),[height(t) 1]));
t=[table(Subject,Group,A1,A2) t];
if i==1
    mt=t;
else
    mt=[mt;t];
end
end

%histfit(mt.CellVolume)
mt=mt(mt.numbranchpts>0,:); % Remove Non-Traced Cells
G = groupsummary(mt,{'Subject','Group','A1','A2'},@(x) median((x)));

% BU=G;
% G=BU(17:34,:);
% BUS=Slice';
% Slice=BUS(17:34,:);
% g(1,1)=gramm('x',mt.Group,'y',mt.CellVolume,'color',mt.Group);
% g(1,2)=gramm('x',mt.Group,'y',mt.numbranchpts,'color',mt.Group);
% g(1,3)=gramm('x',mt.Group,'y',mt.MeanTerminalBranchLength,'color',mt.Group);

g(1,1)=gramm('x',G.Group,'y',G.fun1_FullCellTerritoryVol./Slice,'color',G.Group,'subset',G.Group ~='SN');
g(1,2)=gramm('x',G.Group,'y',G.fun1_numbranchpts./Slice,'color',G.Group,'subset',G.Group ~='SN');
g(1,3)=gramm('x',G.Group,'y',G.fun1_numendpts./Slice,'color',G.Group,'subset',G.Group ~='SN');

% Violin
g(1,1).stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g(1,1).geom_jitter('width',.1,'dodge',1,'alpha',.5);
g(1,1).stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
g(1,1).axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 100000]);

g(1,2).stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g(1,2).geom_jitter('width',.1,'dodge',1,'alpha',.5);
g(1,2).stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
g(1,2).axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 3]);


g(1,3).stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g(1,3).geom_jitter('width',.1,'dodge',1,'alpha',.5);
g(1,3).stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
g(1,3).axe_property('LineWidth',1.5,'FontSize',12,'YLim',[0 3]);

g(1,1).set_names('x','','y','Normalized Volume','color','Groups');
g(1,2).set_names('x','','y','Branch Points Per Slice','color','Groups');
g(1,3).set_names('x','','y','End Points Per Slice','color','Groups');

g(1,1).set_order_options('x',{'SS','MS','MN'});
g(1,2).set_order_options('x',{'SS','MS','MN'});
g(1,3).set_order_options('x',{'SS','MS','MN'});

figure('Position',[100 100 1200 400]);
g.draw();

export_fig('Naloxone Morphology.png','-m5');

% 2 Way ANOVA
[V_p,V_t2,V_stats] = anovan(G.fun1_CellVolume./Slice',{G.A1 G.A2},'model','interaction','varnames',{'Drug','Treatment'});
figure;
[Volume_c,~,~,~] = multcompare(V_stats,'Dimension',[1 2],'CType','dunn-sidak');

[B_p,B_t2,B_stats] = anovan(G.fun1_numbranchpts./Slice',{G.A1 G.A2},'model','interaction','varnames',{'Drug','Treatment'});
figure;
[Branch_c,~,~,~] = multcompare(B_stats,'Dimension',[1 2],'CType','dunn-sidak');

[E_p,E_t2,E_stats] = anovan(G.fun1_numendpts./Slice',{G.A1 G.A2},'model','interaction','varnames',{'Drug','Treatment'});
figure;
[End_c,~,~,~] = multcompare(E_stats,'Dimension',[1 2],'CType','dunn-sidak');

save('Morphology_Stats','V_p','V_t2','V_stats','Volume_c','B_p','B_t2','B_stats','Branch_c','E_p','E_t2','E_stats','End_c');

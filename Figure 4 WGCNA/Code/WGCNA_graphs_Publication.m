%% Analysing WGCNA output using custom w (wgcna) class (Kevin Coffey & Russell Marx - 2019)
% You must be in the main code directory to run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ALWAYS BACKUP YOUR WGCNA OUTPUT FILES BEFORE RUNNING! %%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run IP Samples
clear all

%% Initialize class and load WGCNA Output
w = WGCNA;

%% Run on first pass to save better colors and smaller 'dissTOM.mat' 
% % w.loadGeneTable('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\Gene Modules IP.csv','coolColors', 1);
% % w.loadDissTOM('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\dissTOM IP.csv','saveAsMat', 1);

%% Run on all subsequent  
w.loadGeneTable('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\Gene Modules IP.csv','coolColors', 0);
w.loadDissTOM('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\dissTOM IP.mat','saveAsMat', 0);

%% Figure and Results Directory
cd('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\IP Results');

%% Plot the eigen-genes
w.calculateEigenGenes
w.plotEigenGenes
export_fig('eigenGenes_IP.png','-m3')

%% Merge Eigen-genes
% Don't save colors until you are happy with merging
% Don't run at all if the modules are 
% % w.mergeEigenGenes('mergeThreshold', 1,'saveTable',1);
% % w.calculateEigenGenes
% % w.plotEigenGenes
% % export_fig('eigenGenes_merged.png','-m3')

%% Module Significance w/ DESEQ2
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_2_Morphine Vehicle _Saline Vehicle_IP.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Morphine Module Significance_IP.png','-m3')

w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_3_Morphine Naloxone_Morphine Vehicle_IP.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Withdrawal Module Significance_IP.png','-m3')

w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Data Sheets\DESeq2_9_IP Adj_IN.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Enr Module Significance_IP.png','-m3')

%% MinSpanTree - Module Colors & Morphine Significance
w.loadDESEQ("\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_2_Morphine Vehicle _Saline Vehicle_IP.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 1000 1000]);
h = plot(g)
layout(h,'force','Iterations',120,'WeightEffect','inverse','UseGravity','on');
[colorNames, rgb] = colornames(w.colorMap,categories(removecats(g.Nodes.module)));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
export_fig('MinSpan Module Color_IP.png','-m5')
% Add Morphine Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Morphine Significance_IP.png','-m5')

% % 3D Layout for Dope Rotation
% figure('Color','k','Position',[1 1 850 900]);
% h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
% layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
% [colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
% set(gcf,'Colormap',rgb);
% h.NodeCData = findgroups(g.Nodes.module);
% axis off
% box off
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=30;OptionX.Duration=15;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Module Color_IP',OptionX)
% % Add Morphine Significance
% h.NodeCData = g.Nodes.P_adj;
% set(gcf,'Colormap',flipud(inferno),'Color','k');
% caxis([0 .5]);
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=30;OptionX.Duration=20;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Morphine Significance_IP',OptionX)
% 
%% MinSpanTree - Module Colors & Withdrawal Significance
w.loadDESEQ("\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_3_Morphine Naloxone_Morphine Vehicle_IP.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 1000 1000]);
h = plot(g)
layout(h,'force','Iterations',120,'WeightEffect','inverse','UseGravity','on');
[colorNames, rgb] = colornames(w.colorMap,categories(removecats(g.Nodes.module)));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Withdrawal Significance_IP.png','-m5')

% % 3D Layout for Dope Rotation
% figure('Color','k','Position',[1 1 850 900]);
% h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
% layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
% [colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
% set(gcf,'Colormap',rgb);
% h.NodeCData = findgroups(g.Nodes.module);
% axis off
% box off
% % Add Withdrawal Significance
% h.NodeCData = g.Nodes.P_adj;
% set(gcf,'Colormap',flipud(inferno),'Color','k');
% caxis([0 .5]);
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=30;OptionX.Duration=20;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Withdrawal Significance_IP',OptionX)

%% MinSpanTree - Module Colors & Enrichment Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Data Sheets\DESeq2_9_IP Adj_IN.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 850 900]);
h = plot(g)
layout(h,'force','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Enrichment Significance_IP.png','-m5')

% 3D Layout for Dope Rotation
figure('Color','k','Position',[1 1 850 900]);
h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
% Set up recording parameters (optional), and record
OptionX.FrameRate=30;OptionX.Duration=20;OptionX.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Enrichment Significance_IP',OptionX)

%% Plot Module Network for modules of interest Morphine; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_2_Morphine Vehicle _Saline Vehicle_IP.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=8;
OptionZ.SigNode=1;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

% Centrality
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center,(g.Nodes.Wald_Stats),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Morphine Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Morphine.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Morphine.xlsx']);
end

%% Plot Module Network for modules of interest Withdrawal; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_3_Morphine Naloxone_Morphine Vehicle_IP.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=8;
OptionZ.SigNode=1;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

% Centrality
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center,(g.Nodes.Wald_Stats),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Morphine Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Withdrawal.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Withdrawal.xlsx']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Run IP Minus Noise Samples
clear all

%% Initialize class and load WGCNA Output
w = WGCNA;

%% Run on first pass to save better colors and smaller 'dissTOM.mat' 
% % w.loadGeneTable('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\Gene Modules IP Minus.csv','coolColors', 1);
% % w.loadDissTOM('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\dissTOM IP Minus.csv','saveAsMat', 1);

%% Run on all subsequent  
w.loadGeneTable('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\Gene Modules IP Minus.csv','coolColors', 0);
w.loadDissTOM('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\dissTOM IP Minus.mat','saveAsMat', 0);

%% Figure and Results Directory
cd('\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\IP Minus Noise Results');

%% Plot the eigen-genes
w.calculateEigenGenes
w.plotEigenGenes
export_fig('eigenGenes_IP Minus.png','-m3')

%% Merge Eigen-genes
% Don't save colors until you are happy with merging
% Don't run at all if the modules are 
% % w.mergeEigenGenes('mergeThreshold', .85,'saveTable',1);
% % w.calculateEigenGenes
% % w.plotEigenGenes
% % export_fig('eigenGenes merged_IP Minus.png','-m3')

%% Module Significance w/ DESEQ2
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_7_Morphine Vehicle _Saline Vehicle_IP_Adj.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Morphine Module Significance_IP Minus.png','-m3')

w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_8_Morphine Naloxone_Morphine Vehicle_IP_Adj.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Withdrawal Module Significance_IP Minus.png','-m3')

w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Data Sheets\DESeq2_9_IP Adj_IN.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Enr Module Significance_IP Minus.png','-m3')

%% MinSpanTree - Module Colors & Morphine Significance
w.loadDESEQ("\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_7_Morphine Vehicle _Saline Vehicle_IP_Adj.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 1000 1000]);
h = plot(g)
layout(h,'force','Iterations',120,'WeightEffect','inverse','UseGravity','on');
[colorNames, rgb] = colornames(w.colorMap,categories(removecats(g.Nodes.module)));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
export_fig('MinSpan Module Color_IP Minus.png','-m5')
% Add Morphine Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Morphine Significance_IP Minus.png','-m5')

% 3D Layout for Dope Rotation
% figure('Color','k','Position',[1 1 1080 1080]);
% h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
% layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
% [colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
% set(gcf,'Colormap',rgb);
% h.NodeCData = findgroups(g.Nodes.module);
% axis off
% box off
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=60;OptionX.Duration=20;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Module Color_IP Minus',OptionX)
% % Add Morphine Significance
% h.NodeCData = g.Nodes.P_adj;
% set(gcf,'Colormap',flipud(inferno),'Color','k');
% caxis([0 .5]);
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=60;OptionX.Duration=20;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Morphine Significance_IP Minus',OptionX)

%% MinSpanTree - Module Colors & Withdrawal Significance
w.loadDESEQ("\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_8_Morphine Naloxone_Morphine Vehicle_IP_Adj.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 1000 1000]);
h = plot(g)
layout(h,'force','Iterations',120,'WeightEffect','inverse','UseGravity','on');
[colorNames, rgb] = colornames(w.colorMap,categories(removecats(g.Nodes.module)));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Withdrawal Significance_IP Minus.png','-m5')

% % 3D Layout for Dope Rotation
% figure('Color','k','Position',[1 1 1080 1080]);
% h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
% layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
% [colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
% set(gcf,'Colormap',rgb);
% h.NodeCData = findgroups(g.Nodes.module);
% axis off
% box off
% % Add Withdrawal Significance
% h.NodeCData = g.Nodes.P_adj;
% set(gcf,'Colormap',flipud(inferno),'Color','k');
% caxis([0 .5]);
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=60;OptionX.Duration=20;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Withdrawal Significance_IP Minus',OptionX)

%% MinSpanTree - Module Colors & Enrichment Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Enrichment\DESeq2_IP_Adj_Input - Wave 2.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 850 900]);
h = plot(g)
layout(h,'force','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Enrichment Significance_IP Minus.png','-m5')

% 3D Layout for Dope Rotation
figure('Color','k','Position',[1 1 1080 1080]);
h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
% Set up recording parameters (optional), and record
OptionX.FrameRate=60;OptionX.Duration=20;OptionX.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Enrichment Significance_IP Minus',OptionX)

%% Plot Module Network for modules of interest Morphine; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_7_Morphine Vehicle _Saline Vehicle_IP_Adj.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);

for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=10;
OptionZ.SigNode=1;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

OptionZ.EdgeAlpha=1;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=24;
OptionZ.SigNode=0;
OptionZ.Layout='tree';
g = w.getGraphOfModule(Interest{i,1});
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Tree.png'],'-m5')

% Centrality
idx=~isnan(g.Nodes.Wald_Stats);
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center(idx),(g.Nodes.Wald_Stats(idx)),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Morphine Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Morphine.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Morphine.xlsx']);
end

%% Plot Module Network for modules of interest Withdrawal; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("\\128.95.12.244\kcoffey\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_8_Morphine Naloxone_Morphine Vehicle_IP_Adj.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=8;
OptionZ.SigNode=0;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

OptionZ.EdgeAlpha=1;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=24;
OptionZ.SigNode=0;
OptionZ.Layout='tree';
g = w.getGraphOfModule(Interest{i,1});
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Tree.png'],'-m5')

% Centrality
idx=~isnan(g.Nodes.Wald_Stats);
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center(idx),(g.Nodes.Wald_Stats(idx)),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Withdrawal Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Withdrawal.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Withdrawal.xlsx']);
end

%% Plot Module Network for modules of interest Withdrawal; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Enrichment\DESeq2_IP_Adj_Input - Wave 2.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=8;
OptionZ.SigNode=0;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=12;
OptionZ.SigNode=1;
OptionZ.Layout='force';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Force.png'],'-m5')

% Centrality
idx=~isnan(g.Nodes.Wald_Stats);
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center(idx),(g.Nodes.Wald_Stats(idx)),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Morphine Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Enrichment.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Enrichment.xlsx']);
end

%% Plot graph of individual gene network
g = w.getGraphOfNeighbors('Cdk17', 300);
g = w.pruneEdges(g,50)
h = w.plotGraph(g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Run Input Samples
clear all

%% Initialize class and load WGCNA Output
w = WGCNA;

%% Run on first pass to save better colors and smaller 'dissTOM.mat' 
% % w.loadGeneTable('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\Gene Modules IN.csv','coolColors', 1);
% % w.loadDissTOM('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\dissTOM IN.csv','saveAsMat', 1);

%% Run on all subsequent  
w.loadGeneTable('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\Gene Modules IP Minus.csv','coolColors', 0);
w.loadDissTOM('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\WGCNA Output\dissTOM IP.mat','saveAsMat', 0);

%% Figure and Results Directory
cd('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Input Results');

%% Plot the eigen-genes
w.calculateEigenGenes
w.plotEigenGenes
export_fig('eigenGenes_IN.png','-m3')

%% Merge Eigen-genes
% Don't save colors until you are happy with merging
% Don't run at all if the modules are 
% % w.mergeEigenGenes('mergeThreshold', 1,'saveTable',1);
% % w.calculateEigenGenes
% % w.plotEigenGenes
% % export_fig('eigenGenes merged_IN.png','-m3')

%% Module Significance w/ DESEQ2
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_5_Morphine Vehicle_Saline Vehicle_Input.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Morphine Module Significance_IN.png','-m3')

w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_6_Morphine Naloxone_Morphine Vehicle_Input.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Withdrawal Module Significance_IN.png','-m3')

%% MinSpanTree - Module Colors & Morphine Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_5_Morphine Vehicle_Saline Vehicle_Input.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 850 900]);
h = plot(g)
layout(h,'force','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
export_fig('MinSpan Module Color_IN.png','-m5')
% Add Morphine Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Morphine Significance_IN.png','-m5')

% 3D Layout for Dope Rotation
figure('Color','k','Position',[1 1 850 900]);
h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Set up recording parameters (optional), and record
OptionX.FrameRate=30;OptionX.Duration=15;OptionX.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Module Color_IN',OptionX)
% Add Morphine Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
% Set up recording parameters (optional), and record
OptionX.FrameRate=30;OptionX.Duration=20;OptionX.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Morphine Significance_IN',OptionX)

%% MinSpanTree - Module Colors & Withdrawal Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_6_Morphine Naloxone_Morphine Vehicle_Input.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM);
g = w.removeDisconnectedNodes(g);
figure('Color','k','Position',[1 1 850 900]);
h = plot(g)
layout(h,'force','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
export_fig('MinSpan Withdrawal Significance_IN.png','-m5')

% 3D Layout for Dope Rotation
figure('Color','k','Position',[1 1 850 900]);
h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
[colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
% Add Withdrawal Significance
h.NodeCData = g.Nodes.P_adj;
set(gcf,'Colormap',flipud(inferno),'Color','k');
caxis([0 .5]);
% Set up recording parameters (optional), and record
OptionX.FrameRate=30;OptionX.Duration=20;OptionX.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Withdrawal Significance_IN',OptionX)

%% Plot Module Network for modules of interest Morphine; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_5_Morphine Vehicle_Saline Vehicle_Input.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=8;
OptionZ.SigNode=0;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

% Centrality
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center,(g.Nodes.Wald_Stats),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Morphine Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Morphine.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Morphine.xlsx']);
end

%% Plot Module Network for modules of interest Withdrawal; & Centrality VS Significance;
% Centrality VS Significance
w.loadDESEQ("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Morphine and Withdrawal\DESeq2_6_Morphine Naloxone_Morphine Vehicle_Input.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
g = w.pruneEdges(g,3)
g = w.removeDisconnectedNodes(g);
OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=8;
OptionZ.SigNode=0;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
export_fig([Interest{i} ' Network Circle.png'],'-m5')

% Centrality
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
center=centrality(g,'closeness','Cost',g.Edges.Weight);
g.Nodes.Centrality=center;
f1=figure('color','w','position',[100 100 450 400]);
scatter(center,(g.Nodes.Wald_Stats),6,'filled','o');
p = polyfit(center,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,center);
hold on
scatter(center,y1,4,'filled','o');
ylabel('Morphine Wald Statistic');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=center(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq) ', P= ' num2str(Pval)]);
export_fig([Interest{i} ' Centrality Withdrawal.png'],'-m3');
close all

% Save Table of Module of Interest
g.Nodes=sortrows(g.Nodes,'Centrality','descend');
writetable(g.Nodes,[Interest{i} ' Withdrawal.xlsx']);
end
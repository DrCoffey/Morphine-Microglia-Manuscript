%% DeSeq Graphing Figures for Enrichment
% Sets Directory For Data & Figures
clear all close all

% Generate cmap
% lightness, chroma, hue range
lightness = [65, 65];
chroma = [75, 75];
hue = [205 385];

% colormap resolution
n = 20;
LHC = [
    linspace(lightness(1),lightness(2),n)
    linspace(chroma(1),chroma(2),n)
    linspace(hue(1),hue(2),n)
    ]';
cmap = pa_LCH2RGB(LHC);

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 600 400]);
tmp=importdata('DESeq2_9_IP Adj_IN.xlsx'); % Import seq data
filename=strtok('DESeq2_9_IP Adj_IN.xlsx','.')    
C = strsplit(filename,'_');
tmp.textdata=tmp.textdata(2:end,1);
axes(gca);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
hold on;
scatter(log2(tmp.data(:,1)),tmp.data(:,2),6,[.5 .5 .5],'filled','o');
idx=(tmp.data(:,6))<.10 & (tmp.data(:,4))>0;
scatter(log2(tmp.data(idx,1)),tmp.data(idx,2),6,tmp.data(idx,6),'filled','o','markerfacealpha',.6);
colormap(flipud(cmap))
caxis([0 .1])
ylim([-4 8]);
xlim([-5 15]);
yticks([-4 -2 0 2 4 6 8])
yticklabels({'-4','-2','0','2','4','6','8'})
ylabel('log_2(Fold Change)');
text(0,8,'Enrichment','FontSize',12)
xlabel('log_2(TPM)');
colorbar;

% Adding Bubble Plot
[B,index] = sortrows(tmp.data,4,'descend');
sgene=tmp.textdata(index);
idx=isnan(B(:,4))==0;
B=B(idx,:);
sgene=sgene(idx);
B=B(1:15,:);
sgene=sgene(1:15);
hold on
bp=bubbleplot(log2(B(:,1)),B(:,2), [],(B(:,4)), flipud([1:1:15]'), [],'ColorMap',@spring); 
h=legend(bp,cellstr(sgene),'Location','eastoutside');
h.Box='off';
export_fig('Enrichment IP Adj.png', '-m5'); % Save the Figure


%% GSEA for enrichment

t=readtable("MolecularFunctionENR_ADJ_top10.xlsx")
enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend')
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend')

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 600 400]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
[h,icons]=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');

% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
for i=1:length(icons)
set(icons(i),'MarkerSize',8);
end

export_fig('Enrichment MF GSEA Adj.png', '-m5'); % Save the Figure

%% DeSeq Table PCA

inT=readtable("..\Figure 6-8 WGCNA\Data Sheets\GeneTableIN.csv");
ipT=readtable("..\Figure 6-8 WGCNA\Data Sheets\GeneTableIP_Minus.csv");

desT=[inT ipT(:,2:end)];

Groups=categorical({'IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN','IN',...
                    'IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP','IP'})'

GeneTableVariance=desT(var(desT{:,2:end}')>50,:);
Y = tsne(GeneTableVariance{:,2:end}','Algorithm','exact','Distance','cosine','Perplexity',5);
PC1=Y(:,1);
PC2=Y(:,2);
pcTable=table(Groups,PC1,PC2);
f1=figure('color','w','position',[100 100 280 200]);
g=gramm('x',pcTable.PC1,'y',pcTable.PC2,'color',pcTable.Groups)
g.geom_point();
g.set_names('x','TSNE1','y','TSNE2','color','Fraction')
g.axe_property('FontSize',12,'LineWidth',1.5,'TickDir','out');
g.draw;
g.export('file_name','IP-IN-TSNE.png','file_type','png');

%% Gene Panel For **** Reviewer IP
desT=readtable("..\Figure 6-8 WGCNA\Data Sheets\GeneTableIP_Minus.csv");

Group=categorical({'SS','SS','SS','SS','SN','SN','SN','SN','SN','MS','MS','MS','MS','MN','MN','MN','MN'})';
panel=readtable("41598_2018_27293_MOESM2_ESM-topranked.xlsx",'Sheet','top_mouse_enrich');

GeneTableVariance=desT(var(desT{:,2:end}')>100,:);
Y = tsne(GeneTableVariance{:,2:end}','Algorithm','exact','Distance','cosine','Perplexity',6);
PC1=Y(:,1);
PC2=Y(:,2);
pcTable=table(Group,PC1,PC2);
f1=figure('color','w','position',[100 100 300 200]);
g=gramm('x',pcTable.PC1,'y',pcTable.PC2,'color',pcTable.Group)
g.geom_point();
g.set_names('x','TSNE1','y','TSNE2','color','Fraction')
g.axe_property('FontSize',12,'LineWidth',1.5,'TickDir','out');
g.draw;
g.export('file_name','IP-Group-TSNE.png','file_type','png');

for i=1:height(panel)
    disp(i);
    TF = matches(desT.Genes(:),panel.gene{i},'IgnoreCase',true);
    if sum(TF)==1;
    Count= desT{TF,2:18}';
    Gene=repmat(categorical(desT.Genes((TF))),[17,1]);
    CellType=repmat(categorical(panel.Celltype(i)),[17,1]);
    if i==1
       mt=table(Group,Gene,Count,CellType); 
    else
       mt=[mt;table(Group,Gene,Count,CellType)]; 
    end
    end
end

clear g
f1=figure('color','w','position',[100 100 400 300]);
g=gramm('x',mt.CellType,'y',(mt.Count),'color',mt.Group);
%Boxplots
%g.stat_violin('normalization','width','half',0,'dodge',.75,'fill','transparent')
%g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'point','errorbar'},'type','sem','dodge',.75);
g.axe_property('LineWidth',1.5,'FontSize',12,'Ylim',[0 250]);
g.set_names('x','Presumed Cell Markers','y','Normalized Count','color','Groups');
g.set_order_options('x',{'mic','neu','oli','ast','end','opc'});
g.draw;
g.export('file_name','IP-Panel.png','file_type','png');


%% INPUT
desT=readtable("..\Figure 6-8 WGCNA\Data Sheets\GeneTableIN.csv");

Group=categorical({'SS','SS','SS','SS','SN','SN','SN','SN','SN','MS','MS','MS','MS','MN','MN','MN','MN'})';

GeneTableVariance=desT(var(desT{:,2:end}')>100,:);
Y = tsne(GeneTableVariance{:,2:end}','Algorithm','exact','Distance','cosine','Perplexity',6);
PC1=Y(:,1);
PC2=Y(:,2);
pcTable=table(Group,PC1,PC2);
f1=figure('color','w','position',[100 100 280 200]);
g=gramm('x',pcTable.PC1,'y',pcTable.PC2,'color',pcTable.Group)
g.geom_point();
g.set_names('x','TSNE1','y','TSNE2','color','Fraction')
g.axe_property('FontSize',12,'LineWidth',1.5,'TickDir','out');
g.draw;
g.export('file_name','IN-Group-TSNE.png','file_type','png');

for i=1:height(panel)
    disp(i);
    TF = matches(desT.Genes(:),panel.gene{i},'IgnoreCase',true);
    if sum(TF)==1;
    Count= desT{TF,2:18}';
    Gene=repmat(categorical(desT.Genes((TF))),[17,1]);
    CellType=repmat(categorical(panel.Celltype(i)),[17,1]);
    if i==1
       mt=table(Group,Gene,Count,CellType); 
    else
       mt=[mt;table(Group,Gene,Count,CellType)]; 
    end
    end
end


clear g
f1=figure('color','w','position',[100 100 400 300]);
g=gramm('x',mt.CellType,'y',(mt.Count),'color',mt.Group);
%Boxplots
%g.stat_violin('normalization','width','half',0,'dodge',.75,'fill','transparent')
%g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'point','errorbar'},'type','sem','dodge',.75);
g.axe_property('LineWidth',1.5,'FontSize',12,'Ylim',[0 100]);
g.set_names('x','Presumed Cell Markers','y','Normalized Count','color','Groups');
g.set_order_options('x',{'mic','neu','oli','ast','end','opc'});
g.draw;
g.export('file_name','IN-Panel.png','file_type','png');



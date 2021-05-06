%% GSEA for enrichment
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';

t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Enrichment\MolecularFunctionENR.xlsx")
enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend')
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend')

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 800 700]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');

pngFileName = 'Enrichment MF GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure


%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\IP\Morphine\Morphine IP BP.xlsx")
enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Morphine IP Bio. Process';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Morphine IP BP GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure


%% GSEA for Morphine IP 
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\IP\Morphine\Morphine IP RT.xlsx")
enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Morphine IP Reactome';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Morphine IP RT GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\IP\Morphine\Morphine IP TF.xlsx")
enr=t(t.NES>0,:);
try
enr = sortrows(enr,7,'ascend');
end
    
denr=t(t.NES<0,:);
try
denr = sortrows(denr,7,'ascend');
end

if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.GeneSet{i};
    if length(cellContents)>26
       enr.GeneSet{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    else
       enr.GeneSet{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    end
    
end

for i=1:height(denr)
    cellContents = denr.GeneSet{i};
    if length(cellContents)>26
       denr.GeneSet{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    else
       denr.GeneSet{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')'];
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
try
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
end
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.GeneSet; denr.GeneSet]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Morphine IP Trancription Factors';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Morphine IP TF GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure




%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';

t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\IP\Withdrawal\Withdrawal IP BP.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);

enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end


for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 12]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 12],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Withdrawal IP Bio. Process';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Withdrawal IP BP GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure


%% GSEA for Morphine IP 
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';

t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\IP\Withdrawal\Withdrawal IP RT.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);

enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12
   enr=enr(1:12,:);
end
if height(denr)>12
denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Withdrawal IP Reactome';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Withdrawal IP RT GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\IP\Withdrawal\Withdrawal IP TF.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);

enr=t(t.NES>0,:);
try
enr = sortrows(enr,7,'ascend');
end

denr=t(t.NES<0,:);
denr = sortrows(denr,7,'ascend');
if height(enr)>12
   enr=enr(1:12,:);
end


if height(denr)>12
    try
        denr=denr(1:(24-height(enr)),:);
    end
end

for i=1:height(enr)
    cellContents = enr.GeneSet{i};
    if length(cellContents)>26
       enr.GeneSet{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    else
       enr.GeneSet{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    end
    
end

for i=1:height(denr)
    cellContents = denr.GeneSet{i};
    if length(cellContents)>26
       denr.GeneSet{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    else
       denr.GeneSet{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')'];
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
try
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
end
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.GeneSet; denr.GeneSet]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Withdrawal IP Trancription Factors';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Withdrawal IP TF GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';

t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Input\Morphine\Morphine Input BP.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);
enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
try
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
end
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Morphine IN Bio. Process';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Morphine IN BP GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure


%% GSEA for Morphine IP 
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Input\Morphine\Morphine Input RT.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);
enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 12]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 12],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Morphine IN Reactome';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Morphine IN RT GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Input\Morphine\Morphine Input TF.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);
enr=t(t.NES>0,:);
try
enr = sortrows(enr,7,'ascend');
end
    
denr=t(t.NES<0,:);
try
denr = sortrows(denr,7,'ascend');
end

if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.GeneSet{i};
    if length(cellContents)>26
       enr.GeneSet{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    else
       enr.GeneSet{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    end
    
end

for i=1:height(denr)
    cellContents = denr.GeneSet{i};
    if length(cellContents)>26
       denr.GeneSet{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    else
       denr.GeneSet{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')'];
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
try
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
end
ylim([0 12]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 12],'--k');
h=legend(gca,cellstr([enr.GeneSet; denr.GeneSet]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Morphine IN Trancription Factors';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Morphine IN TF GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure




%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';

t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Input\Withdrawal\Withdrawal Input BP.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);

enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12 & height(denr)>12
   enr=enr(1:12,:);
   denr=denr(1:12,:);
elseif height(enr)>12 & height(enr)+height(denr)>24
   enr=enr(1:(24-height(denr)),:);
elseif height(denr)>12 & height(enr)+height(denr)>24  
   denr=denr(1:(24-height(enr)),:);
end


for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 12]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 12],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Withdrawal IN Bio. Process';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Withdrawal IN BP GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure


%% GSEA for Morphine IP 
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';

t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Input\Withdrawal\Withdrawal Input RT.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);

enr=t(t.NES>0,:);
enr = sortrows(enr,8,'ascend');
denr=t(t.NES<0,:);
denr = sortrows(denr,8,'ascend');
if height(enr)>12
   enr=enr(1:12,:);
end
if height(denr)>12
denr=denr(1:(24-height(enr)),:);
end

for i=1:height(enr)
    cellContents = enr.Description{i};
    if length(cellContents)>26
       enr.Description{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
    else
       enr.Description{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
    end
    
end

for i=1:height(denr)
    cellContents = denr.Description{i};
    if length(cellContents)>26
       denr.Description{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
    else
       denr.Description{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')']; 
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
ylim([0 12]);
xlim([-2.5 2.5]);
hold on
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
plot([0 0],[0 12],'--k');
h=legend(gca,cellstr([enr.Description; denr.Description]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Withdrawal IN Reactome';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Withdrawal IN RT GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

%% GSEA for Morphine IP
clear all
fig_fold='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Figures';


t=readtable("Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\GSEA\Input\Withdrawal\Withdrawal Input TF.xlsx")
idx=t.FDR==0;
t.FDR(idx) = 10^-4./t.NES(idx);

enr=t(t.NES>0,:);
try
enr = sortrows(enr,7,'ascend');
end

denr=t(t.NES<0,:);
denr = sortrows(denr,7,'ascend');
if height(enr)>12
   enr=enr(1:12,:);
end


if height(denr)>12
    try
        denr=denr(1:(24-height(enr)),:);
    end
end

for i=1:height(enr)
    cellContents = enr.GeneSet{i};
    if length(cellContents)>26
       enr.GeneSet{i} = [cellContents(1:23) '... (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    else
       enr.GeneSet{i} = [cellContents(1:end) ' (' num2str(enr.Size(i)) ')']; 
       enr.GeneSet{i} = strrep(enr.GeneSet{i},'_','-');
    end
    
end

for i=1:height(denr)
    cellContents = denr.GeneSet{i};
    if length(cellContents)>26
       denr.GeneSet{i} = [cellContents(1:23) '... (' num2str(denr.Size(i)) ')']; 
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    else
       denr.GeneSet{i} = [cellContents(1:end) ' (' num2str(denr.Size(i)) ')'];
       denr.GeneSet{i} = strrep(denr.GeneSet{i},'_','-');
    end
end

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 450 500]);
try
bp=bubbleplot(enr.NES,-log(enr.FDR),[],enr.Size, 1-(enr.FDR), [],'ColorMap',@spring); 
end
ylim([0 8]);
xlim([-2.5 2.5]);
hold on
try
bp2=bubbleplot(denr.NES,-log(denr.FDR),[],denr.Size, (denr.FDR), [],'ColorMap',@cool); 
end
plot([0 0],[0 8],'--k');
h=legend(gca,cellstr([enr.GeneSet; denr.GeneSet]),'Location','eastoutside');
h.Box='off';
xlabel('Normalized Enrichment');
ylabel('-log(FDR)');
title({'Withdrawal IN Trancription Factors';'FDR<.05'},'FontWeight','normal','FontSize',8);

pngFileName = 'Withdrawal IN TF GSEA.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure
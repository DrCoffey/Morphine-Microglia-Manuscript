%% Morphine Naloxone DREADDs Behavior Figures with Gramm(You must have Gramm installed!)
% Navigate to the folder containing the code!

t=readtable("Statistics-Morphine CX3-HM4Di-Ribotag.xlsx");

%% Normalize to Recording Time

t.Distance=t.Distance./t.RelativeLength;
t.Normal=t.Normal./t.RelativeLength;
t.Contracted=(t.Contracted./t.RelativeLength);
t.Moving=t.Moving./t.RelativeLength;
t.Immobility=t.Immobility./t.RelativeLength;
t.Group=categorical(t.Group);
t.Genotype=categorical(t.Genotype);
t.Treatment=categorical(t.Treatment);

g=gramm('x',t.Genotype,'y',t.Distance,'color',t.Treatment,'Subset',t.Group=='Morphine Naloxone');
g.stat_violin('normalization','width','dodge',1,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
figure('Position',[100 100 325 300]);
g.set_names('x',[],'y','Distance (cm)');
%g.set_title('Morphine Naloxone');
g.set_color_options('map','lch','hue_range',[25 385]+180);
g.axe_property('ylim',[0 4000],'LineWidth',1.5,'FontSize',12);
g.set_order_options('color',{'Vehicle','CNO'});
g.draw();

g=gramm('x',t.Genotype,'y',t.Contracted,'color',t.Treatment,'Subset',t.Group=='Morphine Naloxone');
g.stat_violin('normalization','width','dodge',1,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
figure('Position',[100 100 325 300]);
g.set_names('x',[],'y','Contracted (s)');
%g.set_title('Morphine Naloxone');
g.set_color_options('map','lch','hue_range',[25 385]+180);
g.axe_property('ylim',[0 2000],'LineWidth',1.5,'FontSize',12);
g.set_order_options('color',{'Vehicle','CNO'});
g.draw();

g=gramm('x',t.Genotype,'y',t.Immobility,'color',t.Treatment,'Subset',t.Group=='Morphine Naloxone');
g.stat_violin('normalization','width','dodge',1,'fill','transparent')
g.geom_jitter('width',.1,'dodge',1,'alpha',.5);
g.stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
figure('Position',[100 100 325 300]);
g.set_names('x',[],'y','Immobility (s)');
%g.set_title('Naloxone');
g.set_color_options('map','lch','hue_range',[25 385]+180);
g.axe_property('ylim',[0 2000],'LineWidth',1.5,'FontSize',12);
g.set_order_options('color',{'Vehicle','CNO'});
g.draw();

t2=t(t.Group=='Morphine Naloxone',:);

[p,percentExpress_Anova,stats] = anova1([t2.Distance(t2.Genotype=='Ribotag') t2.Distance(t2.Genotype=='HM4Di' & t2.Treatment=='Vehicle') t2.Distance(t2.Genotype=='HM4Di' & t2.Treatment=='CNO')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');
save('Dist_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');


[p,percentExpress_Anova,stats] = anova1([t2.Contracted(t2.Genotype=='Ribotag') t2.Contracted(t2.Genotype=='HM4Di' & t2.Treatment=='Vehicle') t2.Contracted(t2.Genotype=='HM4Di' & t2.Treatment=='CNO')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');
save('Contract_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');

[p,percentExpress_Anova,stats] = anova1([t2.Immobility(t2.Genotype=='Ribotag') t2.Immobility(t2.Genotype=='HM4Di' & t2.Treatment=='Vehicle') t2.Immobility(t2.Genotype=='HM4Di' & t2.Treatment=='CNO')]);
figure;
[percentExpress_Mult,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');
save('Immobile_stats.m','percentExpress_Anova','percentExpress_Mult','normColoc_Anova','normColoc_Mult');





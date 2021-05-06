%% Morphine Naloxone Behavior Figures with Gramm(You must have Gramm installed!)
% Navigate to the folder containing the code!

close all 
clear all

t=readtable('NaloxoneBehavior.xlsx');
t.Treatment=categorical(t.Treatment);

clear g

g(1,1)=gramm('x',t.Label,'y',t.Distance,'color',t.Treatment);
g(1,2)=gramm('x',t.Label,'y',t.Contracted,'color',t.Treatment);
g(1,3)=gramm('x',t.Label,'y',t.Immobile,'color',t.Treatment);

%Boxplots
g(1,1).stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g(1,1).geom_jitter('width',.1,'dodge',1,'alpha',.5);
g(1,1).stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
g(1,1).axe_property('LineWidth',1.5,'FontSize',12);

g(1,2).stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g(1,2).geom_jitter('width',.1,'dodge',1,'alpha',.5);
g(1,2).stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
g(1,2).axe_property('LineWidth',1.5,'FontSize',12);

g(1,3).stat_violin('normalization','width','half',0,'dodge',0,'fill','transparent')
g(1,3).geom_jitter('width',.1,'dodge',1,'alpha',.5);
g(1,3).stat_summary('geom',{'black_errorbar'},'type','sem','dodge',1);
g(1,3).axe_property('LineWidth',1.5,'FontSize',12);


g(1,1).axe_property('YLim',[0 18000]);
g(1,2).axe_property('YLim',[0 1500]);
g(1,3).axe_property('YLim',[0 1500]);

g(1,1).set_names('x','','y','Distance (cm)','color','Groups');
g(1,2).set_names('x','','y','Contraction (s)','color','Groups');
g(1,3).set_names('x','','y','Immobility (cms)','color','Groups');

g(1,1).set_order_options('x',{'SS','SN','MS','MN'});
g(1,2).set_order_options('x',{'SS','SN','MS','MN'});
g(1,3).set_order_options('x',{'SS','SN','MS','MN'});

figure('Position',[100 100 1200 350]);
g.draw();

export_fig('Naloxone Behavior.png','-m5');
%Some methods can be called on all objects at the same time !

[p,t2,stats] = anovan(t.Distance,{t.A1 t.A2},'model','interaction','varnames',{'Drug','Treatment'});
figure;
[Dist_c,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,t2,stats] = anovan(t.Contracted,{t.A1 t.A2},'model','interaction','varnames',{'Drug','Treatment'});
figure;
[Contract_c,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

[p,t2,stats] = anovan(t.Immobile,{t.A1 t.A2},'model','interaction','varnames',{'Drug','Treatment'});
figure;
[Immobile_c,~,~,~] = multcompare(stats,'Dimension',[1 2],'CType','dunn-sidak');

save('Figure_2bcd_stats.m','Dist_c','Contract_c','Immobile_c');
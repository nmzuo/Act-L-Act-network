function plot_fig_format_revise_overlayBoxplot(thr)
%% this function is to re-draw and format the desired figures

outfig='/datc/dynNet/mod_comp_fig_format_09_8_avg/';
outdat='/datc/dynNet/code/data_09_8_avg/';
tshort = {'REST1','GAMBLING','MOTOR','SOCIAL',  ...
         'EMOTION',  'LANGUAGE',  'RELATIONAL', 'WM'};   
tshort2= {'Rest', 'Gambling', 'Motor', 'Social', ...
          'Emotion', 'Language', 'Relational', 'WM'};
      
preout='revise_avg_';
thr=15;
% sex, age
%sex_age = load('hcp_S465_sex_age_RL.txt');
%sex_age = load('hcp_S463_sex_age.txt');
sex_age = load('hcp_S453_sex_age_avg.txt');

cpool=[1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 0.5 0.5 0; 0.5 0 0.5; ...
           0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 1 0 0.5; 0 0.5 1];

       
       
       
%% plot errobar the weit-interaction change between Task-Rest, for 7 tasks
%% the same plot contents fo the above,
%% but we use a more informative plot style
%% https://github.com/raacampbell/notBoxPlot 
if 0
    addpath(genpath('/datc/dynNet/code/notBoxPlot'));
    load([outdat preout 'interActRestTask_thr' num2str(thr) '.mat']); %% 'interActRest', 'interActTask'
    interTask = squeeze(interActTask(:,:,2));
    interRest = squeeze(interActRest(:,:,2));
    [h,p]=ttest(interTask, interRest); disp(p);
    
    xloc = [1:5:35];
    figure('position', [100,100, 900,500]);
%     H1=notBoxPlot(interTask, xloc, 0.5);
%     for i=1:7
%         set([H1(i).sdPtch],	'FaceColor',cpool(i,:), 'EdgeColor','none');
%     end
%     hold on;
%     H1=notBoxPlot(interRest, xloc+1, 0.5);
%     for i=1:7
%         set([H1(i).sdPtch],	'FaceColor',cpool(i,:), 'EdgeColor','none');
%     end
    notBoxPlot(interTask, xloc, 0.5); hold on; boxplot(interTask,'positions', xloc, 'Widths',0.4); hold on;
    set(gca, 'xticklabel', {''});
    notBoxPlot(interRest, xloc+1, 0.5); hold on; boxplot(interRest,'positions', xloc+1, 'Widths',0.4);
    set(gca, 'xticklabel', {''});
    
    xlim([-2,36]);
    
    tickpos=xloc+0.5;
    set(gca, 'xtick',tickpos, 'xticklabel', tshort(2:8)); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',14);
    title(['Interaction in Task and Rest'], 'interpreter', 'none'); 
%    xticklabel_rotate([],30,[]);
    ylim([0.1,0.5]);
    set(gca, 'ytick', [0.15,0.25,0.35, 0.45],  'yticklabel',{'0.15', '0.25', '0.35', '0.45'});
%    xlabel('Act-nonAct pairs separately defined by 7 Tasks', 'Fontsize',18);
    ylabel('Interaction strength between Act and L-Act', 'Fontsize',16);
   
%    legend({'Task', 'Rest'});
    box off;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1100 500]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'interAct_s500_bin_notBoxPlot_thr' num2str(thr) '_v2.tif']);
end


%% plot the totE (Efficiency) between the Rest and Task
%% by notBoxPlot
if 0
    addpath(genpath('/datc/dynNet/code/notBoxPlot'));
    load([outdat preout 'in_out_E_thr' num2str(thr) '.mat']); %% 'totQ', 'totS'
    %% for t-testing
    EG=EG';
    [h, p]=ttest(EG(:,2:8), repmat(EG(:,1),[1,7])); disp(p);
    
    figure('position', [100,100, 900,500]);
    %notboxPlot
%     H1=notBoxPlot(EG', [], 0.5);
%     for i=2:8
%         set([H1(i).sdPtch],	'FaceColor',cpool(i-1,:), 'EdgeColor','none');
%     end
%     set([H1(1).sdPtch],	'FaceColor',cpool(8,:), 'EdgeColor','none');
    notBoxPlot(EG); hold on; boxplot(EG); hold off;
    
    set(gca, 'xtick', [1:8], 'xticklabel', tshort(1:8)); 
%    xticklabel_rotate([],30,[]);
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',12);
    title(['Efficiency comparisons'], 'interpreter', 'none'); 
    set(gca, 'ytick', [0.35,0.45,0.55]);
%    xlabel('Different states: 1 resting-state and 7 task-states', 'Fontsize',14);
    ylabel('Efficiency index', 'Fontsize',18);
    
    box off;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 900 500]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'GlobalEfficiency_s500_bin_notBoxPlot_thr' num2str(thr) '_v2.tif']);
    
end

%% plot the totE (Efficiency) between the Rest and Task
%% by notBoxPlot
if 1
    addpath(genpath('/datc/dynNet/code/notBoxPlot'));
    E=load([outdat preout 'in_out_E_thr' num2str(thr) '.mat']); %% 'E', 'global E'
    EN=load([outdat preout 'in_out_E_thr' num2str(thr) '_normalized.mat']); %% 'E', 'global E'
    %% for t-testing
    EG=E.EG';
    [h, p]=ttest(EG(:,2:8), repmat(EG(:,1),[1,7])); disp('Original');disp(p);
    EG=EN.EG';
    [h, p]=ttest(EG(:,2:8), repmat(EG(:,1),[1,7])); disp('Randomized');disp(p);
    EG=E.EG' ./ EN.EG';
    [h, p]=ttest(EG(:,2:8), repmat(EG(:,1),[1,7])); disp('Normalized');disp(p);
    
    % Randomized
    figure('position', [100,100, 900,500]);
    %notboxPlot
    notBoxPlot(EN.EG'); hold on; boxplot(EN.EG'); hold off;   
    set(gca, 'xtick', [1:8], 'xticklabel', tshort(1:8)); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',12);
    title(['Efficiency comparisons'], 'interpreter', 'none'); 
    set(gca, 'ytick', [0.35,0.45,0.55]);
    ylabel('Efficiency index', 'Fontsize',18);   
    box off;
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 900 500]);
    print(gcf, '-dtiff', '-r300', ['/datc/dynNet/figure_rebuttal/revise_globalE_8state_randomized.tif']);
 
    % Normalized by the Randomized
    figure('position', [100,100, 900,500]);
    %notboxPlot
    notBoxPlot(E.EG' ./ EN.EG'); hold on; boxplot(E.EG' ./ EN.EG'); hold off;   
    set(gca, 'xtick', [1:8], 'xticklabel', tshort(1:8)); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',12);
    title(['Efficiency comparisons'], 'interpreter', 'none'); 
%    set(gca, 'ytick', [0.35,0.45,0.55]);
    ylabel('Efficiency index', 'Fontsize',18);   
    box off;
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 900 500]);
    print(gcf, '-dtiff', '-r300', ['/datc/dynNet/figure_rebuttal/revise_globalE_8state_randomized_normalized.tif']);    
    
    
end

%% Efficiency in Act and nonAct regions, compare between Task and Rest
if 0
    E=load([outdat preout 'in_out_E_thr' num2str(thr) '.mat']);

    for i=1:7 % 1:7
        mIn=mean(E.EGinTask(i,:)'-E.EGinRest(i,:)');    mOut=mean(E.EGoutTask(i,:)'-E.EGoutRest(i,:)');
        [h3, p3]=ttest(E.EGinTask(i,:)', E.EGinRest(i,:)');
        [h4, p4]=ttest(E.EGoutTask(i,:)', E.EGoutRest(i,:)');
        fprintf('in: mean=%f, p=%e;  out: mean=%f, p=%e;',  mIn, p3, mOut, p4);
        fprintf('\n');
    end
    addpath(genpath('/datc/dynNet/code/notBoxPlot'));
    %% IN
    EGTask = E.EGinTask;
    EGRest = E.EGinRest;
    xloc = [1:5:35];
    figure('position', [100,100, 900,500]);
%     H1=notBoxPlot(EGTask', xloc, 0.5);
%     for i=1:7
%         set([H1(i).sdPtch],	'FaceColor',cpool(i,:), 'EdgeColor','none');
%     end
%     hold on;
%     H1=notBoxPlot(EGRest', xloc+1, 0.5);
%     for i=1:7
%         set([H1(i).sdPtch],	'FaceColor',cpool(i,:), 'EdgeColor','none');
%     end
    notBoxPlot(EGTask', xloc, 0.5); hold on; boxplot(EGTask','positions', xloc, 'Widths',0.4); hold on; 
    set(gca, 'xticklabel', {''});
    notBoxPlot(EGRest', xloc+1, 0.5); hold on; boxplot(EGRest','positions', xloc+1, 'Widths',0.4); hold on; 
    set(gca, 'xticklabel', {''});
    
    xlim([-2,36]);
    
    tickpos=xloc+0.5;
    set(gca, 'xtick',tickpos, 'xticklabel',tshort2(2:8), 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',12); 
%    xticklabel_rotate([],30,[]);
%    title(['Interaction in Task and Rest'], 'interpreter', 'none'); 
%    xticklabel_rotate([],30,[]);
    gety=ylim;
%    [y1, y2]=setTicks(gety);
%    set(gca, 'ytick', y1, 'yticklabel', y2);
%    xlabel('Efficiency comparison in Act communities', 'Fontsize',18);
    ylabel('Act: Efficiency',  'FontName', 'Arial', 'FontWeight', 'bold','Fontsize',16);
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
%    legend({'Task', 'Rest'});
    box off;
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1100 500]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'in_out_efficiency_task_rest_IN_thr' num2str(thr) '_v2.tif']);
    
    %% OUT
    EGTask = E.EGoutTask;
    EGRest = E.EGoutRest;
    xloc = [1:5:35];
    figure('position', [100,100, 900,500]);
%     H1=notBoxPlot(EGTask', xloc, 0.5);
%     for i=1:7
%         set([H1(i).sdPtch],	'FaceColor',cpool(i,:), 'EdgeColor','none');
%     end
%     hold on;
%     H1=notBoxPlot(EGRest', xloc+1, 0.5);
%     for i=1:7
%         set([H1(i).sdPtch],	'FaceColor',cpool(i,:), 'EdgeColor','none');
%     end
    notBoxPlot(EGTask', xloc, 0.5); hold on; boxplot(EGTask','positions', xloc, 'Widths',0.4); hold on; 
    set(gca, 'xticklabel', {''});
    notBoxPlot(EGRest', xloc+1, 0.5); hold on; boxplot(EGRest','positions', xloc+1, 'Widths',0.4); hold on; 
    set(gca, 'xticklabel', {''});
    
    xlim([-2,36]);
    
    tickpos=xloc+0.5;
    set(gca, 'xtick',tickpos, 'xticklabel',tshort2(2:8), 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',12); 
%    xticklabel_rotate([],30,[]);
%    title(['Interaction in Task and Rest'], 'interpreter', 'none'); 
%    xticklabel_rotate([],30,[]);
    gety=ylim;
%    [y1, y2]=setTicks(gety);
%    xlabel('Efficiency comparison in Act communities', 'Fontsize',18);
    ylabel('nonAct: Efficiency',  'FontName', 'Arial', 'FontWeight', 'bold','Fontsize',16);
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
%    legend({'Task', 'Rest'});
    box off;
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1100 500]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'in_out_efficiency_task_rest_OUT_thr' num2str(thr) '_v2.tif']);
    

end % efficiency in Act and nonAct
    


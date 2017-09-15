%% this is to show the connectivity matrix by Act/L-Act blocks
%% these figures are shown in supplemental materials.

function revise_show_net_matrix

preout='revise_avg_';
%thr=15;
outdat='/datc/dynNet/code/data_09_8_avg/';
outfig='/datc/dynNet/mod_comp_fig_format_09_8_avg/';


addpath(genpath('/datc/software/NIfTI_20140122')); %% load_nii

load([outdat preout 'corrR_all_remMActivation.mat']);

roiInd=[1:264];
nSubj = size(corrR_all,2); %463
nNode = size(corrR_all,3); %264
nROI = size(corrR_all,3); %264

%tname = {'rfMRI_REST1_LR','tfMRI_GAMBLING_LR','tfMRI_MOTOR_LR','tfMRI_SOCIAL_LR',  ...
%         'tfMRI_EMOTION_LR',  'tfMRI_LANGUAGE_LR',  'tfMRI_RELATIONAL_LR', 'tfMRI_WM_LR'};
tname = {'rfMRI_REST1_avg','tfMRI_GAMBLING_avg','tfMRI_MOTOR_avg','tfMRI_SOCIAL_avg',  ...
         'tfMRI_EMOTION_avg',  'tfMRI_LANGUAGE_avg',  'tfMRI_RELATIONAL_avg', 'tfMRI_WM_avg'};     
tshort = {'REST1','GAMBLING','MOTOR','SOCIAL',  ...
         'EMOTION',  'LANGUAGE',  'RELATIONAL', 'WM'};     
pathbase='/datc/dynNet/groupActi_avg/';
maskname={'/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% gambling
          '/groupmean.gfeat/cope12345_cluster_mask_zstat1.nii.gz', ...    %% motor
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% social
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% emotion
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% language
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% relational
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% wm
         };

coord264 = load('roi264_Power.txt');

interActRest = nan(nSubj,7,2); % 1.ShareMode, 2.connWeit
interActTask = nan(nSubj,7,2);

interModRest = nan(nSubj,1);
interModTask = nan(nSubj,1);

%% The following for-loop is to compute the task-rest task ROI/Rest interactions by 2 methods
%%('Interaction difference between rest and tasks');
% the following is the output of this section

r0Mat=squeeze(corrR_all(1,:,:,:)); r0Mat=squeeze(mean(r0Mat,1));
if 1
    disp('********compute the task-rest task ROI/Rest interactions by 2 methods********');
    for i=1:7   % 1:7   %% for i=1:7
        taskROI = load_nii([pathbase tname{i+1} maskname{i}]); % 3 is for motor
        taskROI = taskROI.img;
        [inInd, outInd] = sepROIByMask(coord264(roiInd,:), taskROI);  
        Nin=length(inInd);  Nout=length(outInd);
        fprintf('%s: [inInd=%d, outInd=%d]\n', tname{i+1}, Nin, Nout);
        newInd=[inInd; outInd];
        tMat = squeeze(corrR_all(i+1,:,:,:));
        tMat=squeeze(mean(tMat,1));
        rMat=r0Mat(newInd, newInd);
        tMat=tMat(newInd, newInd);
        %drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );  
        figure,
        ax(1)=subplot(1,2,1); imshow(rMat,[-0.3, 1]); colormap(redbluecmap);  
        title(tshort{1},'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',14);
        xc=xlim; yc=ylim;
        line([xc(1)-5, xc(1)-30],[yc(2), yc(2)], 'LineWidth', 2, 'clipping', 'off');hold on;
        line([xc(1)-5, xc(1)-30],[yc(1)+length(inInd), yc(1)+length(inInd)], 'LineWidth', 2, 'clipping', 'off');hold on;
        line([xc(1)-5, xc(1)-30],[yc(1), yc(1)], 'LineWidth', 2, 'clipping', 'off');
        text(xc(1)-40, yc(1)+0.5*length(inInd)+20, 'Act','Rotation',90,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',14);
        text(xc(1)-40, yc(2)-0.5*length(outInd)+35, 'nonAct','Rotation',90,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',14);
        %drawArrow([xc(2)-20, xc(2)+60], [yc(1)+0.5*nROI, yc(1)+0.5*nROI],  'LineWidth', 4, 'clipping', 'off');
        line([xc(2)+10, xc(2)+30],[yc(1)+0.5*nROI, yc(1)+0.5*nROI], 'LineWidth', 2, 'clipping', 'off');
        text(xc(2)+25, yc(1)+0.5*nROI-3.5, '\rightarrow', 'fontsize', 32, 'color', 'b');
        
        ax(2)=subplot(1,2,2); imshow(tMat,[-0.3, 1]); colormap(redbluecmap); 
%         s2Pos=get(s2, 'position');
%         hc=colorbar('location', 'eastoutside');
%         set(s2, 'position', s2Pos);
        % Get positions of all the subplot
        posa = get(ax,'position');
        h    = colorbar;

        % Reset ax(3) position to before colorbar
        set(ax(2),'position',posa{2})

        % Set everything to units pixels (avoids dynamic reposition)
        set([ax h],'units','pix')

        % Widen figure by a factor of 1.1 (tweak it for needs)
        posf = get(gcf,'position');
        set(gcf,'position',[posf(1:2) posf(3)*1.1 posf(4)])
        
        title(tshort{i+1},'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',14); 
        
        print(gcf, '-dtiff', '-r300', [outfig preout 'net_matrix_' tshort{i+1} '.tif']);
    end


end




end % end main function



function [inInd, outInd] = sepROIByMask(roilist, taskROI) % roilist, each row is a ROI
% roilist [N, 3]
% taskROI is in MNI activation ROI mask
    inInd=[];
    outInd=[];
    roiIJK = xyz2ijk_MNI2mm(roilist);
    for i=1:size(roiIJK,1)
        if taskROI(roiIJK(i,1), roiIJK(i,2), roiIJK(i,3)) >0.1
            inInd = [inInd; i];
        else
            outInd = [outInd; i];
        end
    end

end
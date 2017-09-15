%% To examine the network metrics, efficiency, interaction strength, fractional modularity;
%% Includes how to plot the figure 

function network_interact_efficiency_fMod(thr)
% [sec1 sec2]-s to cover the [1-463] subjects

preout='revise_avg_';
%thr=15; 
outdat='/DATA/239/nmzuo/dynNet/code/data_09_8_avg/';
outfig='/datc/dynNet/mod_comp_result/';
load([outdat preout 'totQS_all_bin_thr' num2str(thr) '.mat']); % totQ(numTask,numSubj); totS(numTask,numSubj,nNode);

%addpath(genpath('/DATA/239/nmzuo/dynNet/code/NCT_Bassett')); %% community measurement
addpath(genpath('/DATA/239/nmzuo/dynNet/code/GenLouvain2.0')); %% Mucha methods
addpath(genpath('/DATA/239/nmzuo/software/NIfTI_20140122')); %% load_nii
addpath(genpath('/DATA/239/nmzuo/software/brant/ccm')); %% compute clustering coeficient

%roi_rm = load('roi264_Power_uncertain.txt'); % 28 uncertain ROIs
%roiInd = setdiff([1:264], roi_rm);
roiInd=[1:264];
nSubj = size(totS,2); %463 for LR; 465 for RL
nNode = size(totS,3); %264
nROI = size(totS,3); %264
% 

%load([outdat preout 'corrR_all_remMActivation.mat']); % corrR_all(8,nSubj,nROI,nROI)

% sex, age
%sex_age = load('hcp_S465_sex_age_RL.txt');
%sex_age = load('hcp_S463_sex_age.txt');
sex_age = load('hcp_S453_sex_age_avg.txt');

%tname = {'rfMRI_REST1_RL','tfMRI_GAMBLING_RL','tfMRI_MOTOR_RL','tfMRI_SOCIAL_RL',  ...
         'tfMRI_EMOTION_RL',  'tfMRI_LANGUAGE_RL',  'tfMRI_RELATIONAL_RL', 'tfMRI_WM_RL'};
%tname = {'rfMRI_REST1_LR','tfMRI_GAMBLING_LR','tfMRI_MOTOR_LR','tfMRI_SOCIAL_LR',  ...
%         'tfMRI_EMOTION_LR',  'tfMRI_LANGUAGE_LR',  'tfMRI_RELATIONAL_LR', 'tfMRI_WM_LR'};
tname = {'rfMRI_REST1_avg','tfMRI_GAMBLING_avg','tfMRI_MOTOR_avg','tfMRI_SOCIAL_avg',  ...
%         'tfMRI_EMOTION_avg',  'tfMRI_LANGUAGE_avg',  'tfMRI_RELATIONAL_avg', 'tfMRI_WM_avg'};
tshort = {'REST1','GAMBLING','MOTOR','SOCIAL',  ...
         'EMOTION',  'LANGUAGE',  'RELATIONAL', 'WM'};     
pathbase='../groupActi_RL/';
maskname={'/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% gambling
          '/groupmean.gfeat/cope12345_cluster_mask_zstat1.nii.gz', ...    %% motor
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% social
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% emotion
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% language
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% relational
          '/groupmean.gfeat/cope1.feat/cluster_mask_zstat1.nii.gz', ...   %% wm
         };

coord264 = load('roi264_Power.txt');

interActRest = nan(nSubj,7,2);
interActTask = nan(nSubj,7,2);

interModRest = nan(nSubj,1);
interModTask = nan(nSubj,1);

%% The following for-loop is to compute the task-rest task ROI/Rest interactions by 2 methods
%%('Interaction difference between rest and tasks');
% the following is the output of this section
if 0 % switch to run OR not
    disp('********compute the task-rest task ROI/Rest interactions by 2 methods********');
    for i=1:7  % 1:7   %% for i=1:7
        taskROI = load_nii([pathbase tname{i+1} maskname{i}]); % 3 is for motor
        taskROI = taskROI.img;
        [inInd, outInd] = sepROIByMask(coord264(roiInd,:), taskROI);  
        Nin=length(inInd);  Nout=length(outInd);
        fprintf('%s: [inInd=%d, outInd=%d]\n', tname{i+1}, Nin, Nout);

        for j=1:nSubj
            tmp1 = network_thr_bin(squeeze(corrR_all(1,j,:,:)), thr); 
            tmp2 = network_thr_bin(squeeze(corrR_all(i+1,j,:,:)), thr);
            interModRest(j) = comp_interAct_connWeit(tmp1, inInd, outInd);
            interModTask(j) = comp_interAct_connWeit(tmp2, inInd, outInd);
        end
        interActRest(:,i,2)=interModRest;
        interActTask(:,i,2)=interModTask;
        [h,p2,ci,stats] = ttest(interModTask(:), interModRest(:));
        disp(['    p=', num2str(p2), '  R=', num2str(stats.tstat)]);

    end  %% end for i=1:7
end
%save([outdat preout 'interActRestTask_thr' num2str(thr) '.mat'], 'interActRest', 'interActTask');
load([outdat preout 'interActRestTask_thr' num2str(thr) '.mat']); %, 'interActRest', 'interActTask');

if 0 % whether plot the interAct figures
    addpath('aboxplot');
    tmp1=nan(2,nSubj,7); % construct for the aboxplot()
    tmp1(1,:,:)=interActTask(:,:,1);
    tmp1(2,:,:)=interActRest(:,:,1);
    figure('Position',[0,0, 800,600]);
    aboxplot(tmp1, 'labels', tshort(2:end)); ylim([0,0.5]);
    title(['Task vs. Rest in 7 tasks'], 'interpreter', 'none')    
    print(gcf, '-dtiff', '-r300', [outfig preout 'interAct_s500_FDR_ShareMode_thr' num2str(thr) '.tif']);
    tmp1=nan(2,nSubj,7); % construct for the aboxplot()
    tmp1(1,:,:)=interActTask(:,:,2);
    tmp1(2,:,:)=interActRest(:,:,2);
    figure('Position',[0,0, 800,600]);
    aboxplot(tmp1, 'labels', tshort(2:end)); ylim([0,0.5]);
    title(['Task vs. Rest in 7 tasks'], 'interpreter', 'none')    
    print(gcf, '-dtiff', '-r300', [outfig preout 'interAct_s500_FDR_ConnWeit_thr' num2str(thr) '.tif']);

end % if compute interAct...

load([outdat preout 'in_out_Q_thr' num2str(thr)  '.mat']);
load([outdat preout 'in_out_E_thr' num2str(thr)  '.mat']);   

%% compute in/out modularity and efficiency
recSize=zeros(7,2);
if 0 %switch to run OR not
    disp('*****************  compute in/out modularity and efficiency');
%     QinRest = nan(7,nSubj); QinTask = nan(7,nSubj);
%     QoutRest = nan(7,nSubj); QoutTask = nan(7,nSubj);
%     EGinTask = nan(7,nSubj); EGinRest = nan(7,nSubj);
%     EGoutTask = nan(7,nSubj); EGoutRest = nan(7,nSubj);
%     ELinTask = cell(7,nSubj); ELinRest = cell(7,nSubj);
%     ELoutTask = cell(7,nSubj); ELoutRest = cell(7,nSubj);
%     EG = nan(8,nSubj);
%     EL = nan(8,nSubj,nNode);
  
    
%    matlabpool('open',7);
%    if matlabpool('size')==7
%       disp('********  7 cores applied!  *******'); 
%    end

    for i=1:7   % 1:7   %% for i=1:7
        disp([tname{i+1}]);
        taskROI = load_nii([pathbase tname{i+1} maskname{i}]); % 3 is for motor
        taskROI = taskROI.img;
        [inInd, outInd] = sepROIByMask(coord264(roiInd,:), taskROI);
        recSize(i,:)=[length(inInd), length(outInd)];
%       matlabpool('open',6);
%        if matlabpool('size')==6
%           disp('********  6 cores applied!  *******'); 
%        end
        for j=[] %sec1:sec2  %1:nSubj
%            if mod(j,5) == 0
                fprintf('%d:%d  ', j,i);
%            end
            % inside the ROI
            inTaskMat = squeeze(corrR_all(i+1, j,inInd, inInd));
            inRestMat = squeeze(corrR_all(1, j,inInd, inInd));

            % outside the ROI
            outTaskMat = squeeze(corrR_all(i+1, j,outInd, outInd));
            outRestMat = squeeze(corrR_all(1, j, outInd, outInd));
    %       Thr and bin the matrix (Reviewer 3)
            inTaskMat=network_thr_bin(inTaskMat, thr);
            inRestMat=network_thr_bin(inRestMat, thr);
            outTaskMat=network_thr_bin(outTaskMat, thr);
            outRestMat=network_thr_bin(outRestMat, thr);
            %% Compute the modularity
            if 0
                [Q,S] = Mucha_2D(inTaskMat, 100); QinTask(i,j)=Q; 
                [Q,S] = Mucha_2D(inRestMat, 100); QinRest(i,j)=Q;
                [Q,S] = Mucha_2D(outTaskMat, 100); QoutTask(i,j)=Q; 
                [Q,S] = Mucha_2D(outRestMat, 100); QoutRest(i,j)=Q;
            end

            if 0  %% Compute global efficiency: inside
                [EGinRest(i,j), ELinRest{i,j}] = CCM_GEfficiency(1./inRestMat);
                [EGinTask(i,j), ELinTask{i,j}] = CCM_GEfficiency(1./inTaskMat);
            end

            if 0  %% Compute global efficiency: outside
                [EGoutRest(i,j), ELoutRest{i,j}] = CCM_GEfficiency(1./outRestMat);
                [EGoutTask(i,j), ELoutTask{i,j}] = CCM_GEfficiency(1./outTaskMat);
            end
            if 1  %% Compute global efficiency: full matrix
                if i==1
                    restCorr = network_thr_bin(squeeze(corrR_all(1,j,:,:)), thr);
                    [EG(1,j), EL(1,j,:)] = CCM_GEfficiency(1./restCorr);
                end
                taskCorr = network_thr_bin(squeeze(corrR_all(i+1,j,:,:)), thr);
                [EG(i+1,j), EL(i+1,j,:)] = CCM_GEfficiency(1./taskCorr);
            end

        end % for j=1:nSubj
%        matlabpool('close');
    end
%    matlabpool('close');
    % rest is separate
%     for j=1:nSubj
%         restCorr = network_thr_bin(squeeze(corrR_all(1,j,:,:)), thr);
%         [EG(1,j), EL(1,j,:)] = CCM_GEfficiency(1./restCorr);
%     end
end % whether compute in/out modularity and efficiency
save([outdat preout 'recSize.mat'], 'recSize');
save([outdat preout 'in_out_Q_thr' num2str(thr) '_' num2str(sec1) '_' num2str(sec2) '.mat'], 'QinRest', 'QinTask', 'QoutRest', 'QoutTask');
save([outdat preout 'in_out_E_thr' num2str(thr) '_' num2str(sec1) '_' num2str(sec2) '.mat'], 'EGinTask', 'EGinRest', 'EGoutTask', 'EGoutRest', ...
      'ELinTask', 'ELinRest', 'ELoutTask', 'ELoutRest', 'EG', 'EL');

%% To judge the change of entire Q is relavant ot the interaction? Qin? Qout?
disp('*******To judge the change of entire Q is relavant ot the interaction? Qin? Qout? *******');
%load([outdat preout 'recSize.mat']); %% recSize:7*2
coefQ=nan(7,2);
for i=[]   % i=1:7
    % change of entire Q
    changeQ=squeeze(totQ(i+1,:)-totQ(1,:));
    changeAct=squeeze(interActTask(:,i,1)-interActRest(:,i,1));   
    changeAct=regress_out(changeAct, sex_age(:,2:3));
    [rho, pval]=corr(changeQ',changeAct);
    disp([tname{i+1} ': changeQ vs. Act(ShareMode), ' num2str(pval) '(R=' num2str(rho) ')']);
    
    changeAct2=squeeze(interActTask(:,i,2)-interActRest(:,i,2));
    changeAct2=regress_out(changeAct2, sex_age(:,2:3));
    [rho, pval]=corr(changeQ',changeAct2);
    disp([tname{i+1} ': changeQ vs. Act(connWeit), ' num2str(pval) '(R=' num2str(rho) ')']);
    % Qin
    changeQin = squeeze(QinTask(i,:) - QinRest(i,:));
    changeQin=regress_out(changeQin, sex_age(:,2:3));
    [rho, pval]=corr(changeQ',changeQin');
    disp([tname{i+1} ': changeQ vs. Qin, ' num2str(pval)  '(R=' num2str(rho) ')']);  
    % Qout
    changeQout = squeeze(QoutTask(i,:) - QoutRest(i,:));
    changeQout=regress_out(changeQout, sex_age(:,2:3));
    [rho, pval]=corr(changeQ',changeQout');
    disp([tname{i+1} ': changeQ vs. Qout, ' num2str(pval)  '(R=' num2str(rho) ')']);
    [r, b, s] = regress_out(changeQ', [changeQin'/recSize(i,1), changeQout'/recSize(i,2)]);
    coefQ(i,:) = b(2:3);
%   % similar to the following
%     tt = regstats(changeQ', [changeQin', changeQout'] ,'linear','tstat');
%     disp([tt.tstat.beta, tt.tstat.t, tt.tstat.pval]);
    disp(['beta, coef_Qin, coef_Qout: ' num2str([b' s(3)])]);  % s=[R^2, F, p, error];
        
    if 0 %% if show and save picture
        tifpre='_s500_FDR_remSexAge_';
        figure; plot_regress(changeQ, changeAct); title([tname{i+1} ': ShareMode'], 'interpreter', 'none');
        print(gcf, '-dtiff', '-r300', [outfig preout 'plot_changeQ_ShareMode' tifpre tname{i+1}  '.tif']);
        figure; plot_regress(changeQ, changeAct2); title([tname{i+1} ': connWeit'], 'interpreter', 'none');
        print(gcf, '-dtiff', '-r300', [outfig preout 'plot_changeQ_connWeit' tifpre tname{i+1}  '.tif']);
        figure; plot_regress(changeQ, changeQin); title([tname{i+1} ': Qin'], 'interpreter', 'none');
        print(gcf, '-dtiff', '-r300', [outfig preout 'plot_changeQ_Qin' tifpre tname{i+1}  '.tif']);
        figure; plot_regress(changeQ, changeQout); title([tname{i+1} ': Qout'], 'interpreter', 'none');
        print(gcf, '-dtiff', '-r300', [outfig preout 'plot_changeQ_Qout' tifpre tname{i+1}  '.tif']);
    end
end


%% To judge the change of entire Q is relavant ot Clustering Coeff (Cin? Cout?)?
% change of Q vs. interaction
for i=[] %% i=1:7
    disp('******To judge the change of entire Q is relavant ot Clustering Coeff (Cin? Cout?)? ******');
    % change of entire Q
    changeQ=squeeze(totQ(i+1,:)-totQ(1,:));
    % Cin
    changeCin = squeeze(Cin(i+1,:) - Cin(1,:));
    [rho, pval]=corr(changeQ',changeCin');
    disp([tname{i+1} ': changeQ vs. Cin, ' num2str(pval)  '(R=' num2str(rho) ')']);  
    % Cout
    changeCout = squeeze(Cout(i+1,:) - Cout(1,:));
    [rho, pval]=corr(changeQ',changeCout');
    disp([tname{i+1} ': changeQ vs. Cout, ' num2str(pval)  '(R=' num2str(rho) ')']);
        
    figure; plot_regress(changeQ, changeCin); title([tname{i+1} ': Cin'], 'interpreter', 'none');
    print(gcf, '-dtiff', '-r300', [outpath 'plot_changeQ_Cin_s500_FDR_' tname{i+1}  '.tif']);
    figure; plot_regress(changeQ, changeCout); title([tname{i+1} ': Cout'], 'interpreter', 'none');
    print(gcf, '-dtiff', '-r300', [outpath 'plot_changeQ_Cout_s500_FDR_' tname{i+1}  '.tif']);
end


%% participation coefficient, to find high participation nodes
if 0
    addpath(genpath('/datc/software/BCT/BCT_20150125')); % using participation_coef 
    [M,N,K]= size(totS);
    pCoef = zeros([M,N,K]);
    for i=1:M
        for j=1:N
            tmp =squeeze(corrR_all(i,j,:,:));
            tmp(isnan(tmp))=0;
            [ppos, pneg] = participation_coef_sign(tmp, squeeze(totS(i,j,:)));
            pCoef(i,j,:)=reshape(ppos+pneg, [K,1]);
        end
    end
    %% for each task, test the most participation_coef 
    roi_rm = load('roi264_Power_uncertain.txt'); % 28 uncertain ROIs
    roiInd = setdiff([1:264], roi_rm);
    roi_org = load('roi264_Power.txt'); % 264 ROIs
    refnii='/datc/software/fsl5.0/data/standard/MNI152_T1_2mm_brain.nii.gz';
    for i=1:7
        thisCoef=squeeze(pCoef(i+1,:,:)-pCoef(1,:,:));
        thisCoef = thisCoef - repmat(mean(thisCoef,2),[1,K]);
        [hval, pval, ci, stats] = ttest(thisCoef);
        pvalC = mafdr(pval, 'BHFDR', 'true');
        tstat = zeros(length(pval),1);
        tstat(pvalC<0.05) = stats.tstat(pvalC<0.05);
        xyz=roi_org(roiInd,:);
        ijk=xyz2ijk_MNI2mm(xyz);
        % inT
        prefile=['./data/detect_PartCoef_TpR_' tname{i+1} '.nii.gz'];
        roi2nii_MNI2mm([ijk,tstat],refnii,prefile);
    end

end %% compute participation coefficient

end % main function


function [r, b, stats] = regress_out(Y, counf)
% Y =[n,1], counf=[n,p], there are p counfound variables
    [M,N] = size(Y);
    Y = reshape(Y, [length(Y),1]);
    X = [ones(size(counf,1),1), counf];
    [b,bint,r,rint,stats] = regress(Y, X);
    r = reshape(r, [M,N]);
end

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


% M is node-connection strength matrix
function nCount = comp_interAct_connBin(M, aInd, bInd)
    nCount = 0;
%     for i=1:length(aInd)
%         for j=1:length(bInd)
%             if abs(M(aInd(i),bInd(j))) > 1.0e-6
%                 nCount = nCount+1;
%             end
%         end
%     end
    M(logical(diag(ones(length(M),1))))=0;
    M=logical(M);
    totConn=sum(M(:))/2.0;
    M(aInd,aInd)=0;
    M(bInd,bInd)=0;
    nCount=sum(M(:))/2.0/totConn; %% normalize the conn.
end
% M is node-connection strength matrix
function nCount = comp_interAct_connWeit(M, aInd, bInd)
    nCount = 0;
%     for i=1:length(aInd)
%         for j=1:length(bInd)
%             nCount = nCount + M(aInd(i),bInd(j));
%         end
%     end
    M(logical(diag(ones(length(M),1))))=0;
    totConn=sum(M(:))/2.0;
    M(aInd,aInd)=0;
    M(bInd,bInd)=0; 
    nCount=sum(M(:))/2.0/totConn;  %% normalize the conn Weight.
end
% nPart (N*1) is a partition scheme for N nodes
function nCount = comp_interAct_shareMod(nPart, aInd, bInd)
    nCount = 0;
    for i=1:length(aInd)
        for j=1:length(bInd)
            if nPart(aInd(i)) == nPart(bInd(j))
                nCount = nCount + 1;
            end
        end
    end
    nCount=nCount/(length(aInd) * length(bInd));  %% normalize the conn Weight.
end

% M is a node-connection strength
function Q = comp_roiMod(M, aInd)
    aM = M(aInd, aInd);
    Q = Mucha_2D_signed(aM);
end




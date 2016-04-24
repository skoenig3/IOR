% IOR Script analysis.

% SCM 10 secs
% Social SCM 10 secs
% LIST 10 secs
% SMT 7 secs
% SFT 5 secs
% VPLT variable

% [1] Get Salience maps
% [2] Get Eye Data
% [3] Salience at Fixation Locations
% [4] Salience and IOR
% [5] Salience and IOR Repeat Trials
% [6] Run Seperate Special Salience IOR function for SSCM

addpath('C:\Users\seth.koenig\Documents\MATLAB\IOR')
%% [1] Create Salience Maps...if not created already
Matlab_dir = 'C:\Users\seth.koenig\Documents\MATLAB\';

% for SCM
scm_dir = [Matlab_dir 'BCRW Salience Model\SCM Image Sets\'];
scm_sets = {'Set006','Set007','Set008','Set009','SetE001','SetE002','SetE003','SetE004'};
for set = 1:length(scm_sets)
    cd([scm_dir scm_sets{set}]);
    for img = 1:36
        getSalienceMap([num2str(img) '.bmp']);
    end
end

% for LIST
LIST_dir = [Matlab_dir 'List Task Analysis\List Image Sets\'];
for set = 1:12;
    if set < 10
        cd([LIST_dir 'List0' num2str(set)])
    else
        cd([LIST_dir 'List' num2str(set)])
    end
    for img = 1:90
        getSalienceMap([num2str(img) '.bmp']);
    end
end

%for SMT
SMT_dir = [Matlab_dir 'Scene Memory Task\Image Sets\'];
SMT_sets = {'SMT001','SMT002','SMT003','SMT004','SMT005'};
for set = 1:length(SMT_sets)
    cd([scm_dir SMT_sets{set}]);
    for img = 1:25
        getSalienceMap([num2str(img) '.bmp']);
    end
end

%for SFT 90 images
SFT_dir = [Matlab_dir 'Shift Task\'];
SFT_sets = {'sft22','sft24','sft27'};
for set = 1:length(scm_sets)
    cd([scm_dir SFT_sets{set}]);
    for img = 1:90
        getSalienceMap([num2str(img) '.bmp']);
    end
end

%for SFT 120 images
SFT_dir = [Matlab_dir 'Shift Task\'];
SFT_sets = {'sft5','sft6','sft8','sft9'};
for set = 1:length(scm_sets)
    cd([scm_dir SFT_sets{set}]);
    for img = 1:120
        getSalienceMap([num2str(img) '.bmp']);
    end
end

%for VPLT
VPLT_dir = [Matlab_dir 'VPLT\ImageFiles\'];
VPLT_sets = [126:135 142:144];
for set = VPLT_sets;
    cd([VPLT_dir 'SET' num2str(set)]);
    for img = 1:200
        getSalienceMapVPLT([num2str(img) '.bmp'])
    end
end

% for Social SCM
SSCM_dir = [Matlab_dir '\SSCM\'];
for set = 2:7
    cd([SSCM_dir 'S' num2str(set)]);
    l = ls;
    hid = find(l(:,1) == '.');
    l(hid,:) = [];
    for bmp = 1:size(l,1);
        if ~isempty(strfind(l(bmp,:),'.bmp'));
            getSalienceMap(l(bmp,:));
        end
    end
end

%% [2] Get eye Data
% color change calibration for SMT and SFT tasks using 25 point t2form grid calibration
% color change calibration for SCM, VPLT, and SSCM using affine transformation
Matlab_dir = ['C:\Users\seth.koenig\Documents\MATLAB\'];

IOR_dir = [Matlab_dir 'IOR\'];

% for SCM
scm_dir = [Matlab_dir 'BCRW Salience Model\SCM Image Sets\'];
scm_sets = {'Set006','Set007','Set008','Set009','SetE001','SetE002','SetE003','SetE004'};
imageX = 800;
imageY = 600;

mkdir([IOR_dir 'Eye Data\SCM'])
for set = 1:length(scm_sets)
    cd([scm_dir scm_sets{set}]);
    
    dirData = dir([scm_dir scm_sets{set}]);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        [eyedat,notattended] = getSCMeyedat(cortexfile,0.005,imageX,imageY);
        fixationstats = ClusterFixation_Final(eyedat);
        imgset = scm_sets{set};
        save([IOR_dir 'Eye Data\SCM\' cortexfile(1:end-2) '_' cortexfile(end) '-fixation.mat'],...
            'fixationstats','notattended','imgset')
    end
end

%for LIST
mkdir([Matlab_dir 'IOR\Eye Data\LIST'])
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
imageX = 800;
imageY = 600;
CNDFile = [list_image_dir 'List.cnd'];
ITMFile = [list_image_dir 'List01.itm'];
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        getLISTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY,imset)
    end
end

%for SMT
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory
cd(SMT_dir)
load('allcortexfiles.mat'); %easily accessible data structure from above
imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

mkdir([Matlab_dir 'IOR\Eye Data\SMT'])
for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    imset = allcortexfiles{crow,3};
    if imset < 10
        image_dir = [SMT_dir  'Image Sets\SMT00' num2str(imset)];
        CNDFile = [SMT_dir 'ITMCNDfiles\SMT00' num2str(imset) '.CND'];
        ITMFile = [SMT_dir 'ITMCNDfiles\SMT00' num2str(imset) '.ITM'];
    else
        image_dir = [SMT_dir  'Image Sets\SMT0' num2str(imset)];
        CNDFile = [SMT_dir 'ITMCNDfiles\SMT0' num2str(imset) '.CND'];
        ITMFile = [SMT_dir 'ITMCNDfiles\SMT0' num2str(imset) '.ITM'];
    end
    if strcmpi(cortexfile,'PW131024.1')
        CNDFile = [SMT_dir 'ITMCNDfiles\PWSMT00' num2str(imset) '.CND']; %error compiling cortex smt.sav file so had to rewrite for this trial
    end
    getSMTeyedat2([SMT_dir 'cortexfiles\' cortexfile],ITMFile,CNDFile,imageX,imageY)
    % grab eye data, calibrate, and detect fixations and saccades
end

% for SFT WR
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\';
image_sets = {'sft22','sft24','sft27','sft5'};
imageX = 756;
imageY = 378;
CNDFile = [sft_image_dir 'sft_cnd.cnd'];
mkdir([IOR_dir 'Eye Data\SFT'])
for imset = 1:length(image_sets);
    dirName = [sft_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        ITMFile = [image_sets{imset} '.itm'];
        getSFTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY)
    end
end

% for SFT MP
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\';
image_sets = {'sft6','sft8','sft9',};
imageX = 756;
imageY = 378;
CNDFile = [sft_image_dir 'sft_cnd2.cnd'];
for imset = 1:length(image_sets);
    dirName = [sft_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        ITMFile = [image_sets{imset} '.itm'];
        getSFTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY)
    end
end

%for VPLT
VPLTdir = 'C:\Users\seth.koenig\Documents\MATLAB\VPLT\';
CNDFile = [VPLTdir 'VPCbmp.CND'];
mkdir([IOR_dir 'Eye Data\VPLT'])
imageX=264;
imageY=264;

[NUM,TXT,RAW]=xlsread([VPLTdir 'masterlist2.xls']);%masterlist2 for JN
RAW = RAW(1:4,:);% for JN
for r = 1:size(RAW,1)-1;
    ITMFile = [VPLTdir 'ItemFiles\VPC' num2str(NUM(r)) '.ITM'];
    cortexfile = TXT{r+1,1};
    cortexfile = [VPLTdir 'DataFiles\' cortexfile];
    [eyedat,imageprez] = getVPLTeyedat(cortexfile,ITMFile,CNDFile);
    fixationstats = ClusterFixation_Final(eyedat);
    imgset = NUM(r);
    save([IOR_dir 'Eye Data\VPLT\' cortexfile(end-9:end-2) '_' cortexfile(end)],...
        'fixationstats','imageprez','imgset');
end

[NUM,TXT,RAW]=xlsread([VPLTdir 'masterlist.xls']);
for r = 1:size(RAW,1)-1;
    ITMFile = [VPLTdir 'ItemFiles\VPC' num2str(NUM(r)) '.ITM'];
    cortexfile = TXT{r+1,1};
    cortexfile = [VPLTdir 'DataFiles\' cortexfile];
    [eyedat,imageprez] = getVPLTeyedat(cortexfile,ITMFile,CNDFile);
    fixationstats = ClusterFixation_Final(eyedat);
    imgset = NUM(r);
    save([IOR_dir 'Eye Data\VPLT\' cortexfile(end-9:end-2) '_' cortexfile(end)],...
        'fixationstats','imageprez','imgset');
end

% for Social SCM saline only no oxytocin
% Run SSCMcode.m
%% [3] Salience at Fixation Locations
Matlab_dir = 'C:\Users\seth.koenig\Documents\MATLAB\';
IOR_dir = [Matlab_dir 'IOR\'];

tasks = {'SCM','LIST','SMT','VPLT','SSCM'}; %,'SFT'

image_dirs = {
    'BCRW Salience Model\SCM Image Sets\';
    'List Task Analysis\List Image Sets\';
    'Scene Memory Task\Image Sets\';
    'VPLT\ImageFiles\';
    'SSCM\'
    };     %'Shift Task\';

nov_salience = cell(1,length(image_dirs));
rep_salience = cell(1,length(image_dirs));
rand_salience = cell(1,length(image_dirs));

fixdurlow = 0; fixdurhigh = 5000; %full spectrum

nov_count = ones(1,length(image_dirs));
rep_count = ones(1,length(image_dirs));
for task = 1:length(image_dirs);
    nov_salience{task} = NaN(1750,35);
    rep_salience{task} = NaN(1750,35);
    rand_salience{task} = NaN(1750,35);
    attention{task} = [];
    
    cd([IOR_dir 'Eye Data\' tasks{task}]);
    a = what;
    for m = 1:length(a.mat)
        if isempty(strfind(a.mat{m},'IOR'))
            load(a.mat{m});
            
            if strcmpi(tasks{task},'SCM')||strcmpi(tasks{task},'SSCM')
                novelconditions = 1:2:length(fixationstats);
                repeatconditions = 2:2:length(fixationstats);
            elseif strcmpi(tasks{task},'LIST')
                novelconditions = [];
                for im = 1:length(images)
                    ind = find(images == images(im));
                    if length(ind) == 1 || ind(1) == im
                        novelconditions = [novelconditions im];
                    end
                end
                repeatconditions = 1:length(images);
                [~,ia,~] = intersect(repeatconditions,novelconditions);
                repeatconditions(ia) = [];
                %         elseif strcmpi(tasks{task},'SFT')
                %             novelconditions = 1:length(fixationstats);
            elseif strcmpi(tasks{task},'SMT')
                novelconditions = [];
                for im = 1:ceil(length(images)/2);
                    ind = find(images == im);
                    if ~isempty(ind)
                        novelconditions = [novelconditions ind(1)];
                    end
                end
                repeatconditions = 1:length(images);
                [~,ia,~] = intersect(repeatconditions,novelconditions);
                repeatconditions(ia) = [];
            elseif strcmpi(tasks{task},'VPLT')
                novelconditions = 1:length(fixationstats)/2;
                repeatconditions = novelconditions(end)+1:length(fixationstats);
            end
            
            for cnd = 1:length(novelconditions);
                if strcmpi(tasks{task},'VPLT')
                    load([Matlab_dir image_dirs{task} 'SET' num2str(imgset) '\' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    fixationdata = fixationstats{imageprez(1,cnd)};
                else
                    fixationdata = fixationstats{novelconditions(cnd)};
                end
                
                if strcmpi(tasks{task},'SCM')
                    load([Matlab_dir image_dirs{task} imgset '\' num2str(cnd) '-saliencemap.mat'],'fullmap')
                elseif strcmpi(tasks{task},'LIST')
                    imgsetdir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
                    if imgset < 10
                        load([imgsetdir 'List0' num2str(imgset) '\' num2str(images(novelconditions(cnd))) '-saliencemap.mat'],'fullmap')
                    else
                        load([imgsetdir 'List' num2str(imgset) '\' num2str(images(novelconditions(cnd))) '-saliencemap.mat'],'fullmap')
                    end
                    %             elseif strcmpi(tasks{task},'SFT')
                    %                 if strcmpi(imgset,'sft5') || strcmpi(imgset,'sft22')|| strcmpi(imgset,'sft24') || strcmpi(imgset,'sft27')
                    %                     if cnd < 10
                    %                         load([Matlab_dir image_dirs{task} imgset '\0' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    %                     else
                    %                         load([Matlab_dir image_dirs{task} imgset '\' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    %                     end
                    %                 else
                    %                     if cnd < 10
                    %                         load([Matlab_dir image_dirs{task} imgset '\00' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    %                     elseif cnd < 100
                    %                         load([Matlab_dir image_dirs{task} imgset '\0' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    %                     else
                    %                         load([Matlab_dir image_dirs{task} imgset '\' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    %                     end
                    %                 end
                elseif strcmpi(tasks{task},'SMT')
                    load([Matlab_dir image_dirs{task} imgset(end-5:end) '\' num2str(images(novelconditions(cnd))) '-saliencemap.mat'],'fullmap')
                elseif strcmpi(tasks{task},'SSCM')
                    load([Matlab_dir image_dirs{task} 'S' num2str(imgset) '\' images{novelconditions(cnd)} '-saliencemap.mat'],'fullmap')
                end
                
                if ~isempty(fixationdata.fixationtimes)
                    
                    %                 if strcmpi(tasks{task},'SFT')
                    %                     fixationdata.XY(2,:) = size(fullmap,1)-fixationdata.XY(2,:);
                    %                     crosspos = initialcross(cnd,:);
                    %                     if fixationdata.fixations(1,1) > crosspos(1)-100 && fixationdata.fixations(1,1) < crosspos(1)+100 &&...
                    %                             fixationdata.fixations(2,1) < crosspos(2)+100 && fixationdata.fixations(2,1) > crosspos(2)-100
                    %                         fixationdata.fixations(:,1) = [];
                    %                         fixationdata.fixationtimes(:,1) = [];
                    %                     end
                    %                 else
                    imageX = size(fullmap,2);
                    imageY = size(fullmap,1);
                    if fixationdata.fixations(1,1) > imageX/2-100 && fixationdata.fixations(1,1) < imageX/2+100 &&...
                            fixationdata.fixations(2,1) < imageY/2+100 && fixationdata.fixations(2,1) > imageY/2-100
                        fixationdata.fixations(:,1) = [];
                        fixationdata.fixationtimes(:,1) = [];
                    end
                    %                 end
                    
                    if ~isempty(fixationdata.fixations)
                        data = Salience_at_fixations(fixationdata,fullmap,fixdurlow,fixdurhigh);
                        
                        if length(data.salience) > 35
                            nov_salience{task}(nov_count(task),:) = data.salience(1:35);
                            randnov_salience{task}(nov_count(task),:) = data.shuffled(1:35);
                        else
                            nov_salience{task}(nov_count(task),1:length(data.salience)) = data.salience;
                            rand_salience{task}(nov_count(task),1:length(data.salience)) = data.shuffled;
                        end
                    end
                    nov_count(task) = nov_count(task)+1;
                end
            end
            if isempty(strfind(a.mat{m},'TT'))
                for cnd = 1:length(repeatconditions);
                    if strcmpi(tasks{task},'VPLT')
                        load([Matlab_dir image_dirs{task} 'SET' num2str(imgset) '\' num2str(cnd) '-saliencemap.mat'],'fullmap')
                        fixationdata = fixationstats{imageprez(2,cnd)};
                    else
                        fixationdata = fixationstats{repeatconditions(cnd)};
                    end
                    
                    if strcmpi(tasks{task},'SCM')
                        load([Matlab_dir image_dirs{task} imgset '\' num2str(cnd) '-saliencemap.mat'],'fullmap')
                    elseif strcmpi(tasks{task},'LIST')
                        imgsetdir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
                        if imgset < 10
                            load([imgsetdir 'List0' num2str(imgset) '\' num2str(images(repeatconditions(cnd))) '-saliencemap.mat'],'fullmap')
                        else
                            load([imgsetdir 'List' num2str(imgset) '\' num2str(images(repeatconditions(cnd))) '-saliencemap.mat'],'fullmap')
                        end
                    elseif strcmpi(tasks{task},'SMT')
                        load([Matlab_dir image_dirs{task} imgset(end-5:end) '\' num2str(images(repeatconditions(cnd))) '-saliencemap.mat'],'fullmap')
                    elseif strcmpi(tasks{task},'SSCM')
                        load([Matlab_dir image_dirs{task} 'S' num2str(imgset) '\' images{repeatconditions(cnd)} '-saliencemap.mat'],'fullmap')
                    end
                    
                    if ~isempty(fixationdata.fixationtimes)
                        
                        imageX = size(fullmap,2);
                        imageY = size(fullmap,1);
                        if fixationdata.fixations(1,1) > imageX/2-100 && fixationdata.fixations(1,1) < imageX/2+100 &&...
                                fixationdata.fixations(2,1) < imageY/2+100 && fixationdata.fixations(2,1) > imageY/2-100
                            fixationdata.fixations(:,1) = [];
                            fixationdata.fixationtimes(:,1) = [];
                        end
                        
                        if ~isempty(fixationdata.fixations)
                            data = Salience_at_fixations(fixationdata,fullmap,fixdurlow,fixdurhigh);
                            
                            if length(data.salience) > 35
                                salience{task}(rep_count(task),:) = data.salience(1:35);
                            else
                                rep_salience{task}(rep_count(task),1:length(data.salience)) = data.salience;
                            end
                        end
                        rep_count(task) = rep_count(task)+1;
                    end
                end
            end
        end
    end
end

numpoints = zeros(2,length(tasks));
for task = 1:length(tasks)
    temp = sum(~isnan(nov_salience{task}),2);
    temp(temp == 0) = [];
    numpoints(1,task) = median(temp);
    temp = sum(~isnan(rep_salience{task}),2);
    temp(temp == 0) = [];
    numpoints(2,task) = median(temp);
end
numpoints(numpoints < 10) = 10;

randsalcI = NaN(1,length(tasks));
for task = 1:length(tasks)
    rs = rand_salience{task}(1:end);
    [~,p,ci] = ztest(rs,nanmean(nov_salience{task}(:,1)),nanstd(rs),...
        0.05);
    randsalcI(task) = ci(2);
end

nov_vs_rep_sal = [];
for task = 1:length(tasks)
    [~,p] = ttest2(nov_salience{task}(1:end),rep_salience{task}(1:end));
    nov_vs_rep_sal(task) = p;
end
%%
splots = [1 2 0 3 4];
figure
for task = 1:length(tasks)
    if task ~= 3
        subplot(2,2,splots(task))
        hold on
        plot(nanmean(nov_salience{task}(:,1:numpoints(1,task))),'b')
        plot(nanmean(rep_salience{task}(:,1:numpoints(2,task))),'r')
        plot(0:numpoints(1,task),randsalcI(task)*ones(1,numpoints(1,task)+1),'k--')
        hold off
        title(tasks{task})
        legend('Novel','Repeat','Chance','Location','Best')
        xlabel('Fixation Number')
        ylabel('Salience')
        ylim([0.24 0.45])
        xlim([0 max(numpoints(:,task))])
    end
end
%% Test if significant difference in salience at fixation locations
salience_vals =[];
taskgroup = [];
for i = 1:length(tasks);
    salience_vals = [salience_vals;nov_salience{i}(:,1)];
    taskgroup =[taskgroup; i*ones(length(nov_salience{task}),1)];
end
[P,ANOVATAB,STATS] = anova1(salience_vals,taskgroup);
multcompare(STATS)
%%
for task = 1:length(tasks)
    if task ~= 4%vplt not exponential decay
        sal = nanmean(nov_salience{task}(:,1:numpoints(1,task)));
        figure
        subplot(1,2,1)
        plot(sal)
        hold on
        sal2 = filtfilt(1/5*ones(1,5),1,sal);
        plot(sal2,'r')
        ds = abs(diff(sal2));
        [~,ds] = min(ds);
        ds = ds-1;
        if ds > 18 %17 produces the same result
            ds = 18;
        end
        sal = sal(1:ds);
        sal = sal - min(sal);
        zeroind = find(sal == 0);
        sal = sal(1:zeroind-1);
        p = polyfit(1:length(sal),log(sal),1);
        subplot(1,2,2)
        plot(sal)
        hold on
        plot(1:length(sal),exp(p(2))*exp(p(1)*[1:length(sal)]),'r')
        subtitle([tasks{task} ' with ~ tau_{IOR} = 1/' num2str(-1/(p(1)/3.71))]);
    end
end
%% [4] Salience and IOR %%
% determine number of fixations between prior and return
Matlab_dir = ['C:\Users\seth.koenig\Documents\MATLAB\'];
IOR_dir = [Matlab_dir 'IOR\'];
cd(IOR_dir)

tasks = {'SCM','LIST','VPLT','SSCM'};

image_dirs = {
    'BCRW Salience Model\SCM Image Sets\';
    'List Task Analysis\List Image Sets\';
    'VPLT\ImageFiles\';
    'SSCM\'
    };

% img_sets = {
%     {'Set006','Set007','Set008','Set009','SetE001','SetE002','SetE003','SetE004'};
%     {'LIST01','LIST02','LIST03','LIST04','LIST05','LIST06','LIST07','LIST08','LIST09','LIST10','LIST11','LIST12'};
%     {'SET126','SET126','SET127','SET128','SET129','SET130','SET131','SET132','SET133','SET134','SET135','SET142','SET143','SET144'};
%     {'S2','S3','S4','S5','S6','S7'};
%     };

imageX =[800 800 264 800];
imageY =[600 600 264 600];

distancethreshold = [0  24 48 72 96  120  144 168 200 400];%in pixels 24 pixels/dva
distancethreshold = [distancethreshold;
    [24 48 72 96 120 144  168 200 400 800]];%in pixels 24 pixels/dva

for task = 1:length(tasks)
    cd([IOR_dir 'Eye Data\' tasks{task}])
    a = what;
    a = a.mat;
    for aa = 1:size(a,1);
        load(a{aa})
        if length(a{aa}) == 14 || ~isempty(strfind(a{aa},'fixation'))
            if task == 1
                set_dir = [Matlab_dir image_dirs{task} imgset '\'];
            elseif task == 2
                if imgset < 10
                    set_dir = [Matlab_dir image_dirs{task} 'List0' num2str(imgset) '\'];
                else
                    set_dir = [Matlab_dir image_dirs{task} 'List' num2str(imgset) '\'];
                end
            elseif task == 3
                set_dir = [Matlab_dir image_dirs{task} 'SET' num2str(imgset) '\'];
            elseif task == 4
                set_dir = [Matlab_dir image_dirs{task} 'S' num2str(imgset) '\'];
            end
            
            if task == 1
                images = 1:36;
                novelconditions = 1:2:72;
                min_dist = 10;
            elseif task == 2
                imgs = [];
                novelconditions = [];
                for img = 1:max(images);
                    ind = find(images == img);
                    if ~isempty(ind)
                        novelconditions =[novelconditions ind(1)];
                        imgs = [imgs img];
                    end
                end
                images = imgs;
                min_dist = 10;
            elseif task == 3
                novelconditions = imageprez(1,:);
                images = 1:200;
                min_dist = 5;
            elseif task == 4
                novelconditions = 1:2:length(images);
                images = images(novelconditions);
                min_dist = 10;
            end
            disp(['Running Task ' tasks{task} ' fixiation file ' num2str(aa)])
            SalienceIOR(a{aa},distancethreshold,images,set_dir,imageX(task),imageY(task),novelconditions,min_dist)
        end
    end
end

tasks = {'SCM','LIST','VPLT','SSCM'};

task_dirs = {
    'BCRW Salience Model\SCM Image Sets\';
    'List Task Analysis\List Image Sets\';
    'VPLT\ImageFiles\';
    'SSCM\'
    };

Matlab_dir = ['C:\Users\seth.koenig\Documents\MATLAB\'];
IOR_dir = [Matlab_dir 'IOR\Eye Data\'];
cd(IOR_dir)

tasks = {'SCM','LIST','VPLT','SSCM'};

distancethreshold = [0  24 48 72 96  120  144 168 200 400];%in pixels 24 pixels/dva
distancethreshold = [distancethreshold;
    [24 48 72 96 120 144  168 200 400 800]];%in pixels 24 pixels/dva

totalfixations = zeros(length(tasks),size(distancethreshold,2));
allSalIOR = cell(length(tasks),size(distancethreshold,2));
for task = 1:length(tasks);
    cd([IOR_dir tasks{task}])
    matfiles = what;
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'SalienceIOR'))
            load(matfiles.mat{i});
            for ii = 1:size(allSalIOR,2);
                allSalIOR{task,ii} = [allSalIOR{task,ii}; returnfixsal{ii}];
                totalfixations(task,ii) = totalfixations(task,ii)+size(returnfixsal{ii},1);
            end
        end
    end
end

for task = 1:length(tasks);
    labels = {};
    for d = 1:size(allSalIOR,2)
        return_sal_means(2,d) = nanmean(allSalIOR{task,d}(:,8));
        return_sal_stds(2,d) = nanstd(allSalIOR{task,d}(:,8))/sqrt(sum(~isnan(allSalIOR{task,d}(:,8))));
        return_sal_means(1,d) = nanmean(allSalIOR{task,d}(:,4));
        return_sal_stds(1,d) = nanstd(allSalIOR{task,d}(:,4))/sqrt(sum(~isnan(allSalIOR{task,d}(:,4))));
        labels{d} = [num2str(distancethreshold(1,d)) '-' num2str(distancethreshold(2,d))];
    end
    
    figure
    hold on
    errorbar(return_sal_means(2,:),return_sal_stds(2,:))
    errorbar(return_sal_means(1,:),return_sal_stds(1,:),'g')
    title(tasks{task})
    xlabel('Distance Tresholds')
    ylabel('Salience')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    hold off
    ylim([0.25 0.45])
    legend('Salience @ Return Fixation','Salience @ Prior Fixation')
end

for task = 1:length(tasks);
    labels = {};
    for d = 1:size(allSalIOR,2)
        return_dur_means(1,d) = nanmean(5*allSalIOR{task,d}(:,10));
        return_dur_stds(1,d) = nanstd(5*allSalIOR{task,d}(:,10))/sqrt(sum(~isnan(allSalIOR{task,d}(:,10))));
        return_dur_means(2,d) = nanmean(5*allSalIOR{task,d}(:,11));
        return_dur_stds(2,d) = nanstd(5*allSalIOR{task,d}(:,11))/sqrt(sum(~isnan(allSalIOR{task,d}(:,11))));
        labels{d} = [num2str(distancethreshold(1,d)) '-' num2str(distancethreshold(2,d))];
    end
    
    figure
    hold on
    errorbar(return_dur_means(2,:),return_dur_stds(2,:))
    errorbar(return_dur_means(1,:),return_dur_stds(1,:),'g')
    title(tasks{task})
    xlabel('Distance Tresholds')
    ylabel('Fixation Duration (ms)')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    hold off
    legend('Fix dur @ Return Fixation','Fix dur @ Prior Fixation')
end

rs = [];
ps = [];
slopes = [];
for task = 1:length(tasks);
    dist = [];
    sal = [];
    for i = 1:2
        for ii = 1:size(allSalIOR{task,i},1)
            dist = [dist ; allSalIOR{task,i}(ii,9)];
            sal = [sal; allSalIOR{task,i}(ii,8)];
        end
    end
    
    P = polyfit(dist,sal,1);
    figure
    hold on
    plot(dist,sal,'.')
    plot(1:50,[1:50]*P(1)+P(2),'r')
    hold off
    title(tasks{task})
    xlabel('Distance (pixels)')
    ylabel('Normalized Salience')
    [r,p]= corrcoef(dist,sal);
    
    slopes = [slopes P(1)];
    rs = [rs r(2)];
    ps = [ps p(2)];
end

for task = 1:length(tasks);
    labels = {};
    for d = 1:size(allSalIOR,2)
        mean_time_btwn_return(d) = 5*mean(allSalIOR{task,d}(:,7)-allSalIOR{task,d}(:,3)-1);
        mean_numfixations_btwn_return(d) = mean(allSalIOR{task,d}(:,13)-allSalIOR{task,d}(:,12)-1);
        
        std_time_btwn_return(d) = 5*std(allSalIOR{task,d}(:,7)-allSalIOR{task,d}(:,13)-1);
        std_numfixations_btwn_return(d) = std(allSalIOR{task,d}(:,13)-allSalIOR{task,d}(:,12)-1);
        
        numreturns(d) = sum(~isnan(allSalIOR{task,d}(:,13)));
        
        labels{d} = [num2str(distancethreshold(1,d)) '-' num2str(distancethreshold(2,d))];
    end
    
    figure
    subplot(1,2,1)
    errorbar(mean_time_btwn_return,std_time_btwn_return./sqrt(numreturns));
    xlabel('Distance Tresholds')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    ylabel('Duration between prior and return fixation')
    
    subplot(1,2,2)
    errorbar(mean_numfixations_btwn_return,std_numfixations_btwn_return./sqrt(numreturns))
    xlabel('Distance Tresholds')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    ylabel('Number of fixations in-between prior and return fixation')
    subtitle(tasks{task})
end
%% [5] Salience and IOR Repeat Trials
% determine number of fixations between prior and return
Matlab_dir = ['C:\Users\seth.koenig\Documents\MATLAB\'];
IOR_dir = [Matlab_dir 'IOR\'];
cd(IOR_dir)

tasks = {'SCM','LIST','VPLT','SSCM'};

image_dirs = {
    'BCRW Salience Model\SCM Image Sets\';
    'List Task Analysis\List Image Sets\';
    'VPLT\ImageFiles\';
    'SSCM\'
    };

imageX =[800 800 264 800];
imageY =[600 600 264 600];

distancethreshold = [0  24 48 72 96  120  144 168 200 400];%in pixels 24 pixels/dva
distancethreshold = [distancethreshold;
    [24 48 72 96 120 144  168 200 400 800]];%in pixels 24 pixels/dva

for task = 1:length(tasks)
    cd([IOR_dir 'Eye Data\' tasks{task}])
    a = what;
    a = a.mat;
    for aa = 1:size(a,1);
        load(a{aa})
        if length(a{aa}) == 14 || ~isempty(strfind(a{aa},'fixation'))
            if task == 1
                set_dir = [Matlab_dir image_dirs{task} imgset '\'];
            elseif task == 2
                if imgset < 10
                    set_dir = [Matlab_dir image_dirs{task} 'List0' num2str(imgset) '\'];
                else
                    set_dir = [Matlab_dir image_dirs{task} 'List' num2str(imgset) '\'];
                end
            elseif task == 3
                set_dir = [Matlab_dir image_dirs{task} 'SET' num2str(imgset) '\'];
            elseif task == 4
                set_dir = [Matlab_dir image_dirs{task} 'S' num2str(imgset) '\'];
            end
            
            if task == 1
                images = 1:floor(length(fixationstats)/2);
                repeatconditions = 2:2:length(fixationstats);
                min_dist = 10;
            elseif task == 2
                imgs = [];
                repeatconditions = [];
                for img = 1:max(images);
                    ind = find(images == img);
                    if length(ind) == 2
                        repeatconditions =[repeatconditions ind(2)];
                        imgs = [imgs img];
                    end
                end
                images = imgs;
                min_dist = 10;
            elseif task == 3
                repeatconditions = imageprez(2,:);
                images = 1:200;
                min_dist = 5;
            elseif task == 4
                repeatconditions = 2:2:length(images);
                images = images(repeatconditions);
                min_dist = 10;
            end
            disp(['Running Task ' tasks{task} ' fixiation file ' num2str(aa)])
            SalienceIOR(a{aa},distancethreshold,images,set_dir,imageX(task),imageY(task),repeatconditions,min_dist)
        end
    end
end

tasks = {'SCM','LIST','VPLT','SSCM'};

task_dirs = {
    'BCRW Salience Model\SCM Image Sets\';
    'List Task Analysis\List Image Sets\';
    'VPLT\ImageFiles\';
    'SSCM\'
    };

Matlab_dir = ['C:\Users\seth.koenig\Documents\MATLAB\'];
IOR_dir = [Matlab_dir 'IOR\Eye Data\'];
cd(IOR_dir)

tasks = {'SCM','LIST','VPLT','SSCM'};

distancethreshold = [0  24 48 72 96  120  144 168 200 400];%in pixels 24 pixels/dva
distancethreshold = [distancethreshold;
    [24 48 72 96 120 144  168 200 400 800]];%in pixels 24 pixels/dva

totalfixations = zeros(length(tasks),size(distancethreshold,2));
allSalIOR = cell(length(tasks),size(distancethreshold,2));
for task = 1:length(tasks);
    cd([IOR_dir tasks{task}])
    matfiles = what;
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'SalienceIOR'))
            load(matfiles.mat{i});
            for ii = 1:size(allSalIOR,2);
                allSalIOR{task,ii} = [allSalIOR{task,ii}; returnfixsal{ii}];
                totalfixations(task,ii) = totalfixations(task,ii)+size(returnfixsal{ii},1);
            end
        end
    end
end

for task = 1:length(tasks);
    labels = {};
    for d = 1:size(allSalIOR,2)
        return_sal_means(2,d) = nanmean(allSalIOR{task,d}(:,8));
        return_sal_stds(2,d) = nanstd(allSalIOR{task,d}(:,8))/sqrt(sum(~isnan(allSalIOR{task,d}(:,8))));
        return_sal_means(1,d) = nanmean(allSalIOR{task,d}(:,4));
        return_sal_stds(1,d) = nanstd(allSalIOR{task,d}(:,4))/sqrt(sum(~isnan(allSalIOR{task,d}(:,4))));
        labels{d} = [num2str(distancethreshold(1,d)) '-' num2str(distancethreshold(2,d))];
    end
    
    figure
    hold on
    errorbar(return_sal_means(2,:),return_sal_stds(2,:))
    errorbar(return_sal_means(1,:),return_sal_stds(1,:),'g')
    title(tasks{task})
    xlabel('Distance Tresholds')
    ylabel('Salience')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    hold off
    ylim([0.25 0.45])
    legend('Salience @ Return Fixation','Salience @ Prior Fixation')
end

for task = 1:length(tasks);
    labels = {};
    for d = 1:size(allSalIOR,2)
        return_dur_means(1,d) = nanmean(5*allSalIOR{task,d}(:,10));
        return_dur_stds(1,d) = nanstd(5*allSalIOR{task,d}(:,10))/sqrt(sum(~isnan(allSalIOR{task,d}(:,10))));
        return_dur_means(2,d) = nanmean(5*allSalIOR{task,d}(:,11));
        return_dur_stds(2,d) = nanstd(5*allSalIOR{task,d}(:,11))/sqrt(sum(~isnan(allSalIOR{task,d}(:,11))));
        labels{d} = [num2str(distancethreshold(1,d)) '-' num2str(distancethreshold(2,d))];
    end
    
    figure
    hold on
    errorbar(return_dur_means(2,:),return_dur_stds(2,:))
    errorbar(return_dur_means(1,:),return_dur_stds(1,:),'g')
    title(tasks{task})
    xlabel('Distance Tresholds')
    ylabel('Fixation Duration (ms)')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    hold off
    legend('Fix dur @ Return Fixation','Fix dur @ Prior Fixation')
end

rs = [];
ps = [];
slopes = [];
for task = 1:length(tasks);
    dist = [];
    sal = [];
    for i = 1:2
        for ii = 1:size(allSalIOR{task,i},1)
            dist = [dist ; allSalIOR{task,i}(ii,9)];
            sal = [sal; allSalIOR{task,i}(ii,8)];
        end
    end
    
    P = polyfit(dist,sal,1);
    figure
    hold on
    plot(dist,sal,'.')
    plot(1:50,[1:50]*P(1)+P(2),'r')
    hold off
    title(tasks{task})
    xlabel('Distance (pixels)')
    ylabel('Normalized Salience')
    [r,p]= corrcoef(dist,sal);
    
    slopes = [slopes P(1)];
    rs = [rs r(2)];
    ps = [ps p(2)];
end

for task = 1:length(tasks);
    labels = {};
    for d = 1:size(allSalIOR,2)
        mean_time_btwn_return(d) = 5*mean(allSalIOR{task,d}(:,7)-allSalIOR{task,d}(:,3)-1);
        mean_numfixations_btwn_return(d) = mean(allSalIOR{task,d}(:,13)-allSalIOR{task,d}(:,12)-1);
        
        std_time_btwn_return(d) = 5*std(allSalIOR{task,d}(:,7)-allSalIOR{task,d}(:,13)-1);
        std_numfixations_btwn_return(d) = std(allSalIOR{task,d}(:,13)-allSalIOR{task,d}(:,12)-1);
        
        numreturns(d) = sum(~isnan(allSalIOR{task,d}(:,13)));
        
        labels{d} = [num2str(distancethreshold(1,d)) '-' num2str(distancethreshold(2,d))];
    end
    
    figure
    subplot(1,2,1)
    errorbar(mean_time_btwn_return,std_time_btwn_return./sqrt(numreturns));
    xlabel('Distance Tresholds')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    ylabel('Duration between prior and return fixation')
    
    subplot(1,2,2)
    errorbar(mean_numfixations_btwn_return,std_numfixations_btwn_return./sqrt(numreturns))
    xlabel('Distance Tresholds')
    set(gca,'XTick',1:length(distancethreshold))
    set(gca,'XTickLabel',labels)
    ylabel('Number of fixations in-between prior and return fixation')
    subtitle(tasks{task})
end
%% --- [6] Run Seperate Special Salience IOR function for SSCM ---%
% done so I can know which image each return belongs to so I can related 
% returns to the content of that image and I am only
% interested in true return fixations (0-2dva apart). 
Matlab_dir = ['C:\Users\seth.koenig\Documents\MATLAB\'];
IOR_dir = [Matlab_dir 'IOR\'];
cd(IOR_dir)

tasks = {'SSCM'};

image_dirs = {
    'SSCM\'
    };

distancethreshold = [0;48];

imageX =[800];
imageY =[600];

for task = 1:length(tasks)
    cd([IOR_dir 'Eye Data\' tasks{task}])
    a = what;
    a = a.mat;
    for aa = 1:size(a,1);
        if length(a{aa}) == 14 || ~isempty(strfind(a{aa},'fixation'))
            load(a{aa})
            set_dir = [Matlab_dir image_dirs{task} 'S' num2str(imgset) '\'];
            disp(['Running Task ' tasks{task} ' fixiation file ' num2str(aa)])
            SSCM_SalienceIOR(a{aa},set_dir,imageX(task),imageY(task));
        end
    end
end

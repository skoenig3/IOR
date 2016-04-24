% central bias test repeat vs novel biases

Matlab_dir = 'C:\Users\seth.koenig\Documents\MATLAB\';
IOR_dir = [Matlab_dir 'IOR\'];

tasks = {'SCM','LIST','SMT','VPLT','SSCM'};

image_dirs = {
    'BCRW Salience Model\SCM Image Sets\';
    'List Task Analysis\List Image Sets\';
    'Scene Memory Task\Image Sets\';
    'VPLT\ImageFiles\';
    'SSCM\';
    };

imageX =[800 800 800 264 800];
imageY =[600 600 600 264 600];

fixation_PDF = cell(2,length(tasks));
for task = 1:length(tasks)
    fixation_PDF{1,task} = zeros(imageY(task),imageX(task));
    fixation_PDF{2,task} = zeros(imageY(task),imageX(task));
end

for task = 1:length(image_dirs);
    
    cd([IOR_dir 'Eye Data\' tasks{task}]);
    a = what;
    for m = 1:length(a.mat)
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
                fixationdata = fixationstats{imageprez(1,cnd)}.fixations;
            else
                fixationdata = fixationstats{novelconditions(cnd)}.fixations;
            end
            if ~isempty(fixationdata)
                if fixationdata(1,1) > imageX(task)/2-100 && fixationdata(1,1) < imageX(task)/2+100 &&...
                        fixationdata(2,1) < imageY(task)/2+100 && fixationdata(2,1) > imageY(task)/2-100
                    fixationdata(:,1) = [];
                end
                
                for f = 1:size(fixationdata,2)
                    fixx = fixationdata(1,f);
                    fixy = fixationdata(2,f);
                    fixx = round(fixx);
                    fixy = round(fixy);
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX(task)) = imageX(task);
                    fixy = imageY(task)-fixy;
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY(task)) = imageY(task);
                    fixation_PDF{1,task}(fixy,fixx) = fixation_PDF{1,task}(fixy,fixx)+1;
                end
            end
        end
        
        
        for cnd = 1:length(repeatconditions);
            if strcmpi(tasks{task},'VPLT')
                fixationdata = fixationstats{imageprez(2,cnd)}.fixations;
            else
                fixationdata = fixationstats{repeatconditions(cnd)}.fixations;
            end
            if ~isempty(fixationdata)
                if fixationdata(1,1) > imageX(task)/2-100 && fixationdata(1,1) < imageX(task)/2+100 &&...
                        fixationdata(2,1) < imageY(task)/2+100 && fixationdata(2,1) > imageY(task)/2-100
                    fixationdata(:,1) = [];
                end
                
                for f = 1:size(fixationdata,2)
                    fixx = fixationdata(1,f);
                    fixy = fixationdata(2,f);
                    fixx = round(fixx);
                    fixy = round(fixy);
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX(task)) = imageX(task);
                    fixy = imageY(task)-fixy;
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY(task)) = imageY(task);
                    fixation_PDF{2,task}(fixy,fixx) = fixation_PDF{2,task}(fixy,fixx)+1;
                end
            end
        end
    end
end

f = fspecial('gaussian',[256,256],24);
for task = 1:length(tasks);
    temp = fixation_PDF{1,task};
    temp = imfilter(temp,f);
    temp = temp/sum(sum(temp));
    fixation_PDF{1,task} = temp;
    
    temp = fixation_PDF{2,task};
    temp = imfilter(temp,f);
    temp = temp/sum(sum(temp));
    fixation_PDF{2,task} = temp;
end

for task  = 1:length(tasks)
    figure
    subplot(1,2,1)
    imagesc(fixation_PDF{1,task})
    axis equal
    axis off
    title('Novel Presentations')
    
    subplot(1,2,2)
    imagesc(fixation_PDF{2,task})
    axis equal
    axis off
    title('Repeat Presentations')
    
    subtitle(tasks{task})
end
    
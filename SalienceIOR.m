function SalienceIOR(FIXATIONFILE,distancethreshold,imagefiles,img_dir,imageX,imageY,novelconditions,min_dist)
% created by Seth Koenig 11/21/2012

% function determines rate of return fixations, the time between return
% fixations, time within trial of return, and salience at returned
% location. Inhibition of return was considered for mean fixation postions
% occruing within 0.5 dva of each other non-consequitively

% Inputs:
%   FIXATIONFILE: Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   ImageX,ImageY: x and y dimensions of images
%   Pairings: take all pairing or only closest unique pairing between first
%   and second fixation
%   img_dir: image director so can pull saliencemaps
%   imagefiles: image numbers so can pull salience maps
%   distancethreshold: distance categories for pairing fixations
%   novelconditions: novel condition numbers in fixationfile that represent novel
%   min_dist: minimum saccade distance out of area in dva   

%Outputs:
%   A .mat file named [FIXATIONFILE(1:end-13) '-SalienceIOR']
%   containg saccade statistics. See variable statvariablenames for
%   detailed explanation of variables in .mat file.

if nargin < 1
    error(['Not enough inputs: function requires FixationFile,'...
        'distance threhsold,imageX, imageY, and pairings.'])
end
if nargin < 2
    distancethreshold = [0  24 48 72 96  120  144 168 200 400];%in pixels 24 pixels/dva
    distancethreshold = [distancethreshold;
        [24 48 72 96 120 144  168 200 400 800]];%in pixels 24 pixels/dva
end
if nargin < 4
    imageX = 800;
    imageY = 600;
end

load(FIXATIONFILE);

returnfixsal = cell(1,size(distancethreshold,2));
%fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2
count = ones(1,size(distancethreshold,2));
for i = 1:size(distancethreshold,2)
    returnfixsal{i} = NaN(250,13);
end

for cndlop=1:length(novelconditions)
    if iscell(imagefiles)
        load([img_dir imagefiles{cndlop} '-saliencemap.mat'],'fullmap');
    else
        load([img_dir num2str(imagefiles(cndlop)) '-saliencemap.mat'],'fullmap');
    end
    saliencemap = fullmap;
    fixations = fixationstats{novelconditions(cndlop)}.fixations;
    if ~isempty(fixations)
        fixationtimes = fixationstats{novelconditions(cndlop)}.fixationtimes;
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixationtimes(:,1) = [];
        end
        
        N=size(fixations,2);
        if N > 1
            [x,y]=meshgrid(1:N);
            i=find(ones(N)-eye(N)); %forms pairs except for self-pairing
            i=[x(i), y(i)];
            i(i(:,1) > i(:,2),:) = []; %repeat pairs
            i(i(:,1)+1 == i(:,2),:) = []; %removes consecutive in time pairs
            dist =sqrt((fixations(1,i(:,1))-fixations(1,i(:,2))).^2 +...
                (fixations(2,i(:,1))-fixations(2,i(:,2))).^2);
            wcount = 1;
            pairs = NaN(ceil(size(fixations,2)/2),3);
            while ~isempty(i);
                [minn,mind] = min(dist);
                minn = minn(1);
                mind = mind(1);
                
                middlefixes = i(mind,1)+1:i(mind,2)-1;
                middist = sqrt((fixations(1,i(mind,1))-fixations(1,middlefixes)).^2+...
                    (fixations(2,i(mind,2))-fixations(2,middlefixes)).^2);
                
                if any(middist >= min_dist*24)
                    tind = find(minn > distancethreshold(1,:) & minn <= distancethreshold(2,:));
                    if ~isempty(tind);
                        pairs(wcount,:) = [i(mind,:) tind];
                    end
                    [rmvind1,~] = find(i(:,1) == i(mind,1));
                    [rmvind2,~] = find(i(:,2) == i(mind,1));
                    [rmvind3,~] = find(i(:,1) == i(mind,2));
                    [rmvind4,~] = find(i(:,2) == i(mind,2));
                    rmvind = [rmvind1; rmvind2; rmvind3; rmvind4];
                    rmvind = unique(rmvind);
                    i(rmvind,:) = [];
                    dist(rmvind) = [];
                    wcount = wcount+1;
                else
                    dist(mind) = [];
                    i(mind,:) = [];
                end
            end
            pairs(isnan(pairs(:,1)),:) = [];
            
            for i = 1:size(pairs,1);
                
                spot = [ceil(fixations(:,pairs(i,1))) ceil(fixations(:,pairs(i,2)))];
                spot(2,:) = imageY-spot(2,:);
                spott = [(fixationtimes(1,pairs(i,1))+fixationtimes(2,pairs(i,1)))/2 ...
                    (fixationtimes(1,pairs(i,2))+fixationtimes(2,pairs(i,2)))/2];
                dist = sqrt((spot(1,1)-spot(1,2))^2+(spot(2,1)-spot(2,2))^2);
                spot(spot < 1) = 1;
                spot(1,spot(1,:) > imageX) = imageX;
                spot(2,spot(2,:) > imageY) = imageY;
                
                returnfixsal{pairs(i,3)}(count(pairs(i,3)),:) = [...
                    spot(1,1) spot(2,1) spott(1) saliencemap(spot(2,1),spot(1,1))...
                    spot(1,2) spot(2,2) spott(2) saliencemap(spot(2,2),spot(1,2))...
                    dist diff(fixationtimes(:,pairs(i,1)))+1 ...
                    diff(fixationtimes(:,pairs(i,2)))+1 pairs(i,1) pairs(i,2)];
                %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2
                
                count(pairs(i,3)) = count(pairs(i,3))+1;
            end
        end
    end
end
for i = 1:size(distancethreshold,2)
    returnfixsal{i}(isnan(returnfixsal{i}(:,1)),:) = [];
end

IORvariablenames = {
    'returnfixsal: [  %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal...';
    'fixdist fix1dur fix2dur fixnum1 fixnum2]';
    };

save([FIXATIONFILE(1:10) '-SalienceIOR.mat'],'returnfixsal','IORvariablenames')
end
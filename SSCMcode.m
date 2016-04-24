mkdir('C:\Users\seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM');

% Load and Analyzed Behavioral Data from SSCM 

dirRoot= 'R:\Buffalo Lab\eblab\Cortex Programs\SSCM\';
dirDat=[dirRoot 'Data Files\'];
dirItm=[dirRoot 'Itm Files\'];
dirCnd=[dirRoot 'Cnd Files\'];

imageX = 800;
imageY = 600;
samprate = 5;

% What subjects have seen each set?
sList={...
    %Set1
    {'JN','MP','IW','TT'};...
    %Set2
    {'JN','MP'};...
    %Set3
    {'JN','MP'}...
    %Set4
    {'JN','MP'}...
    %Set5
    {'JN','MP'}...
    %Set6
    {'JN','MP'}...
    %Set7
    {'JN','MP'}...
    };

% Data Files separated by set and subject
fList={...
    %Set1
    {{'JN120717.4','JN130513.3'},{'MP120716.1','MP130703.4'},{'IW120425.1'},{'TT120717.5'}}...
    %Set2
    {{'JN121108.2','JN121212.2'},{'MP121206.2','MP130704.2'}}...
    %Set3
    {{'JN121115.2','JN121213.2'},{'MP121207.2','MP130705.2','MP130705.4'}}...
    %Set4
    {{'JN121116.2','JN130514.2'},{'MP121210.2','MP130708.2'}}...
    %Set5
    {{'JN121119.2','JN130515.2'},{'MP121211.2','MP130709.2'}}...
    %Set6
    {{'JN121120.2','JN130516.2'},{'MP121212.2','MP130710.3'}}...
    %Set7loo
    {{'JN121121.2','JN130517.2'},{'MP121213.2','MP130711.2'}}...
    };

% What Sets Were Shown After Delivering 48 IU of OT? (24 IU/mL) Delivered
% for 10 Minutes, OT=1, SL=0 NaN=No nebulizer
indOT=cell(7,1);
indOT={...
    %Set1 'JN','MP','IW','TT'
    {{nan,0},{nan,0},{nan},{nan}}...
    %Set2 'JN','MP'
    {{1,0},{0,1}}...
    %Set3 'JN','MP'
    {{0,1},{1,0,0}}...
    %Set4 'JN','MP'
    {{1,0},{0,1}}...
    %Set5 'JN','MP'
    {{0,1},{1,0}}...
    %Set6 'JN','MP'
    {{1,0},{0,1}}...
    %Set7 'JN','MP'
    {{0,1},{1,0}}...
    };

% Files With cchgrid calibration
gridCalList={'JN130513.3';'JN130514.2';'JN130515.2';'JN130516.2';'JN130517.2';...
    'MP130704.2';'MP130705.2';'MP130705.4';'MP130708.2';'MP130709.2';'MP130711.2'};
% Calibration Files
gridCalFiles{1,1}{1}{2}={'JN130513.1','JN130513.4'};
gridCalFiles{4,1}{1}{2}={'JN130514.1'};
gridCalFiles{5,1}{1}{2}={'JN130515.1','JN130515.3'};
gridCalFiles{6,1}{1}{2}={'JN130516.1','JN130516.3'};
gridCalFiles{7,1}{1}{2}={'JN130517.1','JN130517.3'};

gridCalFiles{1,1}{2}{2}={'MP130703.3','MP130703.5'};
gridCalFiles{2,1}{2}{2}={'MP130704.1','MP130704.3'};
gridCalFiles{3,1}{2}{2}={'MP130705.1','MP130705.5'};
gridCalFiles{3,1}{2}{3}={'MP130705.1','MP130705.5'};
gridCalFiles{4,1}{2}{2}={'MP130708.1','MP130708.3'};
gridCalFiles{5,1}{2}{2}={'MP130709.1','MP130709.3'};
gridCalFiles{7,1}{2}{2}={'MP130711.3'};


% What .itm and .cnd files?
itmFileList={...
    %Set1
    'SSCM90.itm';...
    %Set2
    'SSCMS2.itm';...
    %Set3
    'SSCMS3.itm';...
    %Set4
    'SSCMS4.itm';...
    %Set5
    'SSCMS5.itm';...
    %Set6
    'SSCMS6.itm';...
    %Set7
    'SSCMS7.itm';...
    };
cndFileList={...
    %Set1
    'SSCM90_ERROR.cnd';...
    %Set2
    'SSCMS2.cnd';...
    %Set3
    'SSCMS3.cnd';...
    %Set4
    'SSCMS4.cnd';...
    %Set5
    'SSCMS5.cnd';...
    %Set6
    'SSCMS6.cnd';...
    %Set7
    'SSCMS7.cnd';...
    };

calX=[0,0,-375,375,0,0,-750,750,0];
calY=[0,375,0,0,-375,750,0,0,-750];

for setloop=3:size(fList,1); % Loop through every Set
    for subjloop=1:size(fList{setloop,1},2) % Loop through every subject with data for the current Set
        for fileloop=1:size(fList{setloop,1}{subjloop},2);% Loop through every Data file in the current Set
            if indOT{setloop}{subjloop}{fileloop} == 0 %if saline trial
                fileName=fList{setloop,1}{subjloop}{fileloop};
                datfil=[dirDat fileName]; % Select data file
                
                itmFile=[dirItm itmFileList{setloop}];
                
                if strcmp('JN121108.2',fileName);
                    cndFile=[dirCnd,'SSCMS2_ERROR.cnd'];
                elseif strcmp('JN130513.3',fileName);
                    cndFile=[dirCnd,'SSCM90.cnd'];
                else
                    cndFile=[dirCnd cndFileList{setloop}];
                end
                
                itmfil=[];
                [fid,message]=fopen(itmFile, 'r');
                if fid<0
                    disp(message);
                else
                    while 1
                        tline = fgetl(fid);
                        if ~isempty(itmfil)
                            if length(tline)>size(itmfil,2)
                                tline=tline(1:size(itmfil,2));
                            end
                        end
                        tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
                        if ischar(tline)
                            itmfil=[itmfil; tline];
                        else
                            break
                        end
                    end
                end
                fclose(fid);
                
                cndfil=[];
                [fid,message]=fopen(cndFile, 'r');
                if fid<0
                    disp(message);
                else
                    while 1
                        tline = fgetl(fid);
                        if ~isempty(cndfil)
                            if length(tline)>size(cndfil,2)
                                tline=tline(1:size(cndfil,2));
                            end
                        end
                        tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
                        if ischar(tline)
                            cndfil=[cndfil; tline];
                        else
                            break
                        end
                    end
                end
                fclose(fid);
                
                calItm=[dirItm 'cchgrid.itm'];
                calitmfil=[];
                [fid,message]=fopen(calItm, 'r');
                if fid<0
                    disp(message);
                else
                    while 1
                        tline = fgetl(fid);
                        if ~isempty(calitmfil)
                            if length(tline)>size(calitmfil,2)
                                tline=tline(1:size(calitmfil,2));
                            end
                        end
                        tline = [tline ones(1,(size(calitmfil,2)-length(tline)))*char(32)];
                        if ischar(tline)
                            calitmfil=[calitmfil; tline];
                        else
                            break
                        end
                    end
                end
                fclose(fid);
                
                calcndfil = [dirCnd 'cchDrew.cnd'];
                [fid,message]=fopen(calcndfil, 'r');
                if fid<0
                    disp(message);
                else
                    while 1
                        tline = fgetl(fid);
                        if ~isempty(calcndfil)
                            if length(tline)>size(calcndfil,2)
                                tline=tline(1:size(calcndfil,2));
                            end
                        end
                        tline = [tline ones(1,(size(calcndfil,2)-length(tline)))*char(32)];
                        if ischar(tline)
                            calcndfil=[calcndfil; tline];
                        else
                            break
                        end
                    end
                end
                fclose(fid);
                
                if ismember(fileName,gridCalList) %if using grid calibration
                    x = cell(1,63);
                    y = cell(1,63);
                    for calFile = 1:length(gridCalFiles{setloop,1}{subjloop}{fileloop});
                        [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata([dirDat gridCalFiles{setloop,1}{subjloop}{fileloop}{calFile}]);
                        
                        calbeg=500;
                        calend=900;
                        % Get eyedata for calibration with clrchng trials
                        numrpt = size(event_arr,2);
                        valrptcnt = 0;
                        clear per clrchgind
                        for rptlop = 1:numrpt
                            if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) <= 1189)) ~=0
                                if size(find(event_arr(:,rptlop) == 200)) ~=0
                                    perbegind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                                    perendind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                                    cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                                    blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                                    % Sample calbeg ms after appearance of gray fixspot
                                    % until calend ms after appearance of gray fixspot
                                    begtimdum = time_arr(perbegind,rptlop)+calbeg;
                                    endtimdum = time_arr(perendind,rptlop)+calend;
                                    if endtimdum > begtimdum
                                        valrptcnt = valrptcnt + 1;
                                        clrchgind(valrptcnt)=rptlop;
                                        per(valrptcnt).begsmpind = begtimdum;
                                        per(valrptcnt).endsmpind = endtimdum;
                                        per(valrptcnt).begpos = 1;
                                        per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                        per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                        per(valrptcnt).allval = event_arr(:,rptlop);
                                        per(valrptcnt).alltim = time_arr(:,rptlop);
                                    end
                                end
                            end
                        end
                        
                        clear cnd cndlst
                        numrpt = size(per,2);
                        for rptlop = 1:numrpt
                            cnd(rptlop)=per(rptlop).cnd;
                        end
                        cndList{calFile} = cnd;
                        
                        evnnmb=2:2:size(eog_arr,1);
                        oddnmb=1:2:size(eog_arr,1);
                        
                        % Create structures x and y of the corresponding average eye data for each trial
                        % instance (p) of each condition (k)
                        cndlst=unique(cnd);
                        for k=1:length(cndlst)
                            cndind=find(cnd==cndlst(k));
                            allind=clrchgind(cndind);
                            for p=1:length(allind)
                                xi{k}(p)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,oddnmb),allind(p)));
                                yi{k}(p)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,evnnmb),allind(p)));
                            end
                        end
                        for kk = 1:length(xi);
                            x{kk} = [x{kk} xi{kk}];
                            y{kk} = [y{kk} yi{kk}];
                        end
                    end
                    
                    % Values that are more than std away from median are removed
                    clear meanx meany
                    for k=1:numel(x)
                        xss = x{k};
                        low = mean(xss)-2*std(xss);
                        high = mean(xss)+2*std(xss);
                        xss(xss < low) = [];
                        xss(xss > high) = [];
                        meanx(k)=median(xss);
                        
                        yss = y{k};
                        low = mean(yss)-2*std(yss);
                        high = mean(yss)+2*std(yss);
                        yss(yss < low) = [];
                        yss(yss > high) = [];
                        meany(k)=median(y{k});
                    end
                    
                    cndlst = unique(cell2mat(cndList))-1000;
                    control = [];
                    for cnd = 1:length(cndlst);
                        itm  = textscan(calcndfil(cnd + 2,:),'%s');
                        itmline = textscan(calitmfil(str2double(itm{1}(3))+5,:),'%s');
                        control(cnd,1) = str2double(itmline{1}(9));%xpos
                        control(cnd,2) = str2double(itmline{1}(10)); %ypos
                    end
                    
                    % Compute the polynomial transformation function
                    tform = cp2tform(control,[meanx' meany'],'polynomial',4);
                    tform.forward_fcn = tform.inverse_fcn;
                    
                                    figure
                                    hold on
                                    for cnd = 1:size(control,1);
                                        plot(control(cnd,1),control(cnd,2),'+');
                                        [cx,cy] = tformfwd(tform,meanx(cnd),meany(cnd));
                                        plot(cx,cy,'r*')
                                    end
                 
                    [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datfil);
                    
                else %if using 9 point built-in calibration
                    
                    [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datfil);
                    
                    calbeg=500;
                    calend=900;
                    % Get eyedata for calibration with clrchng trials
                    numrpt = size(event_arr,2);
                    valrptcnt = 0;
                    clear per clrchgind
                    for rptlop = 1:numrpt
                        if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 1010)) ~=0
                            if size(find(event_arr(:,rptlop) == 200)) ~=0
                                perbegind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                                perendind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                                % Sample calbeg ms after appearance of gray fixspot
                                % until calend ms after appearance of gray fixspot
                                begtimdum = time_arr(perbegind,rptlop)+calbeg;
                                endtimdum = time_arr(perendind,rptlop)+calend;
                                if endtimdum > begtimdum
                                    valrptcnt = valrptcnt + 1;
                                    clrchgind(valrptcnt)=rptlop;
                                    per(valrptcnt).begsmpind = begtimdum;
                                    per(valrptcnt).endsmpind = endtimdum;
                                    per(valrptcnt).begpos = 1;
                                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                    per(valrptcnt).allval = event_arr(:,rptlop);
                                    per(valrptcnt).alltim = time_arr(:,rptlop);
                                end
                            end
                        end
                    end
                    
                    clear cnd
                    numrpt = size(per,2);
                    for rptlop = 1:numrpt
                        cnd(rptlop)=per(rptlop).cnd;
                    end
                    
                    evnnmb=2:2:size(eog_arr,1);
                    oddnmb=1:2:size(eog_arr,1);
                    
                    clear x y
                    cndlst=unique(cnd);
                    for k=1:length(cndlst)
                        cndind=find(cnd==cndlst(k));
                        allind=clrchgind(cndind);
                        for p=1:length(allind)
                            x{k}(p)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,oddnmb),allind(p)));
                            y{k}(p)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,evnnmb),allind(p)));
                        end
                    end
                    
                    if strmatch('IW120425.1',fileName);
                        %Correct for error in cnd 2 and rename of cnd 5 to cnd 2 in IW120425.1
                        %Remove Condition 2 & 3
                        x={x{1,1},x{1,4},x{1,5},x{1,6},x{1,7},x{1,8},x{1,9}};
                        y={y{1,1},y{1,4},y{1,5},y{1,6},y{1,7},y{1,8},y{1,9}};
                        %Change Calibration Points: Remove Cnd 2 and put 5 at 0,6 degrees
                        calX=[0,  3,  0,   0,  -6, 6,  0];
                        calY=[0,  0,   6, 6,   0,   0,  -6];
                    elseif strmatch('JN121108.2',fileName);
                        %Remove Condition 3
                        x={x{1,1},x{1,2},x{1,4},x{1,5},x{1,6},x{1,7},x{1,8},x{1,9}};
                        y={y{1,1},y{1,2},y{1,4},y{1,5},y{1,6},y{1,7},y{1,8},y{1,9}};
                        %Change Calibration Points: Put Cnd 5 at 0,6 degrees
                        calX=[0,0,3,0,0,-6,6,0];
                        calY=[0,3,0,6,6,0,0,-6];
                    elseif strmatch('JN121212.2',fileName);
                        %Remove Condition 3
                        x={x{1,1},x{1,2},x{1,4},x{1,5},x{1,6},x{1,7},x{1,8},x{1,9}};
                        y={y{1,1},y{1,2},y{1,4},y{1,5},y{1,6},y{1,7},y{1,8},y{1,9}};
                        %Change Calibration Points: Put Cnd 5 at 0,6 degrees
                        calX=[0,0,3,0,0,-6,6,0];
                        calY=[0,3,0,6,6,0,0,-6];
                    else
                        calX=[0,0,-3,3,0,0,-6,6,0];
                        calY=[0,3,0,0,-3,6,0,0,-6];
                    end
                    
                    clear meanx meany
                    for k=1:numel(x)
                        xss = x{k};
                        low = mean(xss)-std(xss);
                        high = mean(xss)+std(xss);
                        xss(xss < low) = [];
                        xss(xss > high) = [];
                        meanx(k)=median(xss);
                        
                        yss = y{k};
                        low = mean(yss)-std(yss);
                        high = mean(yss)+std(yss);
                        yss(yss < low) = [];
                        yss(yss > high) = [];
                        meany(k)=median(y{k});
                    end
                    
                    tform = cp2tform([calX' calY'],[meanx' meany'],'affine');
                    tform.forward_fcn = tform.inverse_fcn;
                    
                    %                     figure
                    %                     hold on
                    %                     for cnd = 1:size(calX,2);
                    %                         plot(calX(cnd),calY(cnd),'+');
                    %                         [cx,cy] = tformfwd(tform,meanx(cnd),meany(cnd));
                    %                         plot(cx,cy,'r*')
                    %                     end

                end
                
                numrpt = size(event_arr,2);
                valrptcnt = 0;
                clear per vpcind
                new_eog_arr=[];
                for rptlop = 1:numrpt
                    % Find all stimulus presentation trials (cnd>=1010)
                    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
                        if size(find(event_arr(:,rptlop) == 200)) ~=0
                            % 23 instead of 24 because 23 is stim onset
                            perbegind = find(event_arr(:,rptlop) == 23,1,'first');
                            perendind = find(event_arr(:,rptlop) == 24,1,'first');
                            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                            begtimdum = time_arr(perbegind,rptlop);
                            endtimdum = time_arr(perendind,rptlop);
                            if endtimdum > begtimdum
                                valrptcnt = valrptcnt + 1;
                                vpcind(valrptcnt)=rptlop;
                                per(valrptcnt).begsmpind = begtimdum;
                                per(valrptcnt).endsmpind = endtimdum;
                                per(valrptcnt).begpos = 1;
                                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                per(valrptcnt).allval = event_arr(:,rptlop);
                                per(valrptcnt).alltim = time_arr(:,rptlop);
                                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                            end
                        end
                    end
                end
                
                clear eyedat
                for trlop=1:size(per,2)
                    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
                    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
                    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
                    picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to iscan start
                    picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to iscan start
                    
                    if picstart == 0;
                        picstart = 5;
                    end
                    if picend/5<=length(horeog)
                        eyedat{trlop}(1,:) = (horeog(ceil(picstart/5):floor(picend/5)));
                        eyedat{trlop}(2,:) = (vrteog(ceil(picstart/5):floor(picend/5)));
                    else
                        eyedat{trlop}(1,:)=nan;
                        eyedat{trlop}(2,:)=nan;
                    end
                end
                
                cnd=[];
                numrpt = size(per,2);
                for rptlop = 1:numrpt
                    cnd(rptlop)=per(rptlop).cnd;
                end
                
                %---Recalibrate and automatically scale eye data---%
                for eye = 1:length(eyedat)
                    x = eyedat{eye}(1,:);
                    y = eyedat{eye}(2,:);
                    [x,y] = tformfwd(tform,x,y);
                    eyedat{eye} = [x;y];
                end
                
                notattended = NaN(1,size(eyedat,2));
                for i = 1:size(eyedat,2);
                    x = 24*eyedat{i}(1,:)+imageX/2;
                    y = 24*eyedat{i}(2,:)+imageY/2;
                    if length(x) > 2120 %500 ms cross hair + 7000 ms for image + 100 ms buffer
                        notattended(i) = (length(x)-2120)/2120;
                        x = x(1:2120);
                        y = y(1:2120);
                    else
                        notattended(i) = 0;
                    end
                    badx = find(x < -50 | x > imageX+50); %~1 dva leave margin of error
                    x(badx) = []; y(badx) = [];
                    bady = find(y < -50 | y > imageY+50); %~1 dva margin of error
                    x(bady) = []; y(bady) = [];
                    eyedat{i} = [x;y];
                end
                
                images = cell(1,length(eyedat)); %complex image names
                for i = 1:length(cnd);
                    C = textscan(cndfil(cnd(i)-1000+1,:),'%s');
                    C = textscan(itmfil(str2double(C{1}{5})+6,:),'%s');
                    img = C{1}{end};
                    slash = strfind(img,'\');
                    img = img(slash(2)+1:end-4);
                    images{i} = img;
                end
                
                imgset = setloop;
                %fixationstats = ClusterFixation_Final(eyedat);
                %save(['C:\Users\seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM\' ...
                  %  fileName(1:end-2) '_' fileName(end) '.mat'],'fixationstats',...
                  % 'images','imgset');
            end
        end
    end
end
function data = Salience_at_fixations(fixationdata,saliencemap,fixdurlow,fixdurhigh)
% created by Seth Koenig 1/21/2014. Simplified and modified from previous
% versions: fixation_Salience_and_Significance.mat.

% function determines normalized salience at fixation locations.

% fixationdata: fixationstats{cnd}

imageX = size(saliencemap,2);
imageY = size(saliencemap,1);

fixationtimes = fixationdata.fixationtimes;
if ~isempty(fixationtimes)
    fixations = fixationdata.fixations;

    numfixs = size(fixations,2);
    salience = NaN(1,numfixs);
    randomsalience = NaN(1,numfixs);
    
    for i = 1:numfixs
        fixdur = 5*(fixationtimes(2,i) - fixationtimes(1,i))+5; %samples to ms
        
        if (fixdur >= fixdurlow) && (fixdur <= fixdurhigh)
            
            spot = ceil(fixations(:,i));
            spot(2) = imageY-spot(2);
            spot(spot < 1) = 1;
            spot(1,spot(1) > imageX) = imageX;
            spot(2,spot(2) > imageY) = imageY;
            
            rspot = [ceil(imageX*rand) ceil(imageY*rand)]; %fake x,y data
            
            salience(i) = saliencemap(spot(2),spot(1));
            randomsalience(i) = saliencemap(rspot(2),rspot(1));
        end
    end
    
    data.salience = salience;
    data.shuffled = randomsalience;
end

end
function Spikecuttrials2

files = dir('datahighpassch*.mat');

filesB = {files.name};

load('events001','trials')

% %-------
% trials(1047) = [];
% %-------

l = 1;
% dubinderror = 0;
for k = 1 : length(trials)
    if trials(1,k).header.success == true
        clear states
        states(:,1) = cell2mat({trials(1,k).states.con});
        states(:,2) = cell2mat({trials(1,k).states.time});
        states(logical(cat(1,0,diff(states(:,1)) == 0)),:) = [];
        if length(states) > 10 
            trialstart(l) = round(states(states(:,1) == 0,2)*1000);
%             hand(l) = round(states(states(:,1) == 1,2)*1000);
% %             motor(l) = round(states(states(:,1) == 4,2)*1000);
%             fix(l) = round(states(states(:,1) == 5,2)*1000);
%             cue1(l) = round(states(states(:,1) == 6,2)*1000);
%             mem1(l) = round(states(states(:,1) == 7,2)*1000);
%             cue2(l) = round(states(states(:,1) == 17,2)*1000);
%             mem2(l) = round(states(states(:,1) == 18,2)*1000);
%             reac(l) = round(states(states(:,1) == 8,2)*1000);
%             move(l) = round(states(states(:,1) == 9,2)*1000);
            inirew(l) = round(states(states(:,1) == 10,2)*1000);
%             type{l} = trials(1,k).header.graspType;
%             grip{l} = trials(1,k).header.selectedGraspType;
%             type2(l) = trials(1,k).header.Cue2;
            %             reward(l) = trials(1,k).header.Reward;
            %             powRat(l) = trials(1,k).header.powRat;
            l = l + 1;
%         else
%             dubinderror = dubinderror +1;
        end
    end
end

clear trials hand motor states fix cue1 mem1 cue2 mem2 reac move type grip type2 dubinderror k l


for i = 1 : length(filesB)
    tic
    load(filesB{i});
    
    
    data = double(data);
    
    namefile = filesB{i}
    
    finalsize = 0;
    for j = 1 : length(trialstart)
        finalsize = finalsize + inirew(j)*30-trialstart(j)*30+1;
    end
    
    datanew = zeros(finalsize,1);
    s = 1;
    for j = 1 : length(trialstart)
        e = s+(inirew(j)*30)-(trialstart(j)*30);
        datanew(s:e) = data(trialstart(j)*30:inirew(j)*30);
        s = e+1;
    end
    
    clear data 
    
    data = int16(datanew);
    
    clear datanew
    
    save(['datacut' namefile(end-8:end-4)],'data')
    clear data
    toc
end
        
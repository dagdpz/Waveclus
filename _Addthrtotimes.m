function Addthrtotimes_

channels = dir('times_datacut*.mat');
channels = {channels.name};

for i = 1 : length(channels)
    filename = channels{1,i};
    load([filename(7:end-4) '_spikes.mat'],'thr')
    save(filename,'thr','-append')
    clear thr
end
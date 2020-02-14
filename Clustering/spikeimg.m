function spikeimg2(spikes,cluster_class)



dummy = max(cluster_class(:,1));
if dummy > 5, dummy = 5; end

spikes2 = spikes;
% clip amplitudes
amplim = round([6*std(min(spikes2,[],2)), 6*std(max(spikes2,[],2))]);
ampmean = round([mean(min(spikes2,[],2)), mean(max(spikes2,[],2))]);
ll = ampmean(1)-amplim(1);
ul = ampmean(2)+amplim(2);
spikes2(spikes2 < ll) = ll;
spikes2(spikes2 > ul) = ul;
% round amplitues to integers and shift into positive values
spikes2 = round(spikes2);
shift = abs(min(spikes2(:))) + 1;
spikes2 = spikes2 + shift;

% define the dimensions of the images
r = max(spikes2(:));
c = size(spikes2,2);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
ls = [0,0,0;colormap('lines')];
lc = {'k','b','r','g','c','m','y'};
colormap([1,1,1;colormap('hot')])


% specify the desired number of YTicks and calculate
% how often we'll need a tick point (to the nearest 10)
numticks = 8;
tickint = (ul+shift)/numticks - mod((ul+shift)/numticks,10);
for j = 0:dummy
    % get the current cluster
    rspikes = spikes2(cluster_class(:,1) == j,:);
    
    % create the image
    cols = repmat(1:c,size(rspikes,1),1);
    inds = sub2ind([r,c],rspikes(:),cols(:));    
    spikesimg = accumarray(inds,1,[r*c,1]);
    spikesimg = reshape(spikesimg,r,c);
    clear inds cols
    
    % calculate interspike intervals
    interspk = cluster_class(cluster_class(:,1) == j,2);
    interspk = diff(interspk);
    interspk(interspk > 100) = [];

    % plot the image
    subplot(2,6,(j+1));
    imagesc(interp2(spikesimg,2))
    title(['Cluster ' num2str(j)],'FontSize',20,'Color',ls(j+1,:))
    % set the XTicks
    sf = get(gca,'XLim');
    sf = sf(2)/(c+sf(1));
    set(gca,'XTick',floor(sf*(0:10:c)));
    set(gca,'XTickLabel',num2cell(-20:10:c-20),'XMinorGrid','on','YDir','normal');
    % calculate and set the YTicks
    yl = get(gca,'ylim');
    imrange = yl(1):yl(2);
    sf = yl(2)/(shift+ul);
    zerotick = round((0-ll+1)*sf);
    zeroval = imrange(zerotick);
    newticks = [fliplr(zeroval:-sf*tickint:1),zeroval+sf*tickint:tickint*sf:yl(2)];
    pivot = find(newticks==zeroval);
    newlabels = [fliplr(0:-tickint:-tickint*(pivot-1)), ...
                 tickint:tickint:tickint*length(newticks)-pivot];
    set(gca,'ytick',newticks,'yticklabel',newlabels,'YGrid','on')

    % plot the interspike intervals
    subplot(2,6,(j+7));
    hist(interspk,100)
    isih = findobj(gca,'Type','patch');
    set(isih,'FaceColor',lc{j+1},'EdgeColor',lc{j+1})
    if ~isempty(max(interspk)), set(gca,'XLim',[0,max(interspk)]), end
    title([num2str(length(rspikes)) ' ' num2str(sum(interspk<2)) ' < 2ms; ' ...
           num2str(sum(interspk<1)) ' < 1ms'],'FontSize',12)
end
% clean up
clear dummy spikes2 amplim ampmean ll ul shift r c lc numticks tickint rspikes spikesimg ...
    interspk sf yl imrange zerotick zeroval newticks pivot newlabels isih
%% INIT and load AVI file

[fn, pn] = uigetfile('*.avi', 'pick a movie file');
vidObj = VideoReader([pn fn]);
temp = readFrame(vidObj);
 

%% read first few frames and get background

nFrames = 250; %about 10 seconds
im = zeros([size(temp, 1) size(temp, 2), nFrames], 'uint8');
vidObj.CurrentTime = 0;
for i = 1:nFrames 
    disp(i)
    temp = readFrame(vidObj);
    im(:,:,i) = temp(:,:,1);
end
background = median(im, 3);


%% draw rectangular MASK around arena (and leave this figure open)

figure, hIm = imagesc(background); colormap gray
ir = imrect;
hold on
hP = plot(0,0, 'r.');
hQ = quiver(0,0,0,0, 0);
hold off


%% playback & tracking

threshold = 5;
minarea = 20;

set(gca, 'CLim', [0 threshold*20])
mask = ir.createMask;
vidObj.CurrentTime = 0;
out = []; %this is going to be the structure that contains all extracted data
out.mask = mask;
out.fn = [pn fn];
background = double(background);
backgroundDecay = 1/300; %set the time constant for udating the background estimate

k = 1; %this will be the frame index

while hasFrame(vidObj)
    frame = readFrame(vidObj);
    frame = frame(:,:,1); %force image to be monochrome not RGB
    background = background * (1-backgroundDecay) + backgroundDecay*double(frame);
    foreground = background - double(frame);
    %fall(:,:,k) = foreground;
    
    %operations on binarized frame to extract regions and properties 
    temp = foreground >= threshold; 
    temp = temp .* mask;
    temp = bwareaopen(temp, minarea, 8);
    temp = bwlabel(temp);
    %temp = temp.*foreground;
    s = regionprops(temp, uint8(foreground),'centroid', 'area', 'orientation', 'Eccentricity', 'Solidity', 'ConvexArea', 'WeightedCentroid');
    
    out(k).centroids = cat(1, s.Centroid);
    out(k).solids = cat(1, s.Solidity);
    out(k).oris = cat(1, s.Orientation);
    out(k).careas = cat(1, s.ConvexArea);
    out(k).ecces = cat(1, s.Eccentricity);
    out(k).wcentroids = cat(1, s.WeightedCentroid);

    %exclude ellipses which have low eccentricity and solidity (or rather, include the others for further processing)
    idx = (out(k).ecces > 0.9) & (out(k).solids > 0.5);
    
    %this finds the direction vector of each fish
    vec = out(k).wcentroids - out(k).centroids; 
    vec = 20 * vec ./ sqrt(vec(:, 1).^2 + vec(:, 2).^2) ;
    
    %show video frame and vectors for preview
    set(hIm, 'CData', foreground) 
    set(hP, 'XData', out(k).centroids(idx,1), 'YData', out(k).centroids(idx,2))
    set(hQ, 'XData', out(k).centroids(idx, 1), 'YData', out(k).centroids(idx,2), 'UData', vec(idx, 1), 'VData', vec(idx,2))

    drawnow, pause(.01)
    k = k+1
end

%save([pn fn '_analysis.mat'], 'out')


%% analysis of spatial clustering (requires presence of 'out' structure (see above)

nnd = [];
allcentroids = [];

m = 6; %number of random NND samples to take (see CSR test, set to half of the fish number)

%collect all centroids across frame (used for shuffling below)
for iFrame = 1:numel(out)
    idx = (out(iFrame).ecces > 0.9) & (out(iFrame).solids > 0.5);
    allcentroids = [allcentroids; out(iFrame).centroids(idx, :)]; 
    iFrame
end

%collect average NND values for each frame (regular and shuffled) 
for iFrame = 1:numel(out)
    idx = (out(iFrame).ecces > 0.9) & (out(iFrame).solids > 0.5);
    
    pd = squareform(pdist(out(iFrame).centroids(idx, :))); %matrix of pairwise distances
    pd(pd==0) = Inf; %this will ignore the zeros along diagonal (distance to self)
    nnd(iFrame) = mean(randsample(min(pd), min([m numel(min(pd))])));
    
    %shuffled data
    pd0 = squareform(pdist(allcentroids(ceil(rand(numel(idx), 1)*size(allcentroids,1)), :)));
    pd0(pd0==0) = Inf;
    nnd0(iFrame) = mean(randsample(min(pd0), min([m numel(min(pd0))])));
    
    iFrame
end


%% CSR test on NND values
%see http://www.seas.upenn.edu/~ese502/NOTEBOOK/Part_I/3_Testing_Spatial_Randomness.pdf

d = nanmean(nnd)/26; %26 pixels per cm
numfish = (size(allcentroids, 1)/15000); %average number of fish per frame (taking into account any missing data due to tracking limitations)
lambda = numfish/(24^2); %average fish density (number per 24 cm sqare arena)
mu = 1/(2*sqrt(lambda)); %expected average nnd assuming CSR
sigma = sqrt((4 - pi)/(numfish*pi*4*lambda));
z = (d-mu)/sigma;
p = normcdf(z); %p-Value of CSR test

figure
edges = 0:.1:8;
stairs(edges, [histc(nnd()/26, edges)' histc(nnd0/26, edges)'])
hold on
plot(mu, 0, 'r.')

title(['CSR nnd: ' num2str(mu) '   actual nnd: ' num2str(d) '   p: ' num2str(p)])
xlabel('distance in cm'), ylabel('frames')

%saveas(gcf, [out(1).fn '_nnd.fig'])
%saveas(gcf, [out(1).fn '_nnd.svg'])


%% analysis of orientations

thetas = {};
allvec = cat(1,out.centroids) - cat(1,out.wcentroids);

for iFrame = 1:numel(out)
    idx = (out(iFrame).ecces > 0.9) & (out(iFrame).solids > 0.5);
    vec = out(iFrame).centroids - out(iFrame).wcentroids;
    vec0 = allvec(ceil(rand(size(vec, 1), 1)*size(allvec, 1)), :); %shuffled across frames
    
    phasor = vec(idx, 1) + 1i * vec(idx,2);
    phasor0 =  vec0(:,1) + 1i * vec0(:,2);
    thetas(iFrame) = {angle(conj(mean(phasor)) .* phasor)};
    thetas0(iFrame) = {angle(conj(mean(phasor0)) .* phasor0)};
    iFrame
end

figure
polarhistogram(cat(1, thetas{:}),36*2)
ax = gca; ax.ThetaZeroLocation = 'top';
title(out(1).fn)

drawnow
%saveas(gcf, [out(1).fn '_ori.fig'])
%saveas(gcf, [out(1).fn '_ori.svg'])

figure
polarhistogram(cat(1, thetas0{:}),36*2)
ax = gca; ax.ThetaZeroLocation = 'top';
title([out(1).fn ' shuffled'])

drawnow
%saveas(gcf, [out(1).fn '_ori_shuff.fig'])
%saveas(gcf, [out(1).fn '_ori_shuff.svg'])


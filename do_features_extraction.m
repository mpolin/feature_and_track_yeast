%Script for the first-level analysis of the images from the experiments of
%yeast digestion. The images are first segmented, then featured and finally
%the features are tracked. The features (see below) of each image are saved
%in individual files for subsequent use. As the program cycles through the
%series of images, the script assembles the array "tracks", which is then
%used as input for the tracking program ("track.m"), adapted from Crocker
%and Grier. After rearranging the output, the array "tracks" is saved in
%the same folder where the feature files are saved. 
%The array "tracks" has one feature per row, with the columns being
%organised as
%(CentroidXposition, CentroidYposition, Area, ..., frame, naiveID,
%trackedID). 
%Here the "..." are due to the fact that one can add in there more
%parameters characterising the features, if useful to have later on in the
%array "tracks".
%The naiveID is the ID of the feature when that frame is featured first,
%while trackedID is the ID of that feature once tracked. 
%For example, in order to follow the object with trackID equal to "M", one
%needs to select the rows in "tracks" that have the value "M" in the last
%column. The corresponding frames will be listed in the third-to-last
%column and the naiveID of that object for those frames will be the
%penultimate column. Notice that an object that has trackedID "M", will
%have in general several different naiveIDs for the different frames it
%exists in. From those rows one can then see the time evolution, for object
%"M", of any of the characteristics saved in "tracks" for object "M". If
%one needs to see other characteristics that are not currently saved in
%"tracks", they can read them from the files "XXXX_features.mat"
%corresponding to the frames in which object "M" exists. For example, if a
%row of "tracks" ends with 
%   (.....,T,Y,M)
%this means that object "M" exists in frame "T", and there it corresponds
%to the object with naiveID "Y". One would then read the
%"XXXX_features.mat" file corresponding to frame "T" and retrieve the
%structure "features" corresponding to that file. The features
%corresponding to the object with trackedID "M" would correspond to those
%of "features(Y)". So, for example "features(Y).outs" would give the
%contour length along the boundary, and "features(Y).theta" would give the
%corresponding tangent angle.
%
%HISTORY: 
%   5 April, 2024: MP. Created.
%
%
%TODO:
%1) Write programs to "pull" a given set of parameters for an object with a
%given trackID, and all of the frames for which trackedID exists. 
%2)try watershed on brightfield images of a single feature. should reveal
%the four cells in the tetrad 

%% Input and output info
imgpath = '/Users/mpolin/experiments/yeast/digestion/codes/images/';
datapath = [imgpath,'data/'];
imgextension = 'tif'; %make sure you use single quotation marks, e.g. 'tif', and not double, e.g. "tif". String concatenation is different

%% Parameter definition
%%Segmentation
invert = 1;
paramSegment.int_threshold = 4.325310e+04;
paramSegment.mode_threshold = 5e4;
paramSegment.arearange = [2000,inf];
paramSegment.morph_close_radius = 3;

%%Featuring
paramFeature.minlen = 4;
paramFeature.maxlen = inf;
paramFeature.b_init = 1;

%%Tracking
paramTracks.maxdisp = 10;
paramTracks.mem = 0;
paramTracks.good = 2;
paramTracks.dim = 2;
paramTracks.quiet = 1;
tracks = [];

%% Gather list of files and perform analysis
if ~isdir(datapath)
    mkdir(datapath)
end
myfiles = dir([imgpath,'*.',imgextension]);

for tt = 1:length(myfiles)
    
    if mod(tt-1,10)==0
        disp(['Extracting features from frame ',num2str(tt),' out of ',num2str(length(myfiles))]);
    end

    %read new image
    img = imread([imgpath,myfiles(tt).name]);
    [~,imgname,~] = fileparts(myfiles(tt).name);

    %invert if needed
    if invert
        img = imcomplement(img);
    end

    [BWimg,maskedImage] = segmentImage(img,paramSegment); %segmentation
    features = feature_connected_components(BWimg,paramFeature); %featuring
    tracks = [tracks;[vertcat(features.Centroid),... %assemble array that will be used to track the features
        vertcat(features.Area),... %repeat this structure in the following line if you want to add an extra column representing a numeric property of the features. This might be useful to have in the array "tracks" later on
        (1:length(features))',tt*ones(length(features),1)]];
    save([datapath,imgname,'_features.mat'],'features'); %save features of the objects in the current frame
end
disp('Done with the extraction. Moving onto tracking.')
tracks = track(tracks,paramTracks.maxdisp,paramTracks); %use the track array assembled earlier to track the objects across frames. This appends at the end of "tracks" a column with the trackedID of that specific object
tracks = tracks(:,[1:size(tracks,2)-3,size(tracks,2)-1,size(tracks,2)-2,size(tracks,2)]); %reorder the columns to get (xpos,ypos,Area,...,frame, naiveID, trackedID)
save([datapath,'tracks.mat'],'tracks');
disp('Done.')




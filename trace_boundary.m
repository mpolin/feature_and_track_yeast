function myborders=trace_boundary(img,minlen,maxlen,b_init)
%Function to order the pixels along a boundary. This is taken from the
%boundary tracing python routine from Sebastian Wallkotter available here
% https://github.com/FirefoxMetzger/ipynb_boundary_tracing/tree/master
%
%Besides porting it to Matlab, I made a slight change in the output. Each
%traced boundary is now outputted without immediately repeating positions.
%For more informations and for a review of the idea behind the algorithm,
%please look at the comment section at the end of this program.
%
%INPUT:
%   img: binary image consisting of the unordered boundary (pixels of value
%   1) and background (pixels of boundary 0). The boundaries are close by
%   construction, so the first and last entries are the same.
%   For better performance, it might be advisable to choose a boundary
%   which is as "tight" as possible, meaning that it doesn't have a lot of
%   side branches. These can be removed with the appropriate set of
%   morphological operations.
%
%   minlen,maxlen: Thresholds for boundary length. Only select the
%   boundaries that have a length >minlen and <maxlen. to get rid of small
%   occasional circling and nothing else, choose minlen=4 and maxlen=inf.
%
%   b_init: should be set equal to either 1 or 2. Given the way the
%   algorithm works, the same boundary will **in general** be 
%   re-traced in different ways. The first one takes in the longest extent
%   of the boundary, with all of the corner pixels and the possible short
%   side branches. Set b_init to 1 to select the first and 2 to select the
%   second. If in doubt, try! or set directly to 1 if you know that in your
%   case there is only one tracing (which can happen for thin objects
%   without an interior).
%
%OUTPUT:
%   myborders: STRUCTURE with field "pxls". 
%   The element myborders(ii).pxls is a Nx2
%   array, where each row gives the row and col coordinates of a pixel of
%   the "iith" boundary. The boundary pixels are in order. Notice that
%   given the way the algorithm works, the same boundary will in general be
%   re-traced in different ways. The first one takes in the longest extent
%   of the boundary, with all of the corner pixels and the possible short
%   side branches. On top of that, occasionally the side pixels can be
%   locally arranged in such a way pattern the algorithm will also run
%   through them. This happens for example when by chance there is locally
%   a "o" shaped pattern. These are not real objects and we need to remove
%   them. For these reasons, at the end of the code we run a quick pruning
%   of the ordered boundaries. The first gets rid of all of the boundaries
%   of length below minlen and above maxlen (maxlen can be set to "inf").
%   Then, we also select only one of the two "versions" of the boundary of
%   each real feature. This is selected by the parameter "b_init" below.
%   Selecting b_init to be equal to 1 or 2, gives the first or second of
%   the two versions of the boundary.
%
%   If in doubt, you can take a look at the outputted boundaries and make
%   up your mind. Here is a code snippet that might turn out to be useful
%   to inspect the outputted boundaries
% 
% for idx = 1:length(myborders)
% hold on
% imagesc(image);
% set(gcf,'Position',[100 100 1200 800]);
% for ii=1:length(myborders(idx).pxls)
%     plot(myborders(idx).pxls(ii,2),myborders(idx).pxls(ii,1),'r*-','LineWidth',5);
%     pause(0.02)
% end
% hold off
% end
%
%HISTORY:
%   25 March, 2024: MP. Created from the python code from Sebastian
%   Wallkotter. I added a section under "if new_pos=start_pos" to remove
%   duplicate pixel entries. Other differences are due to the absence in
%   Matlab of the advanced array slicing options available in NumPy (at
%   least... I couldn't find their equivalent in Matlab).
%
%   26 March, 2024: MP. Modified earlier version in the following ways: i)
%   the previous version set the displacements tensor to type int8. This
%   causes a problem because when adding "displacements" to "type, row,col"
%   in the main loop, these last variables become of type int8. This is a
%   problem for images with a large number of pixels (max value of int8 is
%   too small); ii) added a length threshold to get rid of spurious ordered
%   boundaries due to accidental local cycles at the boundary of real
%   objects; iii) added a selection of only one of the two versions of the
%   boundary of each object. This can be selected with the parameter b_init
%   below. The best of the two boundary versions is probebly the second one
%   (b_init=2); 
%
%   27 March, 2024: MP. Moved the b_init parameter to be an input.
%
%   4 April, 2024: MP. Corrected typo in the description of the output
%   structure in the function preamble. It said that the output structure
%   had a field called "name", whereas it actually has a field called
%   "pxls".

%% Preparatory stage
%Define the boundary types
NORTH = 1;
EAST = 2;
SOUTH = 3;
WEST = 4;

Npxls = []; %we will store the length of the boundary for the different objects, here


%Generate copies of the original image, translated by one pixel in any of
%the 4 possible directions
padded_img = padarray(img,[1,1],"both");
img_north = padded_img(1:end-2, 2:end-1);
img_south = padded_img(3:end, 2:end-1);
img_east = padded_img(2:end-1, 3:end);
img_west = padded_img(2:end-1, 1:end-2);

%Define the boundary types that can be assigned to each of the pixels. Mind
%that pixels can (and many will) have more than one boundary type.
%Pixels of type "X" are those such that a translation "X" will move them
%out of the boundary.
border = int8(zeros([4,size(padded_img)]));
border(NORTH,2:end-1,2:end-1) = (img==1) & (img_north == 0); 
border(EAST,2:end-1,2:end-1) = (img==1) & (img_east == 0);
border(SOUTH,2:end-1,2:end-1) = (img==1) & (img_south == 0);
border(WEST,2:end-1,2:end-1) = (img==1) & (img_west == 0);

%For pixels of type NORTH, find which adjacent boundary pixels are
%available. In order, choose the first between: 1) pixel of type WEST in
%the NORTH-EAST direction; 2) pixel of type NORTH in the EAST direction;
%pixel of type EAST in the same position (which means that the current
%pixel should also be of type EAST)
adjacent = int8(zeros([4,size(img)]));
[~,adjacent(NORTH,:,:)] = max(cat(3,...
    squeeze(border(WEST,1:end-2, 3:end)),...
    squeeze(border(NORTH,2:end-1, 3:end)),...
    squeeze(border(EAST,2:end-1, 2:end-1))),[],3);

%For pixels of type EAST, find which adjacent boundary pixels are
%available. In order, choose the first between: 1) pixel of type NORTH in
%the SOUTH-EAST direction; 2) pixel of type EAST in the SOUTH direction;
%pixel of type SOUTH in the same position (which means that the current
%pixel should also be of type SOUTH)
[~,adjacent(EAST,:,:)] = max(cat(3,...
    squeeze(border(NORTH,3:end,3:end)),...
    squeeze(border(EAST,3:end, 2:end-1)),...
    squeeze(border(SOUTH,2:end-1, 2:end-1))),[],3);

%For pixels of type SOUTH, find which adjacent boundary pixels are
%available. In order, choose the first between: 1) pixel of type EAST in
%the SOUTH-WEST direction; 2) pixel of type SOUTH in the WEST direction;
%pixel of type NORTH in the same position (which means that the current
%pixel should also be of type NORTH)
[~,adjacent(SOUTH,:,:)] = max(cat(3,...
    squeeze(border(EAST,3:end,1:end-2)),...
    squeeze(border(SOUTH,2:end-1, 1:end-2)),...
    squeeze(border(WEST,2:end-1, 2:end-1))),[],3);

%For pixels of type WEST, find which adjacent boundary pixels are
%available. In order, choose the first between: 1) pixel of type SOUTH in
%the NORTH-WEST direction; 2) pixel of type WEST in the NORTH direction;
%pixel of type NORTH in the same position (which means that the current
%pixel should also be of type NORTH)
[~,adjacent(WEST,:,:)] = max(cat(3,...
    squeeze(border(SOUTH,1:end-2,1:end-2)),...
    squeeze(border(WEST,1:end-2,2:end-1)),...
    squeeze(border(NORTH,2:end-1,2:end-1))),[],3);

%combine an array of moves. For example, if you are a pixel NORTH (i.e. type 1) and you
%do a move of type "1", then you will change your type to WEST (i.e. type 1
%plus 3), change your row by "-1" and change you column by "+1". 
directions = int32(zeros([4,size(img),3,3]));
for rr = 1:size(img,1)
    for cc = 1:size(img,2)
        directions(NORTH,rr,cc,:,:) = [[3, -1, 1];[0, 0, 1];[1, 0, 0]];
        directions(EAST,rr,cc,:,:) = [[-1, 1, 1];[0, 1, 0];[1, 0, 0]];
        directions(SOUTH,rr,cc,:,:) = [[-1, 1, -1];[0, 0, -1];[1, 0, 0]];
        directions(WEST,rr,cc,:,:) = [[-1, -1, -1];[0, -1, 0];[-3, 0, 0]];
    end
end


unprocessed_border = border(:,2:end-1,2:end-1); %array containing all of the border pixels for all of the different types.
sz_up = size(unprocessed_border);
myborders = struct; %here we will store the ordered borders
idx = 0; %running index of the border (border id)
allborders = find(unprocessed_border~=0); %list of all of the border pixels. Mind that if a certain pixel is both, say, a type NORTH pixel and a type SOUTH pixel, then it will appear twice 

%% Main loop to walk along the borders

for pp = 1:length(allborders)
    start_pos = allborders(pp); %choose one of the pixels at the border (remember this include the type)
    if  unprocessed_border(start_pos) %has this pixel already been used? If not, go ahead
        idx = idx+1; %we will be tracing boundary idx+1
        current_pos = start_pos; %define the current position
        [type,row,col]=ind2sub(sz_up,current_pos); %triplet of pixel type, pixel row and pixel column for the current pixel
        myborders(idx).pxls=[row,col]; %append current position to the list of ordered pixels of the current boundary
        while 1
            unprocessed_border(current_pos)=0; %set the current position in the unprocessed_border to zero (we are processing it now, so we remove it from the list of unprocessed border pixels)
            displacement = squeeze(directions(type,row,col,adjacent(type,row,col),:)); %select the move to the next pixel
            type = type + displacement(1); %update the type of the new pixel
            row = row + displacement(2); %update the row of the new pixel 
            col = col + displacement(3); %update the col of the new pixel
            myborders(idx).pxls=[myborders(idx).pxls;[row,col]]; %append the new pixel position to the ordered list of pixels
            new_pos = sub2ind(sz_up,type,row,col); %array index for the new position 
            if new_pos == start_pos %have I come back to the beginning yet? If so, do the following
                ww = any(myborders(idx).pxls(2:end,:)-myborders(idx).pxls(1:end-1,:),2); %check where we have back to back repeats
                myborders(idx).pxls=[myborders(idx).pxls(ww,:);myborders(idx).pxls(end,:)]; %remove the repeats and append the final position at the end of the ordered list of boundary pixels
                Npxls(idx) = length(ww);
                break %exit
            else
                current_pos = new_pos; %if we have not gotten back to the beginning, set the current position to be equal to the new position and repeat
            end
        end
    end

end

%select only the boundaries whose length is larger than minlen and smaller
%than maxlen
ww=Npxls>minlen & Npxls<maxlen;
myborders = myborders(ww);

%At this point, each "real" boundary should be providing two entries in
%"myboundaries": one where the boundary is cycled clockwise and the other
%where it is cycled counterclockwise (notice that this direction depends on how you plot the image). 
%We will select the second of the two options, which is usually the one
%that has the least wiggles caused by corner boundary pixels or short side
%branches.

myborders = myborders(b_init:2:end);

end




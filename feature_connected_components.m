function [features,centroids] = feature_connected_components(BWimg,param)
%Function to calculate features of the connected components of a BW image.
%It first determines the connected components and calculates a few region
%properties (more can be added easily). Then, for each feature, it traces
%the boundary and calculates a closed spline fit to it, together with the
%tangent angle and local curvature. The pp-forms for contour, tangent angle
%and curvature are included in the output as well.
%
%INPUT:
%   BWimg : binary image where the features correspond to regions in the
%   image, of value 1.
%   param: structure with the following fields, 
%       >param.minlen;param.maxlen: Thresholds for boundary length. Only
%       select the boundaries that have a length >minlen and <maxlen. to
%       get rid of small occasional circling and nothing else, choose
%       minlen=4 and maxlen=inf. 
%       param.b_init: see preamble to trace_boundary.
%
%OUTPUT: 
%
%   features : structure of length Number-of-features and including the
%   following fields:
%       >Centroid : remember that for this, the first coordinate is the
%       x-position (COLUMN) and the second is the y-position (ROW). This
%       can cause some confusion as it is the OPPOSITE of the ordering for
%       the boundary pixels and the boundary spline fit (these last are
%       ordered as (row,col), following standard matrix convention).
%       >Area, Solidity, BoundingBox: as per specifications in the Matlab
%       documentation
%       >pxlborder : ordered array of border pixels. Notice that the
%       starting point is not predetermined.
%       >outs,outpos,outtheta,outkappa,pp_pos,pp_theta,pp_kappa : as per
%       specifications in the preamble to the function "loop_shape.m"
%
%HISTORY:
%   04 April, 2024: MP. Created
%   05 April, 2024: MP. Modified input to group parameters in a single
%   "param" structure. 
%
%
%TODO:
%MP: It Would be good to include also the quadrupolar moment of the component
%(see also Bohr deformation parameters; Zernike polynomials, see niu22.pdf
%( Kuo Niu and Chao Tian 2022 J. Opt. 24 123001). See also bookstein96 on
%Procrustes and morphometrics. There is a Matlab implementation of
%Procrustes. Look up morphometrics in yeast.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = bwconncomp(BWimg);
features = regionprops(BWimg,"Centroid","Area","Solidity","BoundingBox");

%% for each label, trace the boundary and then run loop shape. add the results to the structure returned from regionprops 
for nn=1:CC.NumObjects
    col = floor(features(nn).BoundingBox(1));
    row = floor(features(nn).BoundingBox(2));
    deltacol = features(nn).BoundingBox(3);
    deltarow = features(nn).BoundingBox(4);
    myborders=trace_boundary(BWimg(row:row+deltarow,col:col+deltacol),param.minlen,param.maxlen,param.b_init);
    features(nn).pxlborder = [myborders(1).pxls(:,1)+(row-1),myborders(1).pxls(:,2)+(col-1)];
    [outs,outpos,outtheta,outkappa,pp_pos,pp_theta,pp_kappa]=loop_shape(features(nn).pxlborder);
    features(nn).outs = outs;
    features(nn).outpos = outpos;
    features(nn).outtheta = outtheta;
    features(nn).outkappa = outkappa;
    features(nn).pp_pos = pp_pos;
    features(nn).pp_theta = pp_theta;
    features(nn).pp_kappa = pp_kappa;
end

end








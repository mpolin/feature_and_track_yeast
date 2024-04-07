function [outs,outpos,outtheta,outkappa,pp_pos,pp_theta,pp_kappa]=filament_shape(rpos)
%This function smoothes the closed contour given by the positions in the
%input variable rpos, resamples the smoothed version of the closed contour
%into a fixed number of positions, outputs the arclength of these points,
%the local tangent angle and the local curvature, together with their
%pp_forms.
%The initial (raw) positions along the filament are first smoothed with a
%cubic spline (smoothing rsmooth). This results in a smoothed open line,
%where the beginning and the end, which should be "smoothly" joined, are
%in fact open. To remedy this, we take the smoothed version, remove points
%1 and "last" and say that the point right after "last-1" is given by point
%2. This is done to remove the location where enforcing a periodic condition
%on the spline would potentially result in spurious tight bends. The result
%is then pline FITTED with periodic boundary conditions. The resulting
%spline is then evaluated at "nsamples" equally spaced samples and is also
%used to calculate the local tangent angle and the local curvature.
%"nsamples" is currently hardwired to 200 but it can easily be moved to an
%input. 
%
%INPUT
%
%   rpos = (npoints x 2) matrix of raw positions along the filament.
%   Columns correspond to x and y coordinates respectively.
%
%OUTPUT
%
%   outs = arclength of the nsamples resampled points along the smoothed
%   filament.
%   outpos = (nsamples x 2) matrix of resampled positions along the
%   filament. Columns correspond to x and y coordinates respectively.
%   outtheta = angles of the local tangent vectors to the filament, at the
%   resampled points positions. 
%   outkappa = curvatures of the local tangent vectors to the filament, at
%   the resampled points positions. The unit of measurement of this
%   curvature is the inverse of the unit of measurement of the positions in
%   the input matrix rpos.
%   pp_pos, pp_theta, pp_kappa = pp-forms for the periodic, smoothed
%   boundary, local tangent angle, local curvature
%
%HISTORY
%
%   28 March 2024: Marco Polin. Created from filament_shape.m


nsamples = 200; %we will resample the filament at 200 equally spaced points
rpos = double(rpos); %just to make sure that we're working with double 
rpos = rpos(1:end-2,:); %remove the repeating endpoint
rs = mycontourlength(rpos); %calculate the approximate contour length for the different points
rsmooth = 1/(1+mean(rs(2:end)-rs(1:end-1))^3/0.06); %smoothing used to iron out small fluctuations in the shape. This seems to be working fine but it needs to be checked further

%% 1) smooth the shape of the filament with a cubic spline fit 
pp_pos = csaps(rs,rpos',rsmooth); %pp form for the cubic spline FIT
tmppos = (fnval(pp_pos,rs))'; %evaluate the new positions
tmppos = [tmppos(2:end-1,:);tmppos(2,:)]; %remove first and last point and append the second point as last. This removes the points immediately close to the "cut" where the initial and final bits of the curve try to meet
tmps = mycontourlength(tmppos); %recalculate the contour length (for good measure)
%at this point we have a smoothed version of the contour, which is closed.

%% 2) use the smoothed points to do a spline INTERPOLATION with periodic boundary conditions
pp_pos = csape(tmps,tmppos','periodic'); %pp form for the cubic spline INTERPOLATION with periodic constraint
tmppos = (fnval(pp_pos,tmps))';%evaluate the new positions
tmps = mycontourlength(tmppos);%recalculate the contour length. This is used to have a better idea of what is the length of the contour

outs = (linspace(0,tmps(end),nsamples))'; %calculate the set of contour lengths at which we want to know the boundary position 
outpos = (fnval(pp_pos,outs))'; %calculate the boundary positions
outs = mycontourlength(outpos); %recalculate the corresponding contour lengths

%% get the local tangent angle and local curvature
tt_pos=fnder(pp_pos); %derivative of the positions spline pp-form
tt = (fnval(tt_pos,outs))'; %these are the local tangent verctors (not normalised)
tt = tt./abs(tt(:,1)+i*tt(:,2));%these are the local tangent verctors (normalised)... in case they might be needed 
outtheta = unwrap(angle(tt(:,1)+i*tt(:,2))); %theta(s) derived directly from the angle of the tangent vectors

pp_theta = csape(outs,outtheta); %gives the ppform for the cubic spline interpolation to out_theta(s). This is calculated to be able to use fnder to perform the derivative and obtain the local curvature 
pp_kappa = fnder(pp_theta); %calculates the derivative. This will be used to evaluate the curvature kappa
outkappa = fnval(pp_kappa,outs);

end

function s = mycontourlength(pos)

%function to calculate the contour length of the curve described by the
%points in pos. Currently it is based only on the distances between
%successive points, but it can be refined (look for arclength estimate in
%matlab)

disp = pos(2:end,:)-pos(1:end-1,:); %point-to-point displacements
ds = sqrt(disp(:,1).^2+disp(:,2).^2); %magnitude of displacements
s = [0;cumsum(ds)]; %contour length of the differenc points for original (raw) sampling. This can be refined by a better evaluation of the contour length from the raw data. See function "actlength" in Matlab File Exchange

end


# Featuring and Tracking yeast cells

First set of analysis of the experiments on yeast ascus digestion. Run "do_features_extraction.m" to:
  - For each image:
      > segment the features in the image 
      > calculate parameters for each feature like centroid position, area, etc.
      > trace the boundary and return a spline fit for the boundary shape (together with other goodies, check functions' preambles)
  - track the features across the frames
  - save everything

## TO DO: 
- improve segmentation
- track the changes in orientation of the objects
- write "pull" functions to pull specific sets of observables for a specific object, along the whole movie.


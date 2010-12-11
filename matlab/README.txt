This is a Matlab implementation of
	Crandall et al. "Spatial Priors for Part-based Recognition using Statistical Models", CVPR'05

The general steps to run the program involve:
1. generate edge maps from background and object images
2. generate object model (using a k-fan for the shape and the edge maps for the appereance)
3. try to detect object in test image


1. GENERATE EDGE MAPS
==========================================================================================================
Use runEdgeInFolder to compute edge maps. For example,
>> runEdgeInFolder('background_images','background_edges')

2. GENERATE OBJECT MODEL
==========================================================================================================

2.1 Appereance model:







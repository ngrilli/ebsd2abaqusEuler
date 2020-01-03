% Nicolo Grilli
% University of Oxford
% AWE project 2019
% 21 November 2019

% Matlab script using MTEX and EBSD2abaqus
% to generate input file from EBSD map for Titanium
% C:\Users\engs1992\YiXiong\EBSD\Ti6Al-before.ctf

% open MTEX and launch
% Import EBSD data
% import the file:
% C:\Users\engs1992\YiXiong\EBSD\Ti6Al-before.ctf
% x axis should be set horizontal towards the right
% y axis should be set vertical towards the bottom
% then import the variable ebsd in the workspace

ipfKey = ipfHSVKey(ebsd('Titanium').CS); % get color key
color = ipfKey.orientation2color(ebsd('Titanium').orientations); % get color map

% data filtering to remove non-indexed points
% However, if large regions with non-indexed locations exist
% then this technique is slow because F.numNeighbours
% must be set to a large number (> 20)
F = medianFilter;
F.numNeighbours = 10;
ebsd_smoothed = smooth(ebsd,F);

% but you can use also the fill function
% first you find a vector where each component indicates if a particular phase
% has been indexed or not
% for instance in this case it will be:
% [1 0]
% because first phase is not indexed and second phase is Titanium
hasNotIndexed = strcmp(ebsd.mineralList,'notIndexed') & strcmp(ebsd.mineralList,'notIndexed');
% the ebsd associated with the list of phases that are not indexed is found
ebsdNotInd = ebsd(ebsd.mineralList(hasNotIndexed));
% then remove the corresponding EBSD measurement
ebsd(ebsd.mineralList(hasNotIndexed)) = [];
% then fill the EBSD map by
% extrapolating spatial EBSD data by nearest neighbour for tetragonal lattice
% if fill is used without removing the non-indexed measurements
% many non-indexed points remain
ebsd = fill(ebsd);

% go to the EBSD2abaqus folder
% and run the function to remove the smaller grains
% 5 degrees misorientation is sufficient to find the same
% grain structure of grains 1 to 7 that is reported in
% C:\Users\engs1992\YiXiong\groupmeeting1209.pptx
grainAngle = 5.0;
ebsd_clean = clean4fem(ebsd,500,0.1,grainAngle);

% uncomment to plot the EBSD map
% plot(ebsd_clean,ebsd_clean.orientation)

% reduce the pixel size of the EBSD map
ebsd_reduced = reduce(ebsd_clean,3); % take every third pixel horiz. and vert.

% build the abaqus input file
% with the same misorientation at 5 degrees
% However, it will be necessary to reduce the number of elements
% for running crystal plasticity simulations
ebsd2abaqusEuler(ebsd_reduced,5);

% define grains with the same angle criterion as before
[grains,ebsd_reduced.grainId,ebsd_reduced.mis2mean] = calcGrains(ebsd_reduced,'angle',grainAngle*degree);

% for instance the Euler angles (in radians) of grain 1 are found as:
grains(1).meanOrientation.phi1
grains(1).meanOrientation.Phi
grains(1).meanOrientation.phi2






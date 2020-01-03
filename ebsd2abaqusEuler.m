function ebsd2abaqusEuler(ebsd,angle)
    % *ebsd2abaqus* generates ABAQUS input file for given EBSD map
    % 
    % % Syntax
    % ebsd2abaqus(ebsd,angle)
    %
    % % Input
    % ebsd  - ebsd data
    % angle - threshold angle for grain reconstruction
    
    % % Output
    % ebsd.inp file which contains
    % - element sets with individual grains
    % - element sets with individual phases
    % - sections with individual phases
    % - node sets of faces for BCs

    % get Euler angles of grains
	% and write rotation matrices to abaqus input file
	% to use with the crystal plasticity UMAT
    
    % --------------------------
    % written by
    % Marat I. Latypov (GT Lorraine) 
    % marat.latypov@georgiatech-metz.fr
    % March 2015, revised July 2016
    % --------------------------
	% Nicolo Grilli
	% University of Oxford
	% AWE project 2019
	% --------------------------
    
    %% preliminarities
    hasNotIndexed = strcmp(ebsd.mineralList,'notIndexed') & strcmp(ebsd.mineralList,'notIndexed');
    % check if ebsd on hex grid and convert to sqr if so
    if length(ebsd.unitCell) == 6
        ebsd = fill(ebsd);
        fprintf('WARNING! EBSD was on hex grid and so was converted to sqr grid using fill function\n')
        figure; plot(ebsd)
    % check if ebsd on hex grid and convert to sqr if so
    elseif any(hasNotIndexed) == 1
        ebsdNotInd = ebsd(ebsd.mineralList(hasNotIndexed));
        ebsd(ebsd.mineralList(hasNotIndexed)) = [];
        ebsd = fill(ebsd);
        numNotIndPx = numel(ebsdNotInd.x);
        fprintf('WARNING! EBSD had %d non-indexed pixels and so was filled using fill function\n', numNotIndPx)
        figure; plot(ebsd)
    end
    
    % get step size from input
    dxy = max(ebsd.unitCell) - min(ebsd.unitCell);
    step(1) = dxy(1);
    step(2) = dxy(2);
    step(3) = min(dxy); 
    % set inp file name
    inpFileName = 'ebsd.inp';
    
    % reconstruct grains
    [grainsReconstructed,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',angle*degree);
    phases = ebsd.phase;
    grains = ebsd.grainId; % grain index for each point of the EBSD
    
    % get coordinates
    xyz = zeros(numel(ebsd.x),3);
    xyz(:,1) = ebsd.x;
    xyz(:,2) = ebsd.y;
    xyz(:,3) = zeros(numel(xyz(:,1)),1);
    
    % sort the coordinate array
    [~,order] = sortrows(xyz,[1,3,2]);
    
    % reorder grain IDs and Euler angles according to ABAQUS convention
	% but same grain ID correspond to same Euler angles
	% so no need to reorder grainsReconstructed
    phases = phases(order);
    grains = grains(order);
    
    % get the number of voxels along x, y, z
    xVox = size(unique(xyz(:,1)),1);
    yVox = size(unique(xyz(:,2)),1);
    zVox = size(unique(xyz(:,3)),1);

    % get step size and boundaries for the mesh
    boxmin = zeros(1,3);
    boxmax = zeros(1,3);
    for ii = 1:3
        boxmin(ii) = min(xyz(:,ii))-step(ii)/2;
        boxmax(ii) = max(xyz(end,ii))+step(ii)/2;
    end

    %% Generate 3D mesh
    % generate nodes 
    [x,y,z] = meshgrid(boxmin(1):step(1):boxmax(1),boxmin(2):step(2):boxmax(2),boxmin(3):step(3):boxmax(3));
    numNodes = numel(x);
    coord = [reshape(x,numNodes,1), reshape(y,numNodes,1), reshape(z,numNodes,1)];
    nodes = [(1:numNodes)', sortrows(coord,[1,3,2])];

    % allocate array for elements
    elem = zeros(size(xyz,1),9);
    count = 1;

    % start loop over voxel dimensions
    for ix = 1:xVox
        for iz = 1:zVox
            for iy = 1:yVox

                % get element label
                elem(count,1) = count;

                % nodes on the plane with lower x
                elem(count,2) = iy + (iz-1)*(yVox+1) + (ix-1)*(yVox+1)*(zVox+1);
                elem(count,3) = elem(count,2) + 1;
                elem(count,4) = elem(count,3) + yVox + 1;
                elem(count,5) = elem(count,2) + yVox + 1;

                % nodes on the plane with higher x
                elem(count,6) = iy + (iz-1)*(yVox+1) + ix*(yVox+1)*(zVox+1);
                elem(count,7) = elem(count,6) + 1;
                elem(count,8) = elem(count,7) + yVox + 1;
                elem(count,9) = elem(count,6) + yVox + 1;

                count = count+1;
            end
        end
    end

    %% Write inp file
    % open inp file and write keywords 
    inpFile = fopen(inpFileName,'wt');
    fprintf(inpFile,'**PARTS\n**\n');
    fprintf(inpFile,'*Part, name=SAMPLE\n');

    % write nodes
    fprintf(inpFile,'*NODE, NSET=AllNodes\n');
    fprintf(inpFile,'%d,\t%e,\t%e, \t%e\n',nodes');

    % write elements
    fprintf(inpFile,'*Element, type=C3D8, ELSET=AllElements\n');
    fprintf(inpFile,'%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',elem');
    
    % create element sets containing grains
    for ii = 1:numel(unique(grains))
        fprintf(inpFile,'\n*Elset, elset=Grain-%d\n',ii);
        fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',elem(grains==ii)');
    end   

    % create element sets containing phases (no need if material is single phase)
    %uniPhases = unique(phases);
    %phaseNames = ebsd.mineralList;
    %phaseNames = phaseNames(~strcmp(phaseNames,'notIndexed') & ~strcmp(phaseNames,'notIndexed'));
    %for ii = 1:numel(unique(phases))
    %    fprintf(inpFile,'\n*Elset, elset=Phase-%s\n',phaseNames{ii});
    %    fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',elem(phases==uniPhases(ii))');
    %end
    
    % create sections for phases (no need if material is single phase)
    %for ii = 1:numel(uniPhases)
    %    fprintf(inpFile,'\n**Section: Section-%s\n*Solid Section, elset=Phase-%s, material=%s\n',phaseNames{ii},phaseNames{ii},phaseNames{ii});
    %end
    %fprintf(inpFile,'\n');
	
	% create sections for grains
	for ii = 1:numel(unique(grains))
		fprintf(inpFile,'\n**Section: Section-%d\n*Solid Section, elset=Grain-%d, material=Grain-%d\n',ii,ii,ii);
	end

    % create node sets containing surface nodes for BCs
    for ii = 1:3
        fprintf(inpFile,'\n**\n*Nset, nset=NODES-%d\n',ii);
        fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',nodes(nodes(:,ii+1)==boxmin(ii))');
        fprintf(inpFile,'\n**\n*Nset, nset=NODES+%d\n',ii);
        fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',nodes(nodes(:,ii+1)==boxmax(ii))');
    end
	
	% write closing 
    fprintf(inpFile,'\n**\n*End Part\n');
	
	% get Euler angles of the grains
	% and write them to input file in materials
	fprintf(inpFile,'\n**\n** MATERIALS\n**');
	for ii = 1:numel(unique(grains))
		tempphi1 = grainsReconstructed(ii).meanOrientation.phi1;
		tempPhi = grainsReconstructed(ii).meanOrientation.Phi;
		tempphi2 = grainsReconstructed(ii).meanOrientation.phi2;
		tempRotMat = bungeMatrix(tempphi1,tempPhi,tempphi2);
		fprintf(inpFile,'\n*Material, name=Grain-%d',ii);
		fprintf(inpFile,'\n*Depvar');
		fprintf(inpFile,'\n    125,'); % set proper number of state variables
		fprintf(inpFile,'\n*User Material, constants=10');
		fprintf(inpFile,'\n0.,%6.5f,%6.5f,%6.5f,%6.5f,%6.5f,%6.5f,%6.5f,',tempRotMat(1,1),tempRotMat(1,2),tempRotMat(1,3),tempRotMat(2,1),tempRotMat(2,2),tempRotMat(2,3),tempRotMat(3,1));
		fprintf(inpFile,'\n%6.5f,%6.5f',tempRotMat(3,2),tempRotMat(3,3));
	end
    
    fprintf('Success: input file written as %s\n', inpFileName);
    
    % close the file
    fclose(inpFile); 
end

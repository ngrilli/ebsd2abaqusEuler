function [ret,indSmallCi] = clean4fem(ebsd,minSize,minCi,angle)
    % *clean4fem* removes data from inaccurate pixels and
    % assigns them with the data from the surrounding grain
    %
    % pixels are considered inaccurate if 
    % - their CI is lower than allowed minimum, or 
    % - they belong to grain whose px size is lower than allowed minimum
    %
    % % Syntax
    % clean4fem(ebsd,minSize,minCi,angle)
    %
    % % Input
    % ebsd    - raw ebsd data to be cleaned up (required)
    % minSize - minimum allowed grain size [in px] (default: 25 px)
    % minCi   - minimum allowed confidence index (default: 0.1)
    % angle   - threshold angle for grain reconstruction (default: 15 deg)

    % --------------------------
    % written by
    % Marat I. Latypov (GT Lorraine) 
    % marat.latypov@georgiatech-metz.fr
    % March 2015
    % --------------------------

    %% remove inaccurate pixels
    % segment grains 
    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',angle*degree);

    figure; plot(ebsd)
    set(gcf,'renderer','zBuffer')

    % remove non-indexed pixels
    cleanEbsd = ebsd;
    cleanEbsd('notIndexed') = [];

    % check if CI or MAD properties exist
    try 
        indSmallCi = cleanEbsd.ci <= minCi;
        catch ME
        if strcmp(ME.identifier,'MATLAB:noSuchMethodOrField')
            try
                indSmallCi = cleanEbsd.confidenceindex <= minCi;
            catch ME
                if strcmp(ME.identifier,'MATLAB:noSuchMethodOrField')
                    try
                        indSmallCi = cleanEbsd.mad <= minCi;
                    catch ME 
                        if strcmp(ME.identifier,'MATLAB:noSuchMethodOrField')
                            fprintf(1,'WARNING: no CI property found! Cleaning of pixels CI will be skipped');
                            indSmallCi = false(numel(ebsd.x),1);
                        end 
                    end 
                end
            end
        end
    end
    
    % remove pixels with low CI
    cleanEbsd(indSmallCi) = [];
    
    % remove pixels of tiny grains
    indSmallSize = grains.grainSize < minSize;
    cleanEbsd(grains(indSmallSize)) = [];

    % plot data without inaccurate points
    figure; plot(cleanEbsd)
    set(gcf,'renderer','zBuffer')

    % identify grains with 15 degree boundary:
    [cleanGrains,cleanEbsd.grainId,cleanEbsd.mis2mean] = calcGrains(cleanEbsd,'angle',angle*degree);

    %% assign good values
    poly = cleanGrains.poly;

    for ii = 1:max(cleanEbsd.grainId)

        v = cleanGrains.V(poly{ii},:);
        inside = inpolygon(ebsd.x,ebsd.y,v(:,1),v(:,2));
        ebsd(inside).grainId = cleanGrains(ii).id;
        ebsd(inside).phase = cleanGrains(ii).phase;
        ebsd(inside).rotations = cleanGrains(ii).meanOrientation;
    end
    
    figure; plot(ebsd)
    set(gcf,'renderer','zBuffer')
    ret = ebsd;
end

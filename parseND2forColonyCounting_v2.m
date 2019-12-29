function [] = parseND2forColonyCounting_v2(scanFile, scanDim, varargin)
    p = inputParser;

    p.addRequired('scanFile', @ischar);
    p.addRequired('scanDim',  @(x)validateattributes(x,{'numeric'},{'size',[1 2]}));
    
    p.addParameter('inDir', '', @ischar);
    p.addParameter('outDir', '', @ischar);
    p.addParameter('divideX', 1, @isnumeric);
    p.addParameter('divideY', 1, @isnumeric);
    p.addParameter('resizeFactor', 1, @isnumeric)
    
    p.parse(scanFile, scanDim, varargin{:});
    
    scanFile = p.Results.scanFile;
    scanDim = p.Results.scanDim;
    divideX = p.Results.divideX;
    divideY = p.Results.divideY;
    resizeFactor = p.Results.resizeFactor;
    
    if ~isempty(p.Results.inDir)
        inDir = p.Results.inDir;
    else
        inDir = pwd;
    end
    
    if ~isempty(p.Results.outDir)
        outDir = p.Results.outDir;
    else
        outDir = fullfile(inDir, 'splitPoints');
    end

%     %For testing
%       scanDim = [40 100];
%       divideX = 1;
%       divideY = 2;
    
    % Make directory for  images
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end
    
    % Read scan file
    reader = bfGetReader(fullfile(inDir, scanFile)); 

    omeMeta = reader.getMetadataStore();
    wavelengths = omeMeta.getPixelsSizeC(0).getValue(); % number of wavelength channels
    
    %Create matrix of how the scan was acquired. This can be updated as a
    %parameter in the future for other scan patterns. 
    scanMatrix = vec2mat(1:scanDim(1)*scanDim(2), scanDim(2));
    for i = 2:2:scanDim(1)
        scanMatrix(i, :) = fliplr(scanMatrix(i, :));
    end
    
    nRegions = divideX * divideY;
    regions = cell(1,nRegions);
    regionDim = idivide(int16(scanDim),int16([divideX divideY]), 'floor');
    
    for i = 1:divideX
        for ii = 1:divideY
            regionTiles = scanMatrix(1+((i-1)*regionDim(1)):((i)*regionDim(1)), 1+((ii-1)*regionDim(2)):((ii)*regionDim(2)));
            regions{i+((ii-1)*divideX)} = regionTiles(:);
            %regions{i+((ii-1)*divideX)} = regionTiles;
        end
    end
    
    for i = 1:numel(regions)
        for ii = 1:numel(regions{i})
            for iii = 1:wavelengths
                reader.setSeries(regions{i}(ii)-1);
                iPlane = reader.getIndex(0, iii - 1, 0) + 1;
                tmpPlane  = bfGetPlane(reader, iPlane);
                tmpPlane = imresize(tmpPlane, 1/resizeFactor);

                imwrite(tmpPlane, fullfile(outDir, sprintf('Scan%03d_w%d_s%d_t1.TIF', i, iii, ii)))
            end
        end
        fprintf('Finished parsing region %d\n', i);
    end
end
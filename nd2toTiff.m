function [] = nd2toTiff(infile,varargin)
% Reads .nd2 files and produces .tif files
% Creates a .tif for each wavelenght and z stack 
%
% Input
% *|infile|*  - filename or full filepath string of the image stack
% Optional:
%          *|outDir|*  - filepath for the Outputn (Current Dirctory)
%          *|nDigits|* - max order of magnitude for naming convention 
%                        ei. 'dapi0001'  nDigits = 4 (3 )
% 
% Output
% *|tiff|* - 3D tiff files named "NAME#" were # coresponds to the Z stack and
%           takes into acount previous files. NAME is taken from the wavelength
%           channel (ie DAPI, CY3,...) and changed so that arjlabimagetools can
%           read it.
%           This is an Example of the naming convention as follows:
%                   Our_names     arj_names
%                   alexa594  ->  alexa
%                   Atto647N  ->  cy
%                   cy3       ->  tmr
%                   700       ->  nir
%          Add to this list in the "channelMap"
%
% Required:
%            Get bfmatlab   
%                    1)Go to:   https://downloads.openmicroscopy.org/bio-formats/6.0.1/artifacts/
%                    2)Download bfmatlab.zip
%                    3)Unzip and move bfmatlab folder to your MATLAB folder
%                    4)Add bfmatlab path to Matlab     
% Usage
%  >> T1 = nd2toTiff('tmr001.nd2');               % read image stack in working directory
%
%  >> T1 = nd2toTiff('/Path/to/file/tmr001.nd2'); % read from specific filepath
%
%  >> T1 = nd2toTiff(infile,"outDir",'/Path/to/output/folder'); % output files in specified folder
%                                                                                    
%  >> T1 = nd2toTiff('tmr001.nd2'); % read from specific filepath
%
% Shortcuts
%  >> T1 = nd2toTiff('/Path/to/file/*.nd2');    % read all nd2 files from specified filepath
%
%  >> T1 = nd2toTiff(''/Path/to/folder/');      % read all nd2 files from specified folder
%
%  >> T1 = readmm('tmr001');                    % file in the current directory
%  
% Last Update 05/15/2019



% The MAP
channelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'GFP', 'CY3', 'A594', 'CY5', 'A647', '700', 'CY7','NIR'},...
                            {'trans'      , 'dapi', 'gfp', 'gfp', 'tmr', 'alexa', 'cy', 'cy'  , 'nir', 'nir','nir'});


% Input check
p = inputParser;
p.addRequired('inDir', @ischar)
p.addParameter('outDir', '', @ischar);
p.addParameter('nDigits', 3, @(x)validateattributes(x,{'numeric'},{'positive','integer'}));

p.parse(infile, varargin{:});



%  ------------------------------------------------------------------------
%                   What Dir to use????  not the nicest but it works
%  ------------------------------------------------------------------------
 infile = p.Results.inDir;

 % Get Info from infile
[Dir,file_name,c] = fileparts(infile);

if file_name == "*" % for the *.nd2 case
    c = [];
end

if numel(c) == 0    % Folder path given
    Files = ls(Dir);
    file_name_nd2 = string(split(Files));
    TF = contains(file_name_nd2,'.nd2');
    file_name_nd2 = file_name_nd2(TF);

    file_name = strings(numel(file_name_nd2));
    for j = 1:numel(file_name_nd2)
        [~,hold,~] = fileparts(file_name_nd2(j));
        file_name(j) = hold;
    end

    fprintf('Reading %d Files:  \n',numel(file_name_nd2));
else % File path given
    file_name = string(file_name);
    file_name_nd2 = string(strcat(file_name,'.nd2'));
    fprintf('Reading one File:  \n');
end      


if isempty(p.Results.outDir)       
    outDir = pwd;
else
    outDir = p.Results.outDir;
end



%  ------------------------------------------------------------------------
%                  Are there any previous files  
%  ------------------------------------------------------------------------

% Assign starting number for output files
if isempty(dir(fullfile(outDir, '*.tif')))
     imageCount = 0;
else
    currentFiles = dir(fullfile(outDir, '*.tif'));
    currentFiles = {currentFiles.name};
    currentFileCount = regexp(currentFiles, '\d+', 'match');
    currentFileCount = vertcat(currentFileCount{:});
    currentFileCount = str2num(cell2mat(currentFileCount));
    imageCount = max(currentFileCount); 
end 

%  ------------------------------------------------------------------------
%                          Read Through nd2  
%  ------------------------------------------------------------------------
nDigits = num2str(p.Results.nDigits);
cnt_mult_files = imageCount;
cnt_stacks = 0;
for f = 1:numel(file_name_nd2)
    
    if mod(f-1,5) == 0
        % file is being read
        fprintf('                 %s\n',file_name(f));
    end
    
    % Read through 
    full_Data_path = fullfile(Dir,file_name_nd2(f));
    reader = bfGetReader(char(full_Data_path));
    omeMeta = reader.getMetadataStore();
    
    Zstacknumb = omeMeta.getImageCount();                    % number of Z stacks
    for i = 1:Zstacknumb % Running through # of z stacks
        cnt_stacks = cnt_stacks + 1;
        
        reader.setSeries(i-1)                                % set ith z stack 
        
        if mod(i,25) == 0
            % stack being read
            fprintf('                       Stack = %03d\n',i);
        end
        
        % Getting info from .nd2 file 
        stackSizeC = omeMeta.getPixelsSizeC(i-1).getValue(); % # of wavelength channels
        stackSizeZ = omeMeta.getPixelsSizeZ(i-1).getValue(); % # of Z slices

        % Arrays
        wave_numb = repmat(1:stackSizeC,1,stackSizeZ);     % repeating vector of channel #
        Z_numb    = repmat(1:stackSizeZ,1,stackSizeC);     % repeating vector of Z #

        for ii = 1:stackSizeC  % (# channels in i)
            cnt = 1;
            for iii = 1:stackSizeZ  % (# z slices in stack)
                
                % Read plane from series iSeries at Z, C, T coordinates (iZ, iC, iT)
                iPlane = reader.getIndex(Z_numb(iii) - 1, wave_numb(ii) - 1, 0) + 1; 
                stack_fig  = bfGetPlane(reader, iPlane);

                %---------------------------------------------------------
                %                         Save Tiff 
                %---------------------------------------------------------
                % Change name 
                channelName = omeMeta.getChannelName(i-1, ii-1);
                channelName = channelName.toCharArray';
 
                outBaseName = strcat('%s%0', nDigits, 'd.tif');
                outputFileName = fullfile(outDir, sprintf(outBaseName,channelMap(channelName),cnt_mult_files + cnt_stacks));

                if cnt == 1
                    imwrite(stack_fig, outputFileName)
                    cnt = 1 + cnt;
                else
                    imwrite(stack_fig, outputFileName, 'writemode', 'append')
                    cnt = 1 + cnt;
                end
                
            end
        end            
     end
end
end
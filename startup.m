%%% Originally written by prof S. Lee

% % Reset path
% temporarily turn off warning; if warning is not turned off above, MATLAB
% might give warning
warning('off','MATLAB:dispatcher:pathWarning');
restoredefaultpath; % reset path to default
warning('on','MATLAB:dispatcher:pathWarning'); % turn on warning back

% % Add to path
fpath = fileparts(mfilename('fullpath')); % the TN directory in which this "startup.m" lies
addpath(fpath);
dirnames = {'linearoptics', 'tensor', 'utils'}; % sub-directories
for it1 = (1:numel(dirnames))
    fpath2 = [fpath,filesep,dirnames{it1}];
    if exist(fpath2,'dir')
        addpath(genpath(fpath2)); % add sub-directories to path
    end
    addpath(genpath(fpath2));
end

% % Clear variables in memory
clear

% % startup messeage
fprintf('startup.m | (Gaussian) boson sampling simulation\n');
chkmem;
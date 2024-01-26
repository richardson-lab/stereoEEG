rootDir = '/users/tnl/Library/CloudStorage/Box-Box/sEEG';
fileSuffix = '*_induction.mat'; 

% Loop through root directory and return all subject folders containing
% induction file
items = dir(rootDir);
dirs = items([items.isdir] & ~startsWith({items.name}, '.'));
subjects = {};
for i = 1:length(dirs)
    currentDir = fullfile(rootDir, dirs(i).name);
    if ~isempty(dir(fullfile(currentDir, fileSuffix)))
        subjects{end+1} = dirs(i).name;
        
    end
end


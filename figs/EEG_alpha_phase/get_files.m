function files = get_files(path, string)
%% Return all filenames in a directory, that contain a given string.

%%
dir_files = dir(path);

[nFiles, ~] = size(dir_files);

files = {};
j = 1;

for i = 1:nFiles
    if findstr(string, dir_files(i).name)
        files{j} = dir_files(i).name;
        j = j+1;
    end
end



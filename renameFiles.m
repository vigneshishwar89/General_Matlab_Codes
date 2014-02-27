function renameFiles(queryFolder,schar)
files = dir(fullfile(queryFolder));
[m n]=size(files);
for i = 1:m
    fname = files(i).name
    fname_new = handleSplChars(fname,schar);
    fname1 = ['"' fname '"' ' '];
    str = strcat({'mv '}, {fname1}, fname_new)
    system(str{1})
end
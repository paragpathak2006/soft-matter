function A = readMatrix(fileName,n)
    fileID = fopen(fileName);
    A = fscanf(fileID,'%e',[n inf]);
    A = A';
end
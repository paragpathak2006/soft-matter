function writematrix(Matrix,fileName)
fid = fopen(fileName, 'w');
[n, m] = size(Matrix);

for i=1:n
    for j=1:m
        fprintf(fid, '%f\t',Matrix(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
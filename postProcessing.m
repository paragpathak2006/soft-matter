function postProcessing()
clc;
% R0 = [1 3;3 4;5 6;7 8];
% writematrix(R0,'myFile.txt');
Matrix = readMatrix('polynomial_Solution_test.txt',14);

getDbVsLdForGammaVariation(Matrix,0.3) ;
getgammaVsLdForDbVariation(Matrix,0.3);

% getDbVsLdForGammaVariation(Matrix,0.7) ;
% getgammaVsLdForDbVariation(Matrix,0.7);

end
function [R0,R2] = getDbVsLdForGammaVariation(Matrix,ca)
disp('Db Vs Ld (r0 = constant)')
['            ca          ' 'r0           ' 'r1           ' 'r2           ' '*Db*           ' '*ld*']
R0 = 0;
for r1=0.8:0.1:1.2
    r0 = 1 - r1;    r2 = 0;
    op = mask(Matrix,  ca, r0,  r1,  r2, inf);
    
    if size(op,1)==1, continue;end
    
    if R0 == 0,    R0 = op;
    else         R0 = [R0;op];
    end
    
end
R2 = 0;

for r1=0.8:0.1:1.2
    r0 = 0;    r2 = 1 - r1;
    op = mask(Matrix,  ca, r0,  r1,  r2, inf);
    if size(op,1)==1, continue;end
    
    if R2 == 0,    R2 = op;
    else         R2 = [R2;op];
    end
    
end

global NAME;
NAME = strcat('Graph-Db-Vs-Ld-For-Gamma-Variation - R0 - ca - ',num2str(ca*100),'.txt');   writematrix(R0,NAME);
NAME = strcat('Graph-Db-Vs-Ld-For-Gamma-Variation - R2 - ca - ',num2str(ca*100),'.txt');   writematrix(R2,NAME);
if size(R0,1)>1
    R0
end
disp('Db Vs Ld (r2 = constant)')
['            ca          ' 'r0           ' 'r1           ' 'r2           ' '*Db*           ' '*ld*']
if size(R2,1)>1
    R2
end

end

function [R0,R2] = getgammaVsLdForDbVariation(Matrix,ca)
disp('r0 Vs Ld (Db = constant)')
['            ca          ' '*r0*           ' 'r1           ' 'r2           ' 'Db           ' '*ld*']
R0 = 0;
for Db=[1 10]
    op = mask(Matrix,  ca, inf,  inf,  0, Db);
    if size(op,1)==1, continue;end
    
    if R0 == 0,    R0 = op;
    else         R0 = [R0;op];
    end
end
R2 = 0;

for Db=[1 10]
    op = mask(Matrix,  ca, 0,  inf,  inf, Db);
    if size(op,1)==1, continue;end
    
    if R2 == 0,    R2 = op;
    else         R2 = [R2;op];
    end
end

global NAME;
NAME = strcat('Graph-gamma-Vs-Ld-For-Db-Variation - R0 - ca - ',num2str(ca*100),'.txt');   writematrix(R0,NAME);
NAME = strcat('Graph-gamma-Vs-Ld-For-Db-Variation - R2 - ca - ',num2str(ca*100),'.txt');   writematrix(R2,NAME);

if size(R0,1)>1
    R0
end

disp('r2 Vs Ld (Db = constant)')
['            ca          ' 'r0           ' 'r1           ' '*r2*           ' 'Db           ' '*ld*']
if size(R2,1)>1
    R2
end

end

function op = tol(a,b)
if abs(a-b)<1e-10
    op = 1;
else
    op = 0;
end
end
function finalMatrix = mask(Matrix,ca,r0,r1,r2,Db)
%finalMatrix(i,:) = [   ca, r0, r1, r2,  Db ,   ldcrMicro360 , ldcrMicro180, Db, k1avg360 ,k1avg180 , k1minMicro180 , k1maxMicro180 , k1minMicro360 , k1maxMicro360 ]';
%finalMatrix(i,:) = [   1,  2,  3,  4,   5 ,        6 ,            7,        8,      9 ,     10 ,          11 ,            12 ,            13 ,            14 ]';
finalMatrix = [0 0 0 0 0 0];

j = 2;
for i = 1:size(Matrix,1)
    C = tol(Matrix(i,1),ca)||ca==inf;    R0 = tol(Matrix(i,2),r0)||r0==inf;
    R1 = tol(Matrix(i,3),r1)||r1==inf;    R2 = tol(Matrix(i,4),r2)||r2==inf;    D = tol(Matrix(i,5),Db)||Db==inf;

    if ca == inf, x = 1;end;    if r1 == inf, x = 3;end
    if r0 == inf, x = 2;end;    if r2 == inf, x = 4;end;    if Db == inf, x = 5;end
        
    if C&&R0&&R1&&R2&&D
        finalMatrix(j,:) = Matrix(i,1:6);
        j = j + 1;
    end
end

end
function polynomial_Solution_test()

clear; clc; close all; format compact; format shortG; warning('off','all');
global plotFolder plotSuperFolder fileName;    plotSuperFolder = 'C:/Users/ppathak3/Documents/MATLAB/';    plotFolder = plotSuperFolder;
fileName = 'polynomial_Solution_test.txt';
O = -2*ones([1 14]);    finalMatrix = O;
for ca = [0.3 0.5 0.7]
    [r0,r1,r2] = getGammas('Normal');
    opMat = runTest('sun',ca, 2.5 , 0.85,r0,r1,r2);
    finalMatrix = [finalMatrix; opMat; O];
    writematrix(finalMatrix,fileName);
end
fileName = 'extreamities.txt';
[r0,r1,r2] = getGammas('High');    opMat = runTest('L-sun',0.6, 2.5 , 0.85,r0,r1,r2);       finalMatrix = [finalMatrix; opMat; O];         writematrix(finalMatrix,fileName);
[r0,r1,r2] = getGammas('Low');     opMat = runTest('U-sun',0.6, 2.5 , 0.85,r0,r1,r2);       finalMatrix = [finalMatrix; opMat; O];         writematrix(finalMatrix,fileName);

end

%%
function [r0,r1,r2] = getGammas(modType)
     switch modType
        case 'Normal',  r0 = [0.2:-0.1:-0.2 0 0 0 0 0];    r1 = [0.8:0.1:1.2 0.8:0.1:1.2];    r2 = [ 0 0 0 0 0  0.2:-0.1:-0.2];
        case 'High',    r0 = [0.7:-0.1:0.3 0 0 0 0 0];    r1 = [0.3:0.1:0.7 0.3:0.1:0.7];    r2 = [ 0 0 0 0 0  0.7:-0.1:0.3];
        case 'Low',     r0 = [-0.3:-0.1:-0.7 0 0 0 0 0];    r1 = [1.3:0.1:1.7 1.3:0.1:1.7];    r2 = [ 0 0 0 0 0  -0.3:-0.1:-0.7];
        otherwise,      r0 = 0;r1 = 0;r2 = 0;% Ideal DE ( Zhao et al., 2007 );
    end
end
%%
function finalMatrix = runTest(prefix,ca, mub , msmu0b,r0,r1,r2)
    global plotFolder pFile plotSuperFolder NAME;
    % Input parameters
    
    model = getModel('NeoHookeanWithGamma'); ref = 'a'; n = size(r1,2);
    Gmb = 10;   setProp(ca,mub,msmu0b,ref,model,Gmb); % ca = 0.6;       mub = 2.5;    msmu0b = 0.85;
    NAME = strcat('ca - ',num2str(ca*100), ' mub - ',num2str(mub*10),' msmu0b - ',num2str(msmu0b*100),prefix);     mkdir(NAME);
%    pFile = strcat(plotSuperFolder,NAME,'.xlsx');   
    pFile = strcat(NAME,'.txt');   
    O = -1*ones([1 14]);     finalMatrix = O;
    for i = 1:10
        mkdir(strcat(plotSuperFolder,NAME),num2str(i));
        plotFolder = strcat(plotSuperFolder,NAME,'/',num2str(i),'/');
        %----------------------
        opMat = runcase(r0(i),r1(i),r2(i)); 
        %----------------------
        finalMatrix = [finalMatrix; opMat; O]; 
        writematrix(finalMatrix,pFile);  
        % XLW(opMat,'1');   XL(opMat(:,7),'C');XL(opMat(:,11),'N');
    end
end

%%
function XL(matrix,range)
    global pFile;
    ran = strcat(range,'5',':',range,'11');
    writematrix(matrix,pFile,'Sheet','11','Range',ran);
end
%
function XLW(matrix,sheet)
    global pFile;
    writematrix(matrix,pFile,'Sheet',sheet);
end

%%
function finalMatrix = runcase(R0,R1,R2)
    global r0 r1 r2 ldIni;
    r0 = R0;    r1 = R1;    r2 = R2;
    %  Instability Analysis
    %   DbPlot = [0.1 0.5 1:1:10];
    DbPlot = [0.1:0.1:10];
%    DbPlot = [1 10];
    ldIni = 0;              ldcrMacro = 0.89;
    ldcrMicro180 = 0.89;     k1minMicro180 = 15;    k1maxMicro180 = 25;
    ldcrMicro360 = 0.891;    k1minMicro360 = 15;    k1maxMicro360 = 16;

    finalMatrix = runInstabilityAnalysis(DbPlot,ldcrMacro, ldcrMicro180, k1minMicro180, k1maxMicro180, ldcrMicro360, k1minMicro360, k1maxMicro360);
end

%%
function finalMatrix = runInstabilityAnalysis(DbPlot,ldcrMacro, ldcrMicro180, k1minMicro180, k1maxMicro180, ldcrMicro360, k1minMicro360, k1maxMicro360)
    global zoomLevel ldIni  r0 r1 r2;
    ldIni = 0;
    n = size(DbPlot,2);     finalMatrix = zeros([n 14]);
    for i = 1:n
        Db = DbPlot(i);
        %         for zoomLevel = 0:0
        %     %          [k1Plot,ldPlot] = DefaultRangeGenerator(Db);
        %                 [k1Plot,ldPlot] = MacroRangeGenerator(ldcrMacro,i);
        %                 [ldcrMacro , k1minMacro , k1maxMacro] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot,360);
        %                 fprintf('Macro(%i) | Db = %1.1f | ldcr = %f | k1(%1.3f , %1.3f)\n',zoomLevel,Db , ldcrMacro , k1minMacro , k1maxMacro);
        %         end
        ldcrMacro = 0;

        for zoomLevel = 0:0
            [k1Plot,ldPlot] = MicroRangeGenerator(ldcrMicro180,k1minMicro180 , k1maxMicro180,i);
            [ldcrMicro180 , k1minMicro180 , k1maxMicro180] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot,180);
            fprintf('Mi180(%i) | Db = %1.1f | ldcr = %f | k1(%1.3f , %1.3f)\n',zoomLevel,Db , ldcrMicro180 , k1minMicro180 , k1maxMicro180);
        end

%        fprintf('\n');
        for zoomLevel = 0:0
            [k1Plot,ldPlot] = MicroRangeGenerator(ldcrMicro360 , k1minMicro360 , k1maxMicro360,i);
            [ldcrMicro360 , k1minMicro360 , k1maxMicro360] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot,360);
            fprintf('Mi360(%i) | Db = %1.1f | ldcr = %f | k1(%1.3f , %1.3f)\n',zoomLevel,Db , ldcrMicro360 , k1minMicro360 , k1maxMicro360);
        end

        %    zoomLevel = 0;
        %     [ldcr , k1min , k1max] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot,360);
        %     fprintf('Zoom(-), (Db,ldcr) = (%1.1f,%f)  , k1(%1.3f,%1.3f) , type(%i) , ca(%1.2f)\n',Db , ldcr , k1min , k1max, type, ca);

        close all;
        %   ldcrMicro360 = 0; k1minMicro360 =0; k1maxMicro360 = 0;
        %   zoomLevel = 1;ldPlot = ldcr - 0.05:0.001:ldcr + 0.05;    k1Plot = k1min - (k1max - k1min)/10:(k1max - k1min)/20:k1max + (k1max - k1min)/10;    [ldcr , k1min , k1max] = solveMicroscopicInstabilities(ca,ldPlot,k1Plot);    fprintf('ca = %1.1f , ldcr = %1.2f , k1min = %1.1f, k1max = %1.1f \n',ca , ldcr , k1min , k1max);
        [Gma , Gmb , msmu0a , msmu0b ,  xia , xib , ca , cb , mua , mub , mu0 , model,ref]  = getProp();
        k1avg360 =  (k1minMicro360 + k1maxMicro360)/2;
        k1avg180 =  (k1minMicro180 + k1maxMicro180)/2;
        
        finalMatrix(i,:) = [ca, r0,r1,r2,  Db , ldcrMicro360 , ldcrMicro180, Db, k1avg360 ,k1avg180 , k1minMicro180 , k1maxMicro180 , k1minMicro360 , k1maxMicro360 ]';
        fprintf('ca = %1.1f | r0 = %1.1f | r1 = %1.1f | r2 = %1.1f | Db = %f\n',ca,r0,r1,r2,Db);
        fprintf('------------------------------------------------------\n');

    end
end
% model = 3;  setProp(ca,mub,msmu0b,ref,model,Gmb); runInstabilityAnalysis(DbPlot,ldcrMacro, ldcrMicro180, k1minMicro180, k1maxMicro180, ldcrMicro360, k1minMicro360, k1maxMicro360);
% writematrix(finalMatrix,'TensorComponent-A2222.xls');

% model = 4;  setProp(ca,mub,msmu0b,ref,model,Gmb); runInstabilityAnalysis(DbPlot,ldcrMacro, ldcrMicro180, k1minMicro180, k1maxMicro180, ldcrMicro360, k1minMicro360, k1maxMicro360);
% writematrix(finalMatrix,'TensorComponent-A1212.xls');

%
% k1Plot =  finalMatrix(:,5);
% ldPlot =  finalMatrix(:,2);
%
% plotWithFigureTitleSave('ldcr plot',DbPlot,'D_b',ldPlot,'\lambda_{cr}','\lambda_{cr} vs c_m','ldcr-plot.png');
% plotWithFigureTitleSave('k1cr plot',DbPlot,'D_b',ldPlot,'k_1', 'k_1 vs c_m','k1cr-plot.png');
%%
function runSampleTest()
    global zoomLevel;
    zoomLevel = 1;
    ldPlot = 0.1 : 0.01 : 1;
    k1Plot = [0.01:.01:.09 0.1:0.1:5];
    Db = 0.1;
    % [ldcr , k1min , k1max] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot);
    % [ldcrMacro , k1minMacro , k1maxMacro] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot,360);
    solveMicroscopicInstabilities(Db,ldPlot,k1Plot,360);
    solveMicroscopicInstabilities(Db,ldPlot,k1Plot,180);
    % fprintf('ca = %1.1f , ldcr = %1.2f , k1min = %1.1f, k1max = %1.1f \n',ca , ldcrMacro , k1minMacro , k1maxMacro);
end

%%
function [ldcr , k1min , k1max] = solveMicroscopicInstabilities(Db,ldPlot,k1Plot,type)

    global figType;
    ldcr = 0;       k1min = 0;      k1max = 0;
    Notfound = 0;   found = 1;      criticalPoint = Notfound;

    ldPlot = removeUnityStrechRatios(ldPlot);
    [k1Surf,ldSurf,eigPoly,n,m] = createMeshSurface(k1Plot,ldPlot);

    if type == 360
        sign = 1;
        figType = 1;
    else
        sign = -1;
        figType = 2;
    end

    
    for i = n : -1 : 1
        %Inner loop begins
        %-------------------------
        for j = 1 : m
            %-------------------------
            [eig,~] = calculateEigenValues(ldPlot(i),Db,k1Plot(j),type);
            eig = log10(min(max(real(sign*eig),1e-10),1e10)) + 10;
            %        eigPoly(i,j) = LogScale(real(eig));
            eigPoly(i,j) = eig;
            %-------------------------

            if (eig > 0) && criticalPoint == Notfound
                if ldcr == 0
                    ldcr = ldPlot(i);
                    k1min = k1Plot(j);
                else
                    k1max = k1Plot(j);
                end
            end
            %-------------------------

        end
        %Inner loop ends
        %-------------------------

        dk1 = k1Plot(2)-k1Plot(1);
        if ldcr ~= 0
            if(k1max-k1min>=dk1)
                criticalPoint = found;
                %                     return;
            else
                criticalPoint = Notfound;
                ldcr = 0;       k1min = 0;      k1max = 0;
            end
        end

    end

    %outter loop ends
    %-------------------------

    if ldcr == 0    % No cr point
        ldcr = (ldPlot(1) + ldPlot(end))/2;       k1min = k1Plot(1);      k1max = k1Plot(end);
    end

    createCountourPlot(Db,  k1Surf,ldSurf,eigPoly,type,ldcr , k1min , k1max);

end
%%

function createCountourPlot(Db,k1Surf,ldSurf,eigPoly ,theta, ldcr , k1min , k1max)
    global zoomLevel figType plotFolder;
    global r0 r1 r2;

    if(figType == 1), type = 'Macro';
    elseif(figType == 2),type = 'Micro180';
    elseif(figType == 3),type = 'Micro360';
    else type = 'NotDefined';
    end

    [Gma , Gmb , msmu0a , msmu0b ,  xia , xib , ca , cb , mua , mub , mu0 , model,ref] = getProp();
    FigureTitle = strcat(type,', Db-',num2str(Db*sqrt(Gma*mu0)),' , ca-',num2str(ca),' , mub-',num2str(mub),' , msmu0b-',num2str(msmu0b),' , zoom-',num2str(zoomLevel));     figure('Name',FigureTitle);
    plotTitle1 = strcat('D_b - ',num2str(Db*sqrt(Gma*mu0)),'  ||  \lambda_{cr} = ',num2str(ldcr), ' , k1_{min} = ',num2str(k1min), ' , k1_{max} = ',num2str(k1max),'   ||  zoom - ',num2str(zoomLevel), ' , \theta_0 - ',num2str(theta),' (',type,')');
    plotTitle2 = strcat('\gamma_{0} = ',num2str(r0), ' , \gamma_{1} = ',num2str(r1), ' , \gamma_{2} = ',num2str(r2),'  ||  c_m - ',num2str(ca),' , \mu_b - ',num2str(mub),' , m_s\mu_{0b} - ',num2str(msmu0b));
    surf(k1Surf,ldSurf,eigPoly,'FaceAlpha',0.5); grid on;grid minor;
    shading interp;     colormap('prism');  colorbar;           box on;
    ax = gca;           ax.ZGrid = 'on';    ax.XGrid = 'on';    ax.YGrid = 'on';

    % xticks(0:1/10^zoomLevel:100/10^zoomLevel);
    xlabel('k_1');

    % yticks(0:0.05/10^zoomLevel:1);
    ylabel('\lambda');

    % zlim([-100 100]);
    set(get(gca,'ylabel'),'rotation',0);        view(0,90);     title({plotTitle1;plotTitle2});
    saveas(gcf,strcat(plotFolder,FigureTitle,'.png'));

end
%%
function [k1Plot,ldPlot] = DefaultRangeGenerator(Db)
    k1Plot = 0.1 :   0.1     :   10 ;
    ldPlot =  0.1  :   0.01     :  1.2;
end

%%
function [k1Plot,ldPlot] = MacroRangeGenerator(ldcr,isFirstEstimate)
    global zoomLevel;
    Dk1 = 0.2/10^zoomLevel;
    Dld = 0.001/10^zoomLevel;

    %    k1Plot = Dk1 :   Dk1     :   10*Dk1 ;

    isFirstEstimate = 1;

    k1Plot = 0.1 :   0.1     :   10 ;
    if isFirstEstimate == 1
        ldPlot = 0.5  :   0.01     :   2.0;
        %         k1Plot = 0.1 :   0.1     :   20 ;
    else
        ldPlot = ldcr - 10*Dld  :   Dld     :   ldcr + 10*Dld;
    end
end
%%
function [k1Plot,ldPlot] = MicroRangeGenerator(ldcr,k1min,k1max,isFirstEstimate)

    global zoomLevel;
    Dld = 0.001/10^zoomLevel;
%    Dld = 0.05/10^zoomLevel;    k1Plot = 0.1 :   1     :   10 ;


%    isFirstEstimate = 1;
    k1Plot = 0.1 :   0.1     :   10 ;
%    ldPlot = 0.5  :   0.01     :   2.0;
%    return;

    if zoomLevel == 0
        if isFirstEstimate == 1
            ldPlot = 0.5  :   Dld     :   1.5;
        else
 %           k1Plot = k1min - 1  :   0.1     :   k1max + 1;
            ldPlot = [  ldcr - 5*Dld   :   Dld     :   ldcr + 50*Dld];%...
            %                          ldcr + 21*Dld   :   0.01    :   ldcr + 21*Dld + 0.1];
        end
    else
        Dk1 = (k1max - k1min)/50;
        k1Plot = k1min - Dk1*5  :   Dk1     :   k1max + Dk1*5;
        ldPlot = ldcr - 20*Dld  :   Dld     :   ldcr + 20*Dld;
    end

end
%% Run analysis
% model = getModel('NeoHookeanWithGamma');
% Gmb = 19 * 0.0628318531;     mub = 2.5;    msmu0b = 0.85;

% ca = 0.96;       setProp(ca,mub,msmu0b,ref,model,Gmb); finalMatrixCase = runcase('Default');   finalMatrix = [z;finalMatrixCase];writematrix(finalMatrixCase,'Dielectric-Default.xls');
% ca = 0.97;       setProp(ca,mub,msmu0b,ref,model,Gmb); finalMatrixCase = runcase('Default');   finalMatrix = [z;finalMatrixCase];writematrix(finalMatrixCase,'Dielectric-Default.xls');
% ca = 0.98;       setProp(ca,mub,msmu0b,ref,model,Gmb); finalMatrixCase = runcase('Default');   finalMatrix = [z;finalMatrixCase];writematrix(finalMatrixCase,'Dielectric-Default.xls');
% ca = 0.99;       setProp(ca,mub,msmu0b,ref,model,Gmb); finalMatrixCase = runcase('Default');   finalMatrix = [z;finalMatrixCase];writematrix(finalMatrixCase,'Dielectric-Default.xls');

% if(Db<2)
%         ldPlot = 0.5:0.005:1.1;
%         k1Plot = [0.01:.01:.09 0.1:0.1:20];
%     elseif (Db<4)
%         ldPlot = 0.8:0.005:1.5;
%         k1Plot = [0.01:.01:.09 0.1:0.1:20];
%     elseif (Db<7)
%         ldPlot = 1.1:0.005:2;
%         k1Plot = [0.01:.01:.09 0.1:0.1:20];
%     else
%         ldPlot = 1.5:0.005:2.5;
%         k1Plot = [0.01:.02:.09 0.1:0.2:20];
%     end
% end
function plotWithFigureTitleSave(FigTitle,caPlot,xLab,ldPlot,yLab,PlotTitle,SaveTitle)

    figure('Name',FigTitle);     plot(caPlot,ldPlot);       title(PlotTitle);
    xlabel(xLab);               ylabel(yLab);
    saveas(gcf,strcat('./TensorComponent-A1212/',SaveTitle));
    %close all;

end

%%
%% CaPlots
% finalMatrix = zeros([18 5]);
% caPlot = [0.1:0.1:0.9 .91:.01:.99];
% n = size(caPlot,2);
%
% for i = 1:n
%
%     ca = caPlot(i); setProp(ca,mub,msmu0b,ref,model,Gmb);
%     zoomLevel = 0;
%
%     if(ca<.5)
%         ldPlot = 0.5:0.005:1.1;
%         k1Plot = [0.01:.01:.09 0.1:0.1:10];
%     elseif (ca<.9)
%         ldPlot = 0.8:0.005:1.1;
%         k1Plot = [0.01:.01:.09 0.1:0.1:10];
%     elseif (ca<.98)
%         ldPlot = .9:0.001:0.98;
%         k1Plot = [0.01:.01:.09 0.1:0.1:20];
%     else
%         ldPlot = 0.9:0.001:0.98;
%         k1Plot = [0.01:.01:.09 0.1:0.1:25];
%     end
%
%     [ldcr , k1min , k1max] = solveMicroscopicInstabilities(ca,ldPlot,k1Plot);    fprintf('ca = %1.1f , ldcr = %1.2f , k1min = %1.1f, k1max = %1.1f \n',ca , ldcr , k1min , k1max);
%
%     %     zoomLevel = 1;ldPlot = ldcr - 0.05:0.001:ldcr + 0.05;    k1Plot = k1min - (k1max - k1min)/10:(k1max - k1min)/20:k1max + (k1max - k1min)/10;    [ldcr , k1min , k1max] = solveMicroscopicInstabilities(ca,ldPlot,k1Plot);    fprintf('ca = %1.1f , ldcr = %1.2f , k1min = %1.1f, k1max = %1.1f \n',ca , ldcr , k1min , k1max);
%     if(k1min<0.5)
%          finalMatrix(i,:) = [ca , ldcr , k1min , k1max , 0.01]';
%     else
%          finalMatrix(i,:) = [ca , ldcr , k1min , k1max , (k1min + k1max)/2]';
%     end
%
% end
function [k1Surf,ldSurf,eigPoly360,n,m] = createMeshSurface(k1Plot,ldPlot)
    n = size(ldPlot,2);
    m = size(k1Plot,2);
    [k1Surf,ldSurf] = meshgrid(k1Plot,ldPlot);
    eigPoly360 = k1Surf;
    eigPoly180  = k1Surf;
end
%%
function ldPlot1 = removeUnityStrechRatios(ldPlot)
    n = size(ldPlot,2);ldPlot1 = ldPlot;
    %mask = ldPlot(:)<1.001 && ldPlot(:)>.999;
    ldPlot1(abs(ldPlot(:) - 1.0) < 0.01 )  = [];

    % for i = n : -1 : 1
    %     if ldPlot(i) < 1.001
    %         if i>1
    %             ldPlot1(i) = (ldPlot(i) + ldPlot(i-1))/2;
    %         else
    %             ldPlot1(i) = (ldPlot(i) + ldPlot(i+1))/2;
    %         end
    %     end
    % end
end
function y = LogScale(x)
    if (x > 0)
        y = log10(1+x);
    else
        y = 0;
        %         y = -log10(1-x);
    end
end
function i = getModel(modType)
    switch modType
        case 'NeoHookeanWithGamma',     i = 0;
        case 'NeoHookean',     i = 1;
        case 'Langivin',     i = 2;
        case 'ModifiedNeoHookean3',     i = 3;
        case 'ModifiedNeoHookean4',     i = 4;
        case 'Tensors',     i = 5;
        case 'Dielectric',     i = 6;
        otherwise,          i = -1;         % Ideal DE ( Zhao et al., 2007 );
    end

end
%%
%sol = log10(1 + abs(sol));
%sol = min(sol,6);
%            sol = min(real(sol),1);
%            sol = log10(1 + abs(sol));
%            sol = min(sol,6);

%     fprintf('\nld = %f\t',ldPlot(i));

% eigPoly(n,m)

%     for j=1:m
%         fprintf('k1 = %1.1f\t',k1Plot(j));
%     end

% fprintf('\n--------------------------------');
% fprintf('--------------------------------');
% fprintf('--------------------------------');
% fprintf('--------------------------------');
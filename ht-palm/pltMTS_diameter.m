%---------------------------------------
% Change end of variable name for comparison!!!
%---------------------------------------

function [fitVariable, variableForCompm, variableForCompW2FtsI, variableAfterBinning, variableForPlotFtsI, timeFromLeng, timeFromLengAllCell, cellData] = pltMTS_diameter(srInMesh,varargin)

    manClass=[];
    doBadCell=false;
    doManualClassify = false;
    doTMax=false;
    doFilterList=false;
    tMax=Inf;
    nArg=numel(varargin);
    ii = 1;
    while ii <= nArg
        if strcmp(varargin{ii},'ManualClassification')
            doManualClassify=true;
            manClass= varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'TMax')
            doTMax = true;
            tMax= varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'FilterList')
            doFilterList= true;
            filterList= varargin{ii+1};
            ii=ii+2;
        else
            ii=ii+1;
        end
    end

    nCell = numel(srInMesh);
    kk = 1;
    for ii =1:nCell
        contMesh = srInMesh{ii}.mesh;
        tExpt(kk) = srInMesh{ii}.tExpt_min;
        maxD = srInMesh{ii}.diameter.maxD;
        minD = srInMesh{ii}.diameter.minD;
  
        % Anna: get bacteria length of all the cells
        bactLengthVect(ii) = srInMesh{ii}.bactLength;
  
        % Grab the central position of all the cells
        c0(ii, :) = srInMesh{ii}.c0;
 
        plotMinD(kk) = minD(1,2);
        plotMaxD(kk) = min(maxD(:,2));
        W(kk) = plotMinD(kk)/plotMaxD(kk);
        fov(kk) =srInMesh{ii}.fovNo;
        kk=kk+1;
    end

    % bring the mean max diameter of the cells
    maxDiamMean = mean(plotMaxD);

    cellOk=ones(nCell,1);

    if tMax 
        tOk = tExpt<tMax;
        cellOk = cellOk & tOk(:);
    end
    if doFilterList
        cellOk = cellOk & filterList;
    end
    if doManualClassify
        %first exclude the bad cell
        isBad=logical(zeros(size(cellOk)));
        isBad(manClass.badCell) = 1;
        cellOk = cellOk & ~isBad;
        %annotate the septated cells as 'bad cells' 
        %and plot separately
        isDiv = logical(zeros(size(cellOk)));
        isDiv(manClass.postDiv)=1;
        %HACK - list the cells including diameter
        cellData = [(1:nCell)',tExpt(:),plotMinD(:),W(:)];
        cellData(isDiv,3:4)=0;
        cellData(~cellOk,:)=[];

        %END HACK
        % Save idx of cell divided
        cellDivIdx = isDiv;
        variableForCompm.cellDivIdxm = cellDivIdx;
  
   
        %separate classification keeping the manual cells
        cellOk = cellOk & ~isDiv;

        %annotate the manual cells
        tDiv = tExpt(isDiv);
        tDiv = tDiv(:);
        dDiv = zeros(size(tDiv));
        lDiv = bactLengthVect(isDiv); 
        dDiv2 = [plotMinD(isDiv)', plotMaxD(isDiv)', W(isDiv)'];
        %add the non-segmented post-div cells
        pdExt = manClass.postDivExtra;
        nFov = size(pdExt,1);
        for ii = 1:nFov
            nExtra= pdExt(ii,3);
            if nExtra>=0
                tExt = pdExt(ii,2)*ones(nExtra,1);
                dExt = zeros(nExtra,1);
                dExt2 = zeros(nExtra,3);
                tDiv = [tDiv;tExt];
                dDiv = [dDiv;dExt];
                dDiv2 = [dDiv2;dExt2];
            end
        end

    end % end if manual classification
%---------------------------------------------------------------------%

    t= [tExpt(cellOk)';tDiv(:)];
    w= [W(cellOk)';dDiv(:)];
    wPriv = W(cellOk)';
    d = [plotMinD(cellOk)';dDiv(:)];
    a = pi*d.^2/4;
    l = [bactLengthVect(cellOk)'; lDiv(:)];
    lPriv = bactLengthVect(cellOk)';
    tPriv = tExpt(cellOk);
    dPriv = plotMinD(cellOk);
    c0Priv = c0(cellOk, :);

    lengthTemp = lPriv';
    diamTemp = plotMinD(cellOk);

    variableForPlotFtsI = cell(1, 7);
    variableForPlotFtsI.t = t;
    variableForPlotFtsI.w = w;
    variableForPlotFtsI.l = l;
    variableForPlotFtsI.d = d;
    variableForPlotFtsI.tPriv = tPriv';
    variableForPlotFtsI.wPriv = wPriv;
    variableForPlotFtsI.lPriv = lPriv;
    variableForPlotFtsI.dPriv = dPriv';
    
%---------------------------------------------------------------------%
% Measure Lc adn Lg: length at Tc when constriction start and Tg when the
% division occurs
%---------------------------------------------------------------------%

    [tc tg wMax fitOutTimeNoBinNoLeng]  = fitSeamusModBound(t, w);
    [tcF tgF wMaxF fitOutTimeNoBinNoLengF]  = fitFeingoldModBound(t, w);

    % Find cell between Tc and Tg and with a waist diameter smaller than 0.7
    idxOfCellBetweenTcTg = find(tc + 10 <= tPriv' & tPriv' <= tg + 10 & wPriv <= 0.7);


    idxCellWithATimeIntAroundTc = find(tc-1 <= tPriv' & tPriv' <= tc+1);
    idxCellWithATimeIntAroundTg = find(tg-1 <= tPriv' & tPriv' <= tg+1);
    idxCellWithATimeIntAroundTcF = find(tcF-1 <= tPriv' & tPriv' <= tcF+1);
    idxCellWithATimeIntAroundTgF = find(tgF-1 <= tPriv' & tPriv' <= tgF+1);
    idxCellWithATimeIntAroundTo = find(tPriv(1) <= tPriv' & tPriv' <=tPriv(1)+1);

    Lb = mean(l(idxCellWithATimeIntAroundTo));
    Lc = mean(l(idxCellWithATimeIntAroundTc));
    Lg = mean(l(idxCellWithATimeIntAroundTg));
    LcF = mean(l(idxCellWithATimeIntAroundTcF));
    LgF = mean(l(idxCellWithATimeIntAroundTgF));
    
%---------------------------------------------------------------------%
% Take contour of cell divided
%---------------------------------------------------------------------%
    figure,
    cellDivIdx = 1;
    if sum(isDiv)>0
        % save time of cell divided
        tCellDiv = tExpt(isDiv);
        for cellIdx = 1: length(isDiv)   
            clear x;
            clear y;
            clear meshCont;
            
            % find cells divided within a time window minor of
            % ten minutes 
            if isDiv(cellIdx) == 1 && tCellDiv(cellDivIdx) <= min(tCellDiv) + 5
    
                % grab the mesh of the bacteria
                meshCont = srInMesh{cellIdx}.mesh;
                % save it
                variableForCompm.meshAfterDivm{cellDivIdx} = meshCont;

                % grab the palm localizations
                xcol = findSRField(srInMesh{cellIdx}.localizations.parInfo,'semantic','position in sample space in X'); 
                ycol = findSRField(srInMesh{cellIdx}.localizations.parInfo,'semantic','position in sample space in Y'); 
                x = srInMesh{cellIdx}.localizations.data(:,xcol);
                y = srInMesh{cellIdx}.localizations.data(:,ycol);
                palmLoc = [x, y];

                variableForCompm.palmLocm{cellDivIdx} = palmLoc;
                cellDivIdx = cellDivIdx +1;
                % Debug
%                 plot(palmLoc(:, 1), palmLoc(:, 2), 'k*');
%                 hold on,
%                 plot(meshCont(:,1), meshCont(:,2),'y-', 'LineWidth', 2);
%                 plot(meshCont(:,3), meshCont(:,4),'y-', 'LineWidth', 2);
%                 legend('PALM data','meshLeft', 'meshRight' )
%                 axis equal
%                 hold off
%                 drawnow

%                 % Save the contour
%                 contour{:, cellIdx} = [meshCont(:,1)', meshCont(:,2)', meshCont(:,3)', meshCont(:,4)'];
                
            end % end if it is divided
        end
%         variableForCompm.meshAfterDivm = contour;
    end % end contour of cell divided
    %---------------------------------------------------------------------%
        

%---------------------------------------------------------------------%
% Take contour of cell just before division (save only the bigger lobe)
%---------------------------------------------------------------------%
    figure,
    idxCellWithATimeIntBeforeTg = find(tg-20 <= tPriv' & tPriv' <= tg & wPriv <= 0.8);
    c0GoodCell = c0(idxCellWithATimeIntBeforeTg, :);
    
    for cellIdx = 1: length(idxCellWithATimeIntBeforeTg)   
        clear x;
        clear y;
        clear meshCont;
        cellRightIdx = idxCellWithATimeIntBeforeTg(cellIdx);

        meshCont = srInMesh{cellRightIdx}.mesh;
        variableForCompm.meshContm{cellIdx} = meshCont;

        xcol = findSRField(srInMesh{cellRightIdx}.localizations.parInfo,'semantic','position in sample space in X'); 
        ycol = findSRField(srInMesh{cellRightIdx}.localizations.parInfo,'semantic','position in sample space in Y'); 
        x = srInMesh{cellRightIdx}.localizations.data(:,xcol);
        y = srInMesh{cellRightIdx}.localizations.data(:,ycol);
        palmLoc = [x, y];

        variableForCompm.palmLocm{cellIdx} = palmLoc;

        % Debug
        plot(palmLoc(:, 1), palmLoc(:, 2), 'k*');
        hold on,
        plot(meshCont(:,1), meshCont(:,2),'y-', 'LineWidth', 2);
        plot(meshCont(:,3), meshCont(:,4),'y-', 'LineWidth', 2);
        plot(c0GoodCell(cellIdx, 1), c0GoodCell(cellIdx, 2), 'r*')
        legend('PALM data','meshLeft', 'meshRight', 'division center' )
        axis equal
        hold off
        drawnow

        % Find the central line
        cLine = zeros(size(meshCont,1),2);
        cLine(:,1) = mean([meshCont(:,1),meshCont(:,3)],2);
        cLine(:,2) = mean([meshCont(:,2),meshCont(:,4)],2);

        nBreak = 10;

        xC =  cLine(:,1);
        yC =  cLine(:,2);
        cLineSmooth=[xC(:),yC(:)];

        % grab center of division
        c0Cell = c0GoodCell(cellIdx, :);

        % Find the two poles
        xMesh = [meshCont(:, 1); meshCont(:, 2)];
        yMesh = [meshCont(:, 3); meshCont(:, 4)];
        numberOfCoords = size(xMesh, 1);

        maxDistance = zeros(numberOfCoords, 1);
        indexOfMax = zeros(numberOfCoords, 1);
        for k = 1 : numberOfCoords
          distances = sqrt((xMesh-xMesh(k)).^2 + (yMesh-yMesh(k)).^2);
          [maxDistance(k, 1), indexOfMax(k, 1)] = max(distances);
        end
        pole1PosIdx = find(maxDistance == max(maxDistance));
        pole2PosIdx = indexOfMax(pole1PosIdx);

        % Find the point of the central line closest to the center c0Cell
        distancesFromC = sqrt((xC - c0Cell(1, 1)).^2 + (yC - c0Cell(1, 2)).^2);
        [minDistance, indexOfMin] = min(distancesFromC);

        c0Mesh = [ xC(indexOfMin); yC(indexOfMin)];

        % Anna: compute bacteria pole-center length normalize by the tot lenght
        totLength = arclength( xC, yC, 'spline');
        p1CenterLength = arclength( xC(1 : indexOfMin), yC(1 : indexOfMin), 'spline')/totLength;
        p2CenterLength = arclength( xC(indexOfMin : end), yC(indexOfMin : end), 'spline')/totLength;

        if p1CenterLength >= p2CenterLength 
            
            contourPriv{:, cellIdx} = [meshCont(1:indexOfMin, 1)', meshCont(1:indexOfMin, 2)', meshCont(1:indexOfMin, 3)', meshCont(1:indexOfMin, 4)'];
            % Debug
%             plot(palmLoc(:, 1), palmLoc(:, 2), 'k*');
%             hold on,
%             plot(meshCont(1:indexOfMin, 1), meshCont(1:indexOfMin, 2),'y-', 'LineWidth', 2);
%             plot(meshCont(1:indexOfMin, 3), meshCont(1:indexOfMin, 4),'y-', 'LineWidth', 2);
%             legend('PALM data','meshLeft', 'meshRight' )
%             axis equal
%             hold off
%             drawnow
            
        else
          contourPriv{:, cellIdx} = [meshCont(indexOfMin:end, 1)', meshCont(indexOfMin:end, 2)', meshCont(indexOfMin:end, 3)', meshCont(indexOfMin:end, 4)'];
        end    

    end % end cylce over the cell just before division
    
    % Save the contour of the biggest lobe of the cell just before the
    % division
    variableForCompm.meshPrivDivm = contourPriv;

%---------------------------------------------------------------------%
% Take the distance pole-center of division
%---------------------------------------------------------------------%
    c0GoodCell = c0(idxOfCellBetweenTcTg, :);
    cellIdxTemp = 1;
    figure, 
  for cellIdx = 1: length(c0GoodCell) 
      
       cellIdxGood = idxOfCellBetweenTcTg(cellIdx);
       clear x;
       clear y;
        xcol = findSRField(srInMesh{cellIdxGood}.localizations.parInfo,'semantic','position in sample space in X'); 
        ycol = findSRField(srInMesh{cellIdxGood}.localizations.parInfo,'semantic','position in sample space in Y'); 
        x = srInMesh{cellIdxGood}.localizations.data(:,xcol);
        y = srInMesh{cellIdxGood}.localizations.data(:,ycol); 
        
        mesh = srInMesh{cellIdxGood}.mesh;
        cLine = zeros(size(mesh,1),2);
        cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
        cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);
        
        nBreak = 10;
        % smooth the centreLine
        [xC, yC] = spline2d(cLine(:,1),cLine(:,2),nBreak);
        cLineSmooth=[xC(:),yC(:)];
        
        % grab center of division
        c0Cell = c0GoodCell(cellIdx, :);
        
        % Find the two poles
        xMesh = [mesh(:, 1); mesh(:, 2)];
        yMesh = [mesh(:, 3); mesh(:, 4)];
        numberOfCoords = size(xMesh, 1);
        
        maxDistance = zeros(numberOfCoords, 1);
        indexOfMax = zeros(numberOfCoords, 1);
        for k = 1 : numberOfCoords
            distances = sqrt((xMesh-xMesh(k)).^2 + (yMesh-yMesh(k)).^2);
            [maxDistance(k, 1), indexOfMax(k, 1)] = max(distances);
        end
        pole1PosIdx = find(maxDistance == max(maxDistance));
        pole2PosIdx = indexOfMax(pole1PosIdx);

        % Find the point of the central line closest to the center c0Cell
        distancesFromC = sqrt((xC - c0Cell(1, 1)).^2 + (yC - c0Cell(1, 2)).^2);
        [minDistance, indexOfMin] = min(distancesFromC);
        
        c0Mesh = [ xC(indexOfMin); yC(indexOfMin)];
        

        % Anna: compute bacteria pole-center length normalize by the tot lenght
        totLength = arclength( xC, yC, 'spline');
        p1CenterLength = arclength( xC(1 : indexOfMin), yC(1 : indexOfMin), 'spline')/totLength;
        p2CenterLength = arclength( xC(indexOfMin : end), yC(indexOfMin : end), 'spline')/totLength;
        
        % Filter out cell that did not start constriction and have a random
        % center
        if ( p1CenterLength <= 0.6 && p1CenterLength >= 0.4 && p2CenterLength <= 0.6 && p2CenterLength >= 0.4 )
        
            % Save p1C and p2C accordingly with their lenght 
            if p1CenterLength >= p2CenterLength
                p1C0(cellIdxTemp) = p1CenterLength;
                p2C0(cellIdxTemp) = p2CenterLength;
               
            else
                p1C0(cellIdxTemp) = p2CenterLength;
                p2C0(cellIdxTemp) = p1CenterLength;
            end     
                 
            tFinal(cellIdxTemp) = tPriv(cellIdxGood);
            cellIdxTemp = cellIdxTemp + 1;
               
          % Debug
          plot(x, y, 'k*'); hold on
          plot(mesh(:,1), mesh(:,2),'y-', 'LineWidth', 2);hold on
          plot(mesh(:,3), mesh(:,4),'y-', 'LineWidth', 2);
          plot(xC, yC, 'b-')
          plot(c0Cell(1,1), c0Cell(1,2),'r*')
          plot(xC(1), yC(1), 'ro')
          plot(xC(end), yC(end), 'ro')
          legend('PALM data','meshLeft', 'meshRight', 'centralLine', 'centerDiv', 'pole1', 'pole2' )
          axis equal
          hold off
          drawnow 
       end
  end
  
        p1c0Mean = mean(p1C0, 2);
        p2c0Mean = mean(p2C0, 2);
        % Save the average value of p1c0 and p2c0
        variableForCompm.p1c0Mean = p1c0Mean;
        variableForCompm.p2c0Mean = p2c0Mean;
        
        % fit data
        [par1,stats1] = robustfit(tFinal,p1C0);
        std_robust1 = stats1.robust_s;
        
        [par2,stats2] = robustfit(tFinal,p2C0);
        std_robust2 = stats2.robust_s;
        
        nSigma = 1;
        figure,
        plot(tFinal,  p1C0, 'r*'), hold on,
        plot(tFinal,  p2C0, 'b*'), 
        plot(tFinal, par1(2)*tFinal  + par1(1), 'r-')
        plot(tFinal, par2(2)*tFinal  + par2(1), 'b-')
        plot(tFinal, par1(2)*tFinal  + par1(1) + std_robust1*nSigma, 'r--')
        plot(tFinal, par1(2)*tFinal  + par1(1) - std_robust1*nSigma, 'r--')
        plot(tFinal, par2(2)*tFinal  + par2(1) + std_robust2*nSigma, 'b--')
        plot(tFinal, par2(2)*tFinal  + par2(1) - std_robust2*nSigma, 'b--')
        xlabel('Time [min]')
        ylabel('Pole-Center Distance normalize by the tot lenght')
        legend('Pole1-Center Distance', 'Pole2-Center Distance','Mean of the pole 1 - center distance', 'Mean of the pole 2 - center distance', 'Conf. bound 68%','Conf. bound 68%','Conf. bound 68%', 'Conf. bound 68%')
        
        
        % fit normalized data
        [par1,stats1] = robustfit(tFinal/tg, p1C0);
        std_robust1 = stats1.robust_s;
        
        [par2,stats2] = robustfit(tFinal/tg,p2C0);
        std_robust2 = stats2.robust_s;
        figure,
        plot(tFinal/tg,  p1C0, 'r*'), hold on,
        plot(tFinal/tg,  p2C0, 'b*'), 
        plot(tFinal/tg, par1(2)*tFinal/tg  + par1(1), 'r-')
        plot(tFinal/tg, par2(2)*tFinal/tg  + par2(1), 'b-')
        plot(tFinal/tg, par1(2)*tFinal/tg  + par1(1) + std_robust1*nSigma, 'r--')
        plot(tFinal/tg, par1(2)*tFinal/tg  + par1(1) - std_robust1*nSigma, 'r--')
        plot(tFinal/tg, par2(2)*tFinal/tg  + par2(1) + std_robust2*nSigma, 'b--')
        plot(tFinal/tg, par2(2)*tFinal/tg  + par2(1) - std_robust2*nSigma, 'b--')
        xlabel('Normalize Time (t/tg) ')
        ylabel('Pole-Center Distance normalize by the tot lenght')
        legend('Pole1-Center Distance', 'Pole2-Center Distance','Mean of the pole 1 - center distance', 'Mean of the pole 2 - center distance', 'Conf. bound 68%','Conf. bound 68%','Conf. bound 68%', 'Conf. bound 68%') 
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
% Take the curvature of the central line
%---------------------------------------------------------------------%
    
    figure,
    cellIdxGood = 1;
    for cellIdx = 1: length(wPriv) 
     
       % if wPriv(cellIdx) < 0.8 
            % Grab the division center
            c0Cell = c0(cellIdx, :);
            % Grab the contour and the central line
            mesh = srInMesh{cellIdx}.mesh;
            cLine = zeros(size(mesh,1),2);
            cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
            cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);

            nBreak = 5;
            % smooth the centreLine
            [xC, yC] = spline2d(cLine(:,1),cLine(:,2),nBreak);
            %cLineSmooth=[xC(:),yC(:)];

            [k, dfdx, d2fdx2, x, y] = getCurvature(xC, yC);
            
            % Grab the mean curvature
            goodCurvIdx = find(k < mean(k) + std(k) & k > mean(k) - std(k));
            meanK = mean(k(goodCurvIdx));
            meanKAll(cellIdxGood) = meanK;
            timeKAll(cellIdxGood) = tPriv(cellIdx);
            lengthKAll(cellIdxGood) = lPriv(cellIdx);
            
            cellIdxGood = cellIdxGood +1;
            % Find the point of the central line closest to the center c0Cell
            distancesFromC = sqrt((xC - c0Cell(1, 1)).^2 + (yC - c0Cell(1, 2)).^2);
            [minDistance, indexOfMin] = min(distancesFromC);

            xCirc = (0:pi/64:2*pi);

            % Plot the central line with the osculating circle
            subplot(1, 2, 1)
            plot(mesh(:,1), mesh(:,2),'y-', 'LineWidth', 2);hold on
            plot(mesh(:,3), mesh(:,4),'y-', 'LineWidth', 2);
            plot(xC, yC, 'b-')
            plot(xC(indexOfMin), yC(indexOfMin), 'r*')
            plot(c0Cell(1,1), c0Cell(1,2), 'c*')
            plot(x, y', 'g+')
            xlabel('[nm]')
            ylabel('[nm]')
            
            % Plot osculating circle
            R = 1/meanK;
            d2fdx2Mean =  mean(d2fdx2);
            dfdxMean = mean(dfdx);
            if d2fdx2Mean > 0
               if  dfdxMean > 0
                    Xc = xC(indexOfMin) - R*sin( atan(dfdxMean) );
                    Yc = yC(indexOfMin) + R*cos( atan(dfdxMean) );
                    plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
                    axis equal
               else
                    Xc = xC(indexOfMin) - R*sin( atan(dfdxMean) );
                    Yc = yC(indexOfMin) + R*cos( atan(dfdxMean) );
                    plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
                    axis equal
               end

            else
                if dfdxMean > 0
                    Xc = xC(indexOfMin) + R*sin( atan(dfdxMean) );
                    Yc = yC(indexOfMin) - R*cos( atan(dfdxMean) );
                    plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
                    axis equal
                else
                    Xc = xC(indexOfMin) + R*sin( atan(dfdxMean) );
                    Yc = yC(indexOfMin) - R*cos( atan(dfdxMean) );
                    plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
                    axis equal

                end
            end
            legend('meshLeft', 'meshRight', 'centralLine', 'Division Center', 'Division Center','centralLineRotInsideFunc', 'Osculating circle' )
            axis equal
            hold off
            
            subplot(1, 2, 2)
            hist(k, 5)
            xlabel('Curvature k [1\nm]')
            ylabel('Frequency #')
            drawnow 
%             pause
            disp('cell:'), disp(cellIdx)
       % end % if the cell has a small diameter w<0.8
    end % end fro cycle for curvature 
    
    % Plot curvature k as a function of time
    figure,
    plot(timeKAll, meanKAll, 'b*')
    xlabel('Time[min]')
    ylabel('Curvature k [1\nm]')
    title('Curvature as a function of time [min]')
    
    % Plot curvature(osculating circle R) as a function of time
    figure,
    plot(timeKAll/tg, meanKAll.^-1, 'b*')
    xlabel('Time[min]')
    ylabel('Curvature R[nm]')
    title('Curvature as a function of normalized time t/tg [min]')
    
    % Plot curvature(osculating circle R) as a function of time
    figure,
    plot(lengthKAll/tg, meanKAll.^-1, 'b*')
    xlabel('Length [nm]')
    ylabel('Curvature R[nm]')
    title('Curvature as a function of normalized time t/tg [min]')
    
    lengthKAll
    
    variableForCompm.curvAndTimem = [timeKAll', meanKAll'];
    
%Do binning for length vs time
% tBinsAllCell=transpose(round(min(t)-2):14:round(max(t)) -20);
% for i=1:length(tBinsAllCell)-1;
%     i
%     targetAllCell= find(tBinsAllCell(i)<=t & t<=tBinsAllCell(i+1));
%     lBinsAllCell(i)=mean(l(targetAllCell));
%     wBinsAllCellNoTimeFromLeng(i)=mean(w(targetAllCell));
%     %errlBinsAllCell(i)= std(l(targetAllCell));
%    % errwBinsAllCell(i)= std(w(targetAllCell));
%     tBinsAllCellCentered(i) = (tBinsAllCell(i)+ tBinsAllCell(i+1))/2;
%     lAfterBinningAllCell(:,i) = ( tBinsAllCell(i)<t & t<tBinsAllCell(i+1) )*lBinsAllCell(i);
%     numElePerBinAll(i) = length(targetAllCell);
%     wAfterBinningAllCell(:,i) = ( tBinsAllCell(i)<t & t<tBinsAllCell(i+1) )*wBinsAllCellNoTimeFromLeng(i);
% end
% lAfterBinningAllCellFin = sum(lAfterBinningAllCell,2);
% % wAfterBinningAllCellNoTimeFromLeng = sum(wAfterBinningAllCell,2);
% 
% disp('Average number of elements for bin (all cell -length vs time):');
% disp(mean(numElePerBinAll));

% Do binning w vs time
twBins=transpose(round(min(t)-2):10:round( max(t)-20) );
wAfterBinningAllCell = zeros(length(t), length(twBins)-1);
for i=1:length(twBins)-1;
    
    target= find(twBins(i)<=t & t<=twBins(i+1));

    wtBins(i)=mean(w(target));
   % errlBins(i)= std(l(target));
   % errwBins(i)= std(w(target));
    twBinsCentered(i) = (twBins(i)+ twBins(i+1))/2;
    numElePerBin(i) = length(target);
end

% Do binning length vs time
tBins=transpose(round(min(tPriv)-2):10:round( max(tPriv)-20) );
lAfterBinning = zeros(length(tPriv), length(tBins)-1);
for i=1:length(tBins)-1;
    
    target= find(tBins(i)<=tPriv & tPriv<=tBins(i+1));

    lBins(i)=mean(lPriv(target)); 
    wBinsNoTimeFromLeng(i)=mean(wPriv(target));
   % errlBins(i)= std(l(target));
   % errwBins(i)= std(w(target));
    tBinsCentered(i) = (tBins(i)+ tBins(i+1))/2;
    lAfterBinning(:,i) = (tBins(i)<tPriv & tPriv<tBins(i+1))*lBins(i);
    wAfterBinning(:,i) = ( tBins(i)<tPriv & tPriv<tBins(i+1) )*wBinsNoTimeFromLeng(i);
    numElePerBin(i) = length(target);
end
lAfterBinningFin = sum(lAfterBinning,2);
wAfterBinningNoTimeFromLeng = sum(wAfterBinning,2);
%wAfterBinningAllCellNoTimeFromLeng = [wAfterBinningNoTimeFromLeng; dDiv(:)];

disp('Average number of elements for bin (length vs time):');
disp(mean(numElePerBin));

fitTimeLengthMod = polyfit(tExpt(cellOk), lengthTemp, 2);
%fitTimeLengthMod = polyfit(tBinsCentered, lBins, 2);
%lFun =polyval(fitTimeLengthMod, tBinsCentered);
lFun =polyval(fitTimeLengthMod, tBinsCentered);
lFunAll =polyval(fitTimeLengthMod, tExpt(cellOk));

% DEBUG plot with the fitting done with "polyfit"
figure
plot(tBinsCentered, lBins, 'c*')
hold on
plot(tExpt(cellOk), lengthTemp, 'b*')
plot(tBinsCentered,lFun, 'r--')
plot(tExpt(cellOk),lFunAll, 'g--')
xlabel('Time [min]');
ylabel('Lenght [nm]');
legend('Cells', 'Fit after binning')

% Linear Fit using the population of cell
fitTimeLengthMod = polyfit(tExpt(cellOk), lengthTemp, 1);
fitTimeLengthMod2 = polyfit(tBinsCentered, lBins, 2);
%fitTimeLengthMod2 = griddedInterpolant(tBinsCentered, lBins,'spline');


timeLegthMod = @(leng) leng./fitTimeLengthMod(1) - fitTimeLengthMod(2)./fitTimeLengthMod(1);
timeLegthMod2 = @(leng)  (- fitTimeLengthMod2(2) + sqrt( fitTimeLengthMod2(2).^2 - 4.*fitTimeLengthMod2(1)*(fitTimeLengthMod2(3)-leng)) )...
    ./(2*fitTimeLengthMod2(1));


% DEBUG plot done with linear fit
figure,
plot(tExpt(cellOk), lengthTemp, 'b*')
hold on
%plot(invFuncLengVsTime(lBins),lBins, 'r--')
plot(timeLegthMod(lBins),lBins, 'r--')
plot(timeLegthMod2(lBins),lBins, 'g--')
xlabel('Time [min]');
ylabel('Lenght [nm]');
legend('Cells', 'Fit after binning')

% fit using the bin data
fitTimeLengthModAfterBin = polyfit(tBinsCentered, lBins, 1);
%fitTimeLengthMod = griddedInterpolant(tExpt(cellOk), lengthTemp);
timeLegthModAfterBin = @(leng) leng./fitTimeLengthModAfterBin(1) - fitTimeLengthModAfterBin(2)./fitTimeLengthModAfterBin(1);


if fitTimeLengthMod2(1)>0 % time computed from length (fit with ax^2+bx+c)from all population
    timeFromLeng = timeLegthMod2(lPriv);
    timeFromLengAllCell = [timeFromLeng; tDiv(:)] ;
else % time computed from length (from lin fit)
    timeFromLeng = timeLegthMod(lPriv);
    timeFromLengAllCell = [timeFromLeng; tDiv(:)] ;
   %timeFromLengAllCell = timeLegthMod(l);
end

% time computed from length (from lin fit) from bin data
timeFromLengAfterBinning = timeLegthModAfterBin(lAfterBinningFin );
timeFromLengAllCellAfterBinning = [timeFromLengAfterBinning; tDiv(:)];

% Do the binning of w vs time where time computed from length
tBinsAllCellFromLeng=transpose(round(min(timeFromLengAllCell)-2):12:round(max(timeFromLengAllCell)+10));
for i=1:length(tBinsAllCellFromLeng)-1;
    targetAllCellFromLeng= find(tBinsAllCellFromLeng(i)<timeFromLengAllCell & timeFromLengAllCell<tBinsAllCellFromLeng(i+1));
    wBinsAllCell(i)=mean(w(targetAllCellFromLeng));
    errwBinsAllCell(i)= std(w(targetAllCellFromLeng));
    numElePerBinWaist(i) = length(targetAllCellFromLeng);
    tBinsAllCellCenteredFromLeng(i) = (tBinsAllCellFromLeng(i)+ tBinsAllCellFromLeng(i+1))/2;
end

% tBinsAllCellFromLeng=transpose(round(min(timeFromLengAllCellAfterBinning)-2):34:round(max(timeFromLengAllCellAfterBinning)+10));
% for i=1:length(tBinsAllCellFromLeng)-1;
%     targetAllCellFromLeng= find(tBinsAllCellFromLeng(i)<timeFromLengAllCellAfterBinning & timeFromLengAllCellAfterBinning<tBinsAllCellFromLeng(i+1));
%     wBinsAllCell(i)=mean(w(targetAllCellFromLeng));
%     errwBinsAllCell(i)= std(w(targetAllCellFromLeng));
%     numElePerBinWaist(i) = length(targetAllCellFromLeng);
%     tBinsAllCellCenteredFromLeng(i) = (tBinsAllCellFromLeng(i)+ tBinsAllCellFromLeng(i+1))/2;
% end


disp('Average number of elements for bin (waist vs time):');
disp(mean(numElePerBinWaist));

tBinsFromLeng=transpose(round(min(timeFromLengAfterBinning)-2):12:round(max(timeFromLengAfterBinning)+10));
for i=1:length(tBinsFromLeng)-1;
    targetFromLeng= find(tBinsFromLeng(i)<timeFromLengAfterBinning & timeFromLengAfterBinning<tBinsFromLeng(i+1));
    wBins(i)=mean(w(targetFromLeng));
    errwBins(i)= std(w(targetFromLeng));
    tBinsCenteredFromLeng(i) = (tBinsFromLeng(i)+ tBinsFromLeng(i+1))/2;
end
% wAfterBinningFin = sum(wAfterBinning,2);
% wAfterBinningAllCellFin = [wAfterBinningFin; dDiv(:)];


variableAfterBinning = cell(1, 6);
variableAfterBinning.timeFromLengAfterBinning = timeFromLengAfterBinning;
variableAfterBinning.timeFromLengAllCellAfterBinning = timeFromLengAllCellAfterBinning;

%variableAfterBinning.lAfterBinningAllCellFin = lAfterBinningAllCellFin;
%variableAfterBinning.wAfterBinningAllCellFin = wAfterBinningAllCellFin; 
variableAfterBinning.lAfterBinningFin = lAfterBinningFin;
%variableAfterBinning.wAfterBinningFin = wAfterBinningFin; 

% Do not include cell with negative time 
negTimeIdx = find(timeFromLengAllCellAfterBinning<0);
negTimeIdxForBin = find(tBinsFromLeng<0);

timeFromLengAfterBinningNoNeg = timeFromLengAfterBinning;
timeFromLengAfterBinningNoNeg(negTimeIdx) = [];

% tBinsFromLengFinNoNeg = tBinsCenteredFromLeng;
% tBinsFromLengFinNoNeg(negTimeIdxForBin) = [];
% wAfterBinningFinNoNeg = wAfterBinningFin; 
% wAfterBinningFinNoNeg(negTimeIdxForBin) = [];
% wBinsFinNoNeg = wBins;
% wBinsFinNoNeg(negTimeIdxForBin) = [];

% Do not include cell with negative time (all cells)
negTimeIdxAllCell = find(timeFromLengAllCellAfterBinning<0);
negTimeIdxAllCellForBin = find(tBinsAllCellFromLeng<0);

timeFromLengAllCellAfterBinningNoNeg = timeFromLengAllCellAfterBinning;
timeFromLengAllCellAfterBinningNoNeg(negTimeIdxAllCell) = [];

tBinsAllCellFromLengFinNoNeg = tBinsAllCellCenteredFromLeng;
tBinsAllCellFromLengFinNoNeg(negTimeIdxAllCellForBin) = [];
% wAfterBinningAllCellFinNoNeg = wAfterBinningAllCellFin; 
% wAfterBinningAllCellFinNoNeg(negTimeIdxAllCell) = [];
wBinsAllCellFinNoNeg = wBinsAllCell;
wBinsAllCellFinNoNeg(negTimeIdxAllCellForBin) = [];
errwBinsAllCellNoNeg =  errwBinsAllCell;
errwBinsAllCellNoNeg(negTimeIdxAllCellForBin) = [];

 variableAfterBinning.timeFromLengAllCellAfterBinningNoNeg = timeFromLengAllCellAfterBinningNoNeg;
 variableAfterBinning.timeFromLengAfterBinningNoNeg = timeFromLengAfterBinningNoNeg;
 variableAfterBinning.tBinsAllCellFromLengFinNoNeg = tBinsAllCellFromLengFinNoNeg;
% variableAfterBinning.tBinsFromLengFinNoNeg = tBinsFromLengFinNoNeg;
 variableAfterBinning.wBinsAllCellFinNoNeg = wBinsAllCellFinNoNeg; 
% variableAfterBinning.wBinsFinNoNeg = wBinsFinNoNeg; 

negTimeIdx = find(timeFromLeng<0);
timeFromLengNoNeg = timeFromLeng;
timeFromLengNoNeg(negTimeIdx) = [];
timeFromLengAllCellNoNeg = [timeFromLengNoNeg; tDiv(:)];

wPrivNoNeg = wPriv; 
wPrivNoNeg(negTimeIdx) = [];
variableForPlotFtsI.timeFromLeng = timeFromLeng;
variableForPlotFtsI.timeFromLengAllCell = timeFromLengAllCell;

wAllCellNoNeg = [wPrivNoNeg; dDiv(:)];

variableForPlotFtsI.wPrivNoNeg = wPrivNoNeg;
variableForPlotFtsI.timeFromLengNoNeg = timeFromLengNoNeg;

variableForPlotFtsI.wAllCellNoNeg = [wPrivNoNeg; dDiv(:)];
variableForPlotFtsI.timeFromLengAllCellNoNeg = [timeFromLengNoNeg; tDiv(:)];

variableForPlotFtsI.tBinsAllCellFromLengFinNoNeg = tBinsAllCellFromLengFinNoNeg;
variableForPlotFtsI.wBinsAllCellFinNoNeg = wBinsAllCellFinNoNeg; 
 
variableForPlotFtsI.twBins = twBinsCentered;
variableForPlotFtsI.wtBins = wtBins';
%-------------------------------------------------------------------------%
%  define rate (diameter costrinction rate respect to time and lenght)
%-------------------------------------------------------------------------%

% [deriv1DiamVsTime, deriv2DiamVsTime] = dfdxc(tDiv,dDiv);
% [deriv1DiamVsLeng, deriv2DiamVsLeng] = dfdxc(lPriv,dDiv);
%[deriv1DiamVsLeng, deriv2DiamVsLeng] = dfdxc(l,dDiv);

%--------------------------------------------------------------------------%
% Anna: plot septum diameter vs baceria length 
figure,
plot(bactLengthVect(cellOk), plotMinD(cellOk), 'm*')
title('Caulobacter septum diameter vs length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(diff(lengthTemp), diamTemp(1: length(diamTemp)-1), 'b*')
title('Caulobacter septum diameter vs elongation', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Elongation [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(lengthTemp(1:length(lengthTemp)-1), diff(diamTemp), 'b*')
title('Variation of the septum diameter vs length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Diff Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(tPriv(1:length(tPriv)-1), diff(diamTemp), 'b*')
title('Variation of the septum diameter vs length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Diff Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(diff(lengthTemp), diff(diamTemp), 'b*')
title('Variation of the septum diameter respect to the length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Elongation [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Diff Septum diameter [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

figure,
plot(tExpt(cellOk), lengthTemp, 'b*')
hold on
title('Caulobacter length vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('Time [min]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Length [nm]', 'FontSize',16, 'FontName', 'Copperplate Gothic Light')

plot(timeFromLeng, lengthTemp, 'r--')
text(400, 3000, ['$$F(L)= ' num2str(round(fitTimeLengthMod(1))) '*t + ' num2str(round(fitTimeLengthMod(2))) '$$'],'interpreter','latex','FontSize',12 )
text(400, 2500, '$$ F_{model}(l) = m*t + q $$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')
legend('data', 'fit')
hold off

%--------------------------------------------------------------------------%
figure;
subplot(2,1,1);
hold on;
plot(tExpt(cellOk),plotMinD(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv,'ko');
end
ylim([0 600]);
legend('Caulobacters before division', 'Caulobacters post division')
title('Caulobacter septum diameter vs time - with post divisional cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Septum diameter [nm]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
subplot(2,1,2);
hold on;
plot(tExpt(cellOk),W(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv,'ko');
end
ylim([0 1.2]);
title('Caulobacter septum waist ratio r/R vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
ylabel('Waist ratio r/R (septum diameter/max diameter)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );


figure;
subplot(2,1,1);
hold on;
plot(tExpt(cellOk),plotMinD(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv2(:,1),'ro');
end
ylim([0 600]);
legend('Caulobacters before division', 'Caulobacters post division')
title('Caulobacter septum diameter vs time - with post divisional cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Septum diameter [nm]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
subplot(2,1,2);
hold on;
plot(tExpt(cellOk),W(cellOk),'ko');
if doManualClassify
  plot(tDiv,dDiv2(:,3),'ro');
end
ylim([0 1.2]);
legend('Caulobacters before division', 'Caulobacters post division')
title('Caulobacter septum waist ratio r/R vs time - with post divisional cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
ylabel('Waist ratio r/R (septum diameter/max diameter)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );



figure;
plot(t,a,'ko');
title('Caulobacter area vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
xlabel('T [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
ylabel('Area [nm^2]','FontSize',16, 'FontName', 'Copperplate Gothic Light' );


% Anna: adding fitting using seamus functions
%
% fitVariable:
% 1)- fitOutTimeNoBinNoLeng: w vs t
% 2)- fitOutLength: w vs t computed from length
% 3)- fitOutTimeBinWandL: w vs t after length and w binning (all cell)
% 4)- fitOutTimeBinW: w vs t after only w binning (all cell) 
% 5)- fitOutTimeBinWandLNoNegAllCell: w vs t after length and w binning (all
%                                   cell) (no negative time)
% 6)- fitOutTimeBinWandLNoCellAftNoNeg: w vs t after length and w binning (no cells after division)(no negative time)
% 7)- fitOutTimeOnlyLengthBin: w vs t after only length binning (no cells after division)
% 8)- fitOutLengthNoNeg: wPriv vs time computed from length ( no negative time)

fitVariable = cell(1, 6);

% 1) Anna: fit w vs t

% [tc tg1 wMax1 fitOutTime1NoBinNoLeng]  = fitFeingold(t, w);
% [tg wMax fitOutTimeNoBinNoLeng]  = fitFeingoldFixedStart(t, w, tc);

[tc tg wMax fitOutTimeNoBinNoLeng]  = fitSeamusModBound(t, w);
[tcF tgF wMaxF fitOutTimeNoBinNoLengF]  = fitFeingoldModBound(t, w);


fitVariable.fitOutTimeNoBinNoLeng = fitOutTimeNoBinNoLeng;
fitVariable.fitOutTimeNoBinNoLengTc = tc;
fitVariable.fitOutTimeNoBinNoLengTg = tg;
fitVariable.fitOutTimeNoBinNoLengWmax = wMax;

variableForCompm.tm=t;
variableForCompm.wm=w;
variableForCompm.fitWvsTm = fitOutTimeNoBinNoLeng;
variableForCompm.fitWvsTparTCm = tc;
variableForCompm.fitWvsTparTGm = tg;
variableForCompm.fitWvsTparWmaxm = wMax;
variableForCompm.binDataTWm = twBinsCentered;
variableForCompm.binDataWTm = wtBins';

figure
plot(fitOutTimeNoBinNoLeng(:, 1), fitOutTimeNoBinNoLeng(:, 2), 'r')
hold on
plot(t, w, 'bo')
plot(twBinsCentered, wtBins', 'g*')

% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')

text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )}$$'],'interpreter','latex', 'FontSize',12 )
text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )}$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Caulobacter septum waist ratio r/R vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with 2D septum area model', 'Caulobacter during division', 'Binned data')
xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off

%-----------------------------------------------------------------------------------
% Area plot
%------------------------------------------------------------------------------

figure
%plot(wMax^2.*( (fitOutAreaTimeNoBinNoLeng(:, 1)-tc)./(tg-tc) ), wMax^2 - fitOutAreaTimeNoBinNoLeng(:, 2), 'r')
plot( fitOutTimeNoBinNoLeng(:, 1), wMax^2 - fitOutTimeNoBinNoLeng(:, 2).^2, 'r')
hold on
%plot(wMaxF^2.*( (fitOutAreaFeingold(:, 1)-tcF)./(tgF -tcF) ) , wMaxF^2 - fitOutAreaFeingold(:, 2), 'c')
plot( (fitOutTimeNoBinNoLengF(:, 1) ), wMaxF^2 - fitOutTimeNoBinNoLengF(:, 2).^2, 'c')
%plot( wMax^2.*(t-tc)./(tg-tc), wMax^2 - w.^2, 'bo')
plot( t, wMax^2 - w.^2, 'bo')
title('Caulobacter septum area vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with 2D septum area model', 'Fitting with Faingold model','Caulobacter during division')
xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum area - using w fit','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off

% 1b) Anna: fit w^2 vs t

% [tc tg1 wMax1 fitOutTime1NoBinNoLeng]  = fitFeingold(t, w);
% [tg wMax fitOutTimeNoBinNoLeng]  = fitFeingoldFixedStart(t, w, tc);

[tc tg wMax fitOutAreaTimeNoBinNoLeng]  = fitAreaSeamusModBound(t, w.^2);
[tcF tgF wMaxF fitOutAreaFeingold]  = fitAreaFeingoldModBound(t, w.^2);

fitVariable.fitOutTimeNoBinNoLeng = fitOutAreaTimeNoBinNoLeng;
fitVariable.fitOutTimeNoBinNoLengTc = tc;
fitVariable.fitOutTimeNoBinNoLengTg = tg;
fitVariable.fitOutTimeNoBinNoLengWmax = wMax;


variableForCompW2FtsI.tm=t;
variableForCompW2FtsI.wm=w;
variableForCompW2FtsI.fitWvsTm = fitOutAreaTimeNoBinNoLeng;
variableForCompW2FtsI.fitWvsTparTCm = tc;
variableForCompW2FtsI.fitWvsTparTGm = tg;
variableForCompW2FtsI.fitWvsTparWmaxm = wMax;
variableForCompW2FtsI.fitWvsTFm = fitOutAreaFeingold;
variableForCompW2FtsI.fitWvsTparTCFm = tcF;
variableForCompW2FtsI.fitWvsTparTGFm = tgF;
variableForCompW2FtsI.fitWvsTparWmaxm = wMaxF;
variableForCompW2FtsI.binDataTWm = twBinsCentered;
variableForCompW2FtsI.binDataWTm = wtBins';


figure
plot(fitOutAreaTimeNoBinNoLeng(:, 1), fitOutAreaTimeNoBinNoLeng(:, 2), 'r')
hold on
plot(fitOutAreaFeingold(:, 1), fitOutAreaFeingold(:, 2), 'c')
plot(t, w.^2, 'bo')

% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')

text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax^2,'%.2f') '\left ( 1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right ) \right )$$'],'interpreter','latex', 'FontSize',12 )
text(300, 0.9, '$$ F_{model}(t) = wMax^2 * \left (  1 - \left (  \frac{t - t_c }{ t_g - t_c} \right ) \right )$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Caulobacter septum waist ratio w^2 vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with 2D septum area model', 'Fitting with Faingold model','Caulobacter during division')
xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum w^2','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off




%-----------------------------------------------------------------------------------
% Area plot 2
%------------------------------------------------------------------------------

figure
plot(wMax^2.*( (fitOutAreaTimeNoBinNoLeng(:, 1)-tc)./(tg-tc) ), wMax^2 - fitOutAreaTimeNoBinNoLeng(:, 2), 'r')
%plot( fitOutTimeNoBinNoLeng(:, 1), wMax^2 - fitOutTimeNoBinNoLeng(:, 2), 'r')
hold on
plot(wMaxF^2.*( (fitOutAreaFeingold(:, 1)-tcF)./(tgF -tcF) ) , wMaxF^2 - fitOutAreaFeingold(:, 2), 'c')
%plot( (fitOutFeingold(:, 1) ), wMaxF^2 - fitOutFeingold(:, 2), 'c')
plot( wMax^2.*(t-tc)./(tg-tc), wMax^2 - w.^2, 'bo')
%plot( t, wMax^2 - w.^2, 'bo')
title('Caulobacter septum area vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with 2D septum area model', 'Fitting with Faingold model','Caulobacter during division')
xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum area - using w^2 fit','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off


% 2) Anna: adding fitting using seamus functions (wPriv vs time computed from length)

% [tc tg1 wMax1 fitOutLength1]  = fitFeingoldForLength( lPriv, wPriv, timeLegthMod);
% [tg wMax fitOutLength]  = fitFeingoldFixedStartForLength( lPriv, wPriv, tc, timeLegthMod);

% [tc tg wMax fitOutLength]  = fitSeamusModBound(timeLegthMod(lPriv), wPriv);
% 
% fitVariable.fitOutLength = fitOutLength;
% fitVariable.fitOutLengthTc = tc;
% fitVariable.fitOutLengthTg = tg;
% fitVariable.fitOutLengthWmax = wMax;
% 
% figure
% plot(fitOutLength(:, 1), fitOutLength(:, 2), 'g')
% hold on
% time = timeLegthMod(lPriv);
% plot(time, wPriv, 'bo')
% 
% % text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% % text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% % set(gca,'box','on')
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time (computed from length)', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with 2D septum area model', 'Caulobacter during division')
% xlabel('Time [min] - Computed from length','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );

% 3) Anna: adding fitting w vs t computed form length

[tc tg wMax fitOutTimeFromL] = fitSeamusModBound(timeFromLengAllCellNoNeg, wAllCellNoNeg);

fitVariable.fitOutTimeBinWandL = fitOutTimeFromL;
fitVariable.fitOutTimeBinWandLTc = tc;
fitVariable.fitOutTimeBinWandLTg = tg;
fitVariable.fitOutTimeBinWandLWmax = wMax;

variableForCompm.tFromLm=timeFromLengAllCellNoNeg;
variableForCompm.wForTfromLm=wAllCellNoNeg;
variableForCompm.fitWvsTfromLm = fitOutTimeFromL;
variableForCompm.fitWvsTfromLparTCm = tc;
variableForCompm.fitWvsTfromLparTGm = tg;
variableForCompm.fitWvsTfromLparWmaxm = wMax;
variableForCompm.binDataTfromLm = tBinsAllCellFromLengFinNoNeg;
variableForCompm.binDataWforTfromLm = wBinsAllCellFinNoNeg;

figure
plot(fitOutTimeFromL(:, 1), fitOutTimeFromL(:, 2), 'r')
hold on
plot(timeFromLengAllCellNoNeg, wAllCellNoNeg, 'bo')
plot(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg, 'g*')

text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )}$$'],'interpreter','latex', 'FontSize',12 )
text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )}$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Caulobacter septum waist ratio r/R vs time - After length and waist ratio binning - All cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with 2D septum area model', 'Caulobacter during division', 'Binned data')
xlabel('Time (from l) [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off

% 3b) Anna: adding fitting w^2 vs t computed form length

[tc tg wMax fitOutTimeFromL] = fitAreaSeamusModBound(timeFromLengAllCellNoNeg, wAllCellNoNeg.^2);

fitVariable.fitOutTimeBinWandL = fitOutTimeFromL;
fitVariable.fitOutTimeBinWandLTc = tc;
fitVariable.fitOutTimeBinWandLTg = tg;
fitVariable.fitOutTimeBinWandLWmax = wMax;

variableForCompW2FtsI.tFromLm=timeFromLengAllCellNoNeg;
variableForCompW2FtsI.wForTfromLm=wAllCellNoNeg;
variableForCompW2FtsI.fitWvsTfromLm = fitOutTimeFromL;
variableForCompW2FtsI.fitWvsTfromLparTCm = tc;
variableForCompW2FtsI.fitWvsTfromLparTGm = tg;
variableForCompW2FtsI.fitWvsTfromLparWmaxm = wMax;
variableForCompW2FtsI.binDataTfromLm = tBinsAllCellFromLengFinNoNeg;
variableForCompW2FtsI.binDataWforTfromLm = wBinsAllCellFinNoNeg;

figure
plot(fitOutTimeFromL(:, 1), fitOutTimeFromL(:, 2), 'r')
hold on
plot(timeFromLengAllCellNoNeg, wAllCellNoNeg.^2, 'bo')


text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax^2,'%.2f') '\left ( 1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right ) \right )$$'],'interpreter','latex', 'FontSize',12 )
text(300, 0.9, '$$ F_{model}(t) = wMax^2 * \left (  1 - \left (  \frac{t - t_c }{ t_g - t_c} \right ) \right )$$' ,'interpreter','latex', 'FontSize',12 )
set(gca,'box','on')

title('Caulobacter septum waist ratio w^2 vs time - After length and waist ratio binning - All cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
legend ('Fitting with 2D septum area model', 'Caulobacter during division')
xlabel('Time (from l) [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
ylabel('Septum w^2','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
hold off


% 3c) Anna: adding fitting w vs t after length and w binning (all cell)

% [tc tg1 wMax fitOutTime1BinWandL]  = fitFeingold(tBinsAllCellCenteredFromLeng, wBinsAllCell);
% [tg wMax fitOutTimeBinWandL]  = fitFeingoldFixedStart(tBinsAllCellCenteredFromLeng, wBinsAllCell, tc);

% [tc tg wMax fitOutTimeBinWandL] = fitSeamusModBound(tBinsAllCellCenteredFromLeng, wBinsAllCell);
% 
% fitVariable.fitOutTimeBinWandL = fitOutTimeBinWandL;
% fitVariable.fitOutTimeBinWandLTc = tc;
% fitVariable.fitOutTimeBinWandLTg = tg;
% fitVariable.fitOutTimeBinWandLWmax = wMax;
% 
% figure
% plot(fitOutTimeBinWandL(:, 1), fitOutTimeBinWandL(:, 2), 'r')
% hold on
% plot(tBinsAllCellCenteredFromLeng, wBinsAllCell, 'bo')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time - After length and waist ratio binning - All cells', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division', 'Binned data')
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off


% 4) Anna: adding fitting w vs t after only w binning (all cell) 

% [tc tg1 wMax fitOutTime1BinW]  = fitFeingold(t, wAfterBinningAllCellNoTimeFromLeng);
% [tg wMax fitOutTimeBinW]  = fitFeingoldFixedStart(t, wAfterBinningAllCellNoTimeFromLeng, tc);
% 
% fitVariable.fitOutTimeBinW = fitOutTimeBinW;
% fitVariable.fitOutTimeBinWTc = tc;
% fitVariable.fitOutTimeBinWTg = tg;
% fitVariable.fitOutTimeBinWWmax = wMax;
% 
% figure
% plot(fitOutTimeBinW(:, 1), fitOutTimeBinW(:, 2), 'r')
% hold on
% plot(t, wAfterBinningAllCellNoTimeFromLeng, 'bo')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time - After only waist ratio binning - All cells ', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off

% 5) Anna: adding fitting w vs t after length and w binning (no negative time)
% 

% % Do not include cell with negative time
% negTimeIdx = find(timeFromLengAllCellAfterBinning<0);
% timeFromLengAfterBinningNoNeg = timeFromLengAfterBinning;
% timeFromLengAfterBinningNoNeg(negTimeIdx) = [];
% wAfterBinningFinNoNeg = wAfterBinningFin; 
% wAfterBinningFinNoNeg(negTimeIdx) = [];

% [tc tg1 wMax fitOutTime1BinWandLNoNegNoPost]  = fitFeingold(timeFromLengAfterBinningNoNeg, wAfterBinningFinNoNeg);
% [tg wMax fitOutTimeBinWandLNoNegNoPost]  = fitFeingoldFixedStart(timeFromLengAfterBinningNoNeg, wAfterBinningFinNoNeg, tc);

% [tc tg1 wMax fitOutTime1BinWandLNoNegAllCell]  = fitFeingoldWeighted(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg, errwBinsAllCellNoNeg);
% [tg wMax fitOutTimeBinWandLNoNegAllCell]  = fitFeingoldFixedStartWeighted(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg, tc, errwBinsAllCellNoNeg);

% [tc tg1 wMax fitOutTime1BinWandLNoNegAllCell]  = fitFeingold(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg);
% [tg wMax fitOutTimeBinWandLNoNegAllCell]  = fitFeingoldFixedStart(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg, tc);

% [tc tg wMax fitOutTimeBinWandLNoNegAllCell] = fitSeamusModBound(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg);
% 
% 
% fitVariable.fitOutTimeBinWandLNoNegAllCell = fitOutTimeBinWandLNoNegAllCell;
% fitVariable.fitOutTimeBinWandLNoNegNoPostTc = tc;
% fitVariable.fitOutTimeBinWandLNoNegNoPostTg = tg;
% fitVariable.fitOutTimeBinWandLNoNegNoPostWmax = wMax;
% fitVariable.errwBinsAllCellNoNeg = errwBinsAllCellNoNeg;
% 
% figure
% plot(fitOutTimeBinWandLNoNegAllCell(:, 1), fitOutTimeBinWandLNoNegAllCell(:, 2), 'r')
% hold on
% plot(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg, 'bo')
% errorbar(tBinsAllCellFromLengFinNoNeg, wBinsAllCellFinNoNeg, errwBinsAllCellNoNeg ,'b*')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time - After length and waist ratio binning - No negative time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off


% 6) Anna:  adding fitting w vs t after length and w binning (no cells
% after division) (no negative time)(no post divisional cells)
% [tc tg wMax fitOutTime1BinWandLNoCellAftNoNeg]  = fitFeingold(tBinsFromLengFinNoNeg, wBinsFinNoNeg );
% [tg wMax fitOutTimeBinWandLNoCellAftNoNeg]  = fitFeingoldFixedStart(tBinsFromLengFinNoNeg, wBinsFinNoNeg , tc);
% 
% fitVariable.fitOutTimeBinWandLNoCellAftNoNeg = fitOutTimeBinWandLNoCellAftNoNeg;
% fitVariable.fitOutTimeBinWandLNoCellAftTc = tc;
% fitVariable.fitOutTimeBinWandLNoCellAftTg = tg;
% fitVariable.fitOutTimeBinWandLNoCellAftWmax = wMax;
% 
% figure
% plot(fitOutTimeBinWandLNoCellAftNoNeg(:, 1), fitOutTimeBinWandLNoCellAftNoNeg(:, 2), 'r')
% hold on
% plot(tBinsFromLengFinNoNeg, wBinsFinNoNeg , 'bo')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time - After length and waist ratio binning - No post div cells - No neg time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off

% 7) Anna:  adding fitting w vs t after only length binning (no cells after division)
% [tc tg wMax fitOutTimeOnlyLengthBin]  = fitFeingold(timeFromLengAfterBinning, wPriv );
% [tg wMax fitOutTimeOnlyLengthBin]  = fitFeingoldFixedStart(timeFromLengAfterBinning, wPriv , tc);
% 
% fitVariable.fitOutTimeOnlyLengthBin = fitOutTimeOnlyLengthBin;
% fitVariable.fitOutTimeOnlyLengthBinTc = tc;
% fitVariable.fitOutTimeOnlyLengthBinTg = tg;
% fitVariable.fitOutTimeOnlyLengthBinWmax = wMax;
% 
% figure
% plot(fitOutTimeOnlyLengthBin(:, 1), fitOutTimeOnlyLengthBin(:, 2), 'r')
% hold on
% plot(timeFromLengAfterBinning, wPriv , 'bo')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time - After only length binning', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off

% 8) Anna: adding fitting wPriv vs time computed from length ( no negative time)

% Do not include cell with negative time
% negTimeIdx = find(timeFromLeng<0);
% timeFromLengNoNeg = timeFromLeng;
% timeFromLengNoNeg(negTimeIdx) = [];
% wPrivNoNeg = wPriv; 
% wPrivNoNeg(negTimeIdx) = [];

% [tc tg1 wMax1 fitOutLength1NoNeg]  = fitFeingold( timeFromLengNoNeg, wPrivNoNeg);
% [tg wMax fitOutLengthNoNeg]  = fitFeingoldFixedStart( timeFromLengNoNeg, wPrivNoNeg, tc);

% [tc tg wMax fitOutLengthNoNeg] = fitSeamusModBound(timeFromLengNoNeg, wPrivNoNeg);
% 
% 
% fitVariable.fitOutLengthNoNeg = fitOutLengthNoNeg;
% fitVariable.fitOutLengthNoNegTc = tc;
% fitVariable.fitOutLengthNoNegTg = tg;
% fitVariable.fitOutLengthNoNegWmax = wMax;
% 
% 
% figure
% plot(fitOutLengthNoNeg(:, 1), fitOutLengthNoNeg(:, 2), 'g')
% hold on
% plot(timeFromLengNoNeg, wPrivNoNeg, 'bo')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time (computed from length) - No neg time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% xlabel('Time [min] - Computed from length','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% 
% 
% variableForComp = cell(1, 7);
% variableForComp.fitOutTimeBinWandLNoNegAllCell = fitOutTimeBinWandLNoNegAllCell;
% variableForComp.fitOutTimeBinWandLNoNegNoPostTc = tc;
% variableForComp.fitOutTimeBinWandLNoNegNoPostTg = tg;
% variableForComp.fitOutTimeBinWandLNoNegNoPostWmax = wMax;
% variableForComp.errwBinsAllCellNoNeg = errwBinsAllCellNoNeg;
% 
% variableForComp.tBinsAllCellFromLengFinNoNeg = tBinsAllCellFromLengFinNoNeg;
% variableForComp.wBinsAllCellFinNoNeg = wBinsAllCellFinNoNeg;

% Anna: fit area
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0 0],... % lower boundary paramters conditions
%                'Upper',[Inf Inf],... % upper boundary paramters conditions
%                'StartPoint',[1 1]);
% ft = fittype('a*(x)^n+b','problem','n','options',fo);
% 
% [fitRes, gof2] = fit(t, a, ft, 'problem', 2)
% 
% figure,
% plot( fitRes,'m')
% hold on
% plot(t, a, 'b*')
% 
% legend('Data','n=2')
% hold off


%%try a contour plot
%t= [tExpt(cellOk)';tDiv(:)];
%w= [W(cellOk)';dDiv(:)];
%tCtr = 0:20:max(t);
%wCtr = 0:0.05:1.2;
%ctrs={tCtr,wCtr}
%n=hist3([t,w],'Ctrs',ctrs);
%NMAX = 20;
%n(n>NMAX)=NMAX;
%figure;
%nCountour=16;
%h = contourf(tCtr,wCtr,n');
%colormap('gray');
%cmap = colormap;
%cmap = flipud(cmap);
%colormap(cmap);
%ylim([0 1.2]);
%ylabel('Waist ratio r/R');
%xlabel('T (min)');
%
%figure;
%h = pcolor(tCtr,wCtr,n');
%set(h,'edgecolor','none');
%colormap('gray');
%cmap = colormap;
%cmap = flipud(cmap);
%colormap(cmap);
%ylim([0 1.2]);
%ylabel('Waist ratio r/R');
%xlabel('T (min)');
%

%try an area plot

% Anna: diff(w) vs time computed from length
% [tc tg wMax fitOut]  = fitFeingoldForDeriv( l(1: length(l)-1), diff(w), timeLegthMod);
% [tg wMax fitOut]  = fitFeingoldFixedStartForDeriv( l(1: length(l)-1), diff(w), tc, timeLegthMod);
% figure
% plot(fitOut(:, 1), fitOut(:, 2), 'g')
% hold on
% time = timeLegthMod(l);
% plot(time(1:length(time)-1), diff(w), 'bo')
% 
% text(350, -0.1,[  num2str(wMax,'%.2f') ['$$\sqrt{(1 - ( ( t - ' num2str(round(tc)) ' )/(' num2str(round(tg)) ' - '  num2str(round(tc)) ') )^2)^{-1}}*1/( ' num2str(round(tg)) ' - ' num2str(round(tc)) ')$$'] ],'interpreter','latex', 'FontSize',12 )
% text(350, -0.3, '$$ F_{model}(t) = wMax*\sqrt{(1 -  ( (t - t_c )/( t_g - t_c) )^2 )^{-1}}*1/(t_g - t_c)$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Variation of septum diameters respect to the length', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% title('Caulobacter septum diameter derivative vs time', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length) derivative','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off

%  4) Anna: adding fitting w vs t after length and w binning (all cell) (no negative time)
% 
% % Do not include cell with negative time
% negTimeIdx = find(timeFromLengAllCellAfterBinning<0);
% timeFromLengAllCellAfterBinningNoNeg = timeFromLengAllCellAfterBinning;
% timeFromLengAllCellAfterBinningNoNeg(negTimeIdx) = [];
% wAfterBinningAllCellFinNoNeg = wAfterBinningAllCellFin; 
% wAfterBinningAllCellFinNoNeg(negTimeIdx) = [];

% [tc tg1 wMax fitOutTime1BinWandLNoNeg]  = fitFeingold(timeFromLengAllCellAfterBinningNoNeg, wAfterBinningAllCellFinNoNeg);
% [tg wMax fitOutTimeBinWandLNoNeg]  = fitFeingoldFixedStart(timeFromLengAllCellAfterBinningNoNeg, wAfterBinningAllCellFinNoNeg, tc);

% fitVariable.fitOutTimeBinWandLNoNeg = fitOutTimeBinWandLNoNeg;
% fitVariable.fitOutTimeBinWandLNoNegTc = tc;
% fitVariable.fitOutTimeBinWandLNoNegTg = tg;
% fitVariable.fitOutTimeBinWandLNoNegWmax = wMax;
% 
% figure
% plot(fitOutTimeBinWandLNoNeg(:, 1), fitOutTimeBinWandLNoNeg(:, 2), 'r')
% hold on
% plot(timeFromLengAllCellAfterBinningNoNeg,wAfterBinningAllCellFinNoNeg, 'bo')
% 
% text(300, 0.7,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
% text(300, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
% set(gca,'box','on')
% 
% title('Caulobacter septum waist ratio r/R vs time - After length and waist ratio binning - All cells - No negative time cell', 'FontSize', 20, 'FontName', 'Copperplate Gothic Light'),
% legend ('Fitting with Feingold model', 'Caulobacter during division')
% xlabel('Time [min]','FontSize',16, 'FontName', 'Copperplate Gothic Light'  );
% ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',16, 'FontName', 'Copperplate Gothic Light' );
% hold off
%-------------------------------------------------------------
function [xOut,yOut]=spline2d(x,y,nBreak);
nPt = numel(x);
nPtOut = 10*nPt;% upsample the output, but keep discrete

t=linspace(0,1,nPt);
tOut=linspace(0,1,nPtOut);
ppX = splinefit(t,x,nBreak);
ppY = splinefit(t,y,nBreak);
xOut = ppval(ppX,tOut);
yOut = ppval(ppY,tOut);
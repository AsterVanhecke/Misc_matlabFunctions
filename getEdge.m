function [diamEdge, edges, width, xLeft, xRight,yLeft,yRight]= getEdge(x,y,xLim,searchStep, slicePos,xCtr)


    % x and y are the histogram points coordinates
    %find first maximum on left
    [posLeft,posLeftMax,idxLeftMax, valLeftMax] = getLeft(x,y,xLim(1),searchStep);
    [posRight,posRightMax,idxRightMax,valRightMax] = getRight(x,y,xLim(2),searchStep);
    
% check if the are almost 1 localization
if isempty(x) || isempty(y) && ii==1
    diamEdge = abs(xLim(2)- xLim(1));
    edges=[xLim(1) xLim(2)];
    
elseif isempty(x) || isempty(y) || length(x)<6|| length(y)<6
    diamEdge = abs(xLim(2)- xLim(1));
    edges=[xLim(1) xLim(2)];
    
    % check if there are localization coordinates that are NaN
elseif sum(isnan(x))>1 || sum(isnan(y))>1
    diamEdge = abs(xLim(2)- xLim(1));
    edges=[xLim(1) xLim(2)];
    
elseif length(x)>5
    if isempty(idxLeftMax) || isempty(idxRightMax)
        xxLeft = x(1:round( (length(x)-1)/2 ) );
        yyLeft = y(round( (length(y)-1)/2 ):end);
        xxRight = x(idxRightMax:end);
        yyRight = y(idxRightMax:end);
        
%         [idxLeftMaxDeriv, posLeftMaxDeriv] = getEdgeLoc(xxLeft, yyLeft, 1);
%         [idxRightMaxDeriv, posRightMaxDeriv] = getEdgeLoc(xxRight, yyRight, -1);

        [idxLeftMaxDeriv, posLeftMaxDeriv] = getEdgeLoc(x, y, 1);
        [idxRightMaxDeriv, posRightMaxDeriv] = getEdgeLoc(x, y, -1);
   
    else
        if idxLeftMax==1
            xxLeft = x(1:round(length(x)/2));
            yyLeft = y(1:round(length(x)/2));
            xxRight = x(round(length(x)/2):end);
            yyRight = y(round(length(x)/2):end);
        elseif idxRightMax==length(x)
            xxLeft = x(1:round(length(x)/2));
            yyLeft = y(1:round(length(x)/2));
            xxRight = x(round(length(x)/2):end);
            yyRight = y(round(length(x)/2):end);
        else
        xxLeft = x(1:idxLeftMax);
        yyLeft = y(1:idxLeftMax);
        xxRight = x(idxRightMax:end);
        yyRight = y(idxRightMax:end);
        end
%         [idxLeftMaxDeriv, posLeftMaxDeriv, xSmooLeft, ySmooLeft] = getEdgeLoc(xxLeft, yyLeft, 1);
%         [idxRightMaxDeriv, posRightMaxDeriv, xSmooRight, ySmooRight] = getEdgeLoc(xxRight, yyRight, -1);

%         [idxLeftMaxDeriv, posLeftMaxDeriv, xSmooLeft, ySmooLeft] = getEdgeLoc(x(1:round(end/3)), y(1:round(end/3)), 1);
%         [idxRightMaxDeriv, posRightMaxDeriv, xSmooRight, ySmooRight] = getEdgeLoc(x(round(2*end/3):end), y(round(2*end/3):end), -1);
        
        [idxLeftMaxDeriv, posLeftMaxDeriv, xSmooLeft, ySmooLeft] = getEdgeLoc(x, y, 1);
        [idxRightMaxDeriv, posRightMaxDeriv, xSmooRight, ySmooRight] = getEdgeLoc(x, y, -1);
    
        
    % DEBUG
%     
%     hist(slicePos,xCtr);
%     hold all
%     plot(x,y, 'r--');
%     plot(xSmooLeft,ySmooLeft, 'm--');
%     plot(xSmooRight,ySmooRight, 'm--');
%     plot(xSmooLeft(idxLeftMaxDeriv) ,ySmooLeft(idxLeftMaxDeriv), 'g*')
%     plot(xSmooRight(idxRightMaxDeriv) ,ySmooRight(idxRightMaxDeriv), 'g*')
%     drawnow
%     hold off
%     pause
%         
    end

    edges = [posLeftMaxDeriv posRightMaxDeriv];
    diamEdge = posRightMaxDeriv - posLeftMaxDeriv;
end
  %%DEBUG
%     hold off;
%     plot(x,y);
%     hold all;
%     plot(posLeftMax,valLeftMax,'o');
%     plot(posLeft,valLeftMax/2,'o')
%     plot(posRightMax,valRightMax,'o');
%     plot(posRight,valRightMax/2,'o');
    %%DEBUG

    width = posRight-posLeft;
    xLeft = posLeft;
    xRight = posRight;
    yLeft =valLeftMax/2;
    yRight =valRightMax/2;
    

    



%---------------------------------------------------------
function [idxMaxDeriv, posMaxDeriv, xx, yy] = getEdgeLoc(x, y, signe)
 
    fitMod = griddedInterpolant(x, y, 'spline');
    % not more than 4nm precision expected
    xx = (min(x):4:max(x));
    yi = fitMod(xx);

%     fitPar = polyfit(x, y, 4);
%     % not more than 4nm precision expected
%     xx = (min(x):4:max(x)-4);
%     yy = polyval(fitPar, xx);
%     
   yy = smooth(xx, yi, 'rloess', 4);
%   [yy, xx] = ksdensity(slicePos,'width', 0.8);
    % Compute the first and second derivative
   [dfdxTemp, d2fdx2Temp] = dfdxc(xx,yy');
   if signe==1
       dfdxPivot = dfdxTemp(3:round(end/3));
       d2fdx2Pivot = d2fdx2Temp(3:round(end/3));
       minFoundIdxTemp = find(dfdxPivot<0);
       if isempty(minFoundIdxTemp)
            minFoundIdx = length(dfdxPivot);
       else
            minFoundIdx = min(minFoundIdxTemp);
       end
       dfdx = dfdxPivot(1:minFoundIdx);
       d2fdx2 = d2fdx2Pivot(1:minFoundIdx);
      
   elseif signe==-1
       dfdxPivot = dfdxTemp(round(2*end/3):end-3);
       d2fdx2Pivot = d2fdx2Temp(round(2*end/3):end-3);
       minFoundIdxTemp = find(dfdxPivot>0);
       if isempty(minFoundIdxTemp)
            minFoundIdx = 1;
       else
            minFoundIdx = max(minFoundIdxTemp);
       end
       dfdx = dfdxPivot(minFoundIdx:end);
       d2fdx2 = d2fdx2Pivot(minFoundIdx:end);
       
   end
   stdDer = std(dfdx);
   meanDer = mean(dfdx);
   maxRangeDerIdx = find(dfdx >= meanDer  & dfdx < meanDer + 1*stdDer); 
   minRangeDerIdx = find(dfdx <= meanDer  & dfdx > meanDer - 1*stdDer); 
   if signe==1
       if isempty(maxRangeDerIdx)
           
           if length(find(max(dfdx)))==1
               idxMaxDeriv = find(max(dfdx));
           elseif isempty ( find(max(dfdx)) )
               idxMaxDeriv = 2;
           else   
               idxMaxTemp = find(max(dfdx));
               idxMaxDeriv = idxMaxTemp(2);    
           end
           
       elseif sum( dfdx == max(dfdx(maxRangeDerIdx)) )==1
           idxMaxDeriv = find( dfdx == max(dfdx(maxRangeDerIdx)));
       else   
           idxTemp = find( dfdx == max(dfdx(maxRangeDerIdx)) );
           idxMaxDeriv = idxTemp(2);
       end
        idxMaxDeriv =  idxMaxDeriv  + 3 ;
       % dfdxTemp(idxMaxDeriv)
   elseif signe==-1
       if isempty(minRangeDerIdx)
           if length(find(min(dfdx)))==1
                idxMaxDeriv = find(min(dfdx));
           elseif isempty ( find(min(dfdx)) )
               idxMaxDeriv = length(x)-1;
           else
                idxMaxTemp = find(min(dfdx));
                idxMaxDeriv = idxMaxTemp(end-1);
          
           end
       elseif sum( dfdx == min(dfdx(minRangeDerIdx)) )==1
           idxMaxDeriv = find( dfdx == min(dfdx(minRangeDerIdx)) );
       else
           idxTemp = find( dfdx == min(dfdx(minRangeDerIdx)) );
           idxMaxDeriv = idxTemp(end);
       end
        idxMaxDeriv =  idxMaxDeriv  + round(2*length(dfdxTemp)/3) + minFoundIdx -1;
        % dfdxTemp(idxMaxDeriv)
   end
   posMaxDeriv = xx(idxMaxDeriv);
   
   % DEBUG
%    figure,
%    plot(x, y, 'r--')
%    hold all
%    plot(xx,yy, 'k-')
%    plot(xx(idxMaxDeriv) ,yy(idxMaxDeriv), 'r*' )
%    pause




%---------------------------------------------------------
function [posLeft,posLeftMax,idxLeftMax, valLeftMax] = getLeft(x,y,xLim,searchStep);
nX = numel(x);
ii =2;
foundLeftPk = false;
idxLeftMax=[];
valLeftMax=[];
posLeftMax=[];
while ~foundLeftPk && ii<=nX-1
  if (y(ii)-y(ii-1)>0)&&(y(ii+1)-y(ii)<=0)
    idxLeftMax = ii;
    valLeftMax = y(ii);
    posLeftMax=x(ii);
    foundLeftPk=true;
  elseif (y(ii)-y(ii-1)<0)%if we're descending we're past the peak
    idxLeftMax = 1;
    valLeftMax = y(1);
    posLeftMax=x(1);
    foundLeftPk=true;
  else
    ii=ii+1;
  end
end

if isempty(idxLeftMax)
  posLeft = NaN;
  posLeftMax=NaN;
  valLeftMax=NaN;
else
  %find the half max value
  posLeft = x(idxLeftMax);
  valCur = Inf;
  while posLeft>=x(1) && valCur > valLeftMax/2
    posLeft =posLeft-searchStep;
    valCur = interp1(x,y,posLeft);
  end
end

%---------------------------------------------------------
function [posRight,posRightMax,idxRightMax,valRightMax] = getRight(x,y,xLim,searchStep);
nX = numel(x);
ii =nX-1;
foundRightPk = false;
idxRightMax=[];
valRightMax=[];
posRightMax=[];
while ~foundRightPk && ii>=2
  if (y(ii)-y(ii-1)>=0)&&(y(ii+1)-y(ii)<0)
    idxRightMax = ii;
    valRightMax = y(ii);
    posRightMax=x(ii);
    foundRightPk=true;
  elseif (y(ii)-y(ii-1)>0)%if we're descending we're past the peak
    idxRightMax =nX;
    valRightMax = y(nX);
    posRightMax=x(nX);
    foundRightPk=true;
  else
    ii=ii-1;
  end
end

if isempty(idxRightMax)
  posRight = NaN;
  posRightMax=NaN;
  valRightMax=NaN;
else
  %find the half max value
  posRight = x(idxRightMax);
  valCur = Inf;
  while posRight<=x(end) && valCur > valRightMax/2
    posRight =posRight+searchStep;
    valCur = interp1(x,y,posRight);
  end
end


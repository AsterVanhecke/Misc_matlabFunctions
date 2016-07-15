finished=false;
cellIdx=0;
cellMat=[];
polygonia=cell(1,1);
while ~finished
    % get file path and name, read rapidstorm file
    disp('Choose a file')
    [FileName,PathName] = uigetfile;
    if ~FileName
        finished=true;
    else
        locs=loadLoc([PathName, FileName],'rapidStorm');
        x=locs(:,1);
        y=locs(:,3);
        figure, plot(x,y,'.','MarkerSize',1)
        nCells=input('How many cells do you see in this FOV?');
        for cellFOV=1:nCells
            cellIdx=cellIdx+1;
            [xPoly, yPoly]=ginput;
            IN=inpolygon(x, y, xPoly, yPoly);
            hold on
            plot(xPoly,yPoly)
            Area=polyarea(xPoly,yPoly);
            inNr=sum(IN);
            dens=inNr/Area;

    %         fprintf('grouping time!') %let's ignore the grouping for now
    %         % do grouping, try to get get indices of unique emitters, evaluate
    %         % if they are inpolygon, count.

            cellMat(cellIdx,:)=[inNr, dens]
            polygonia{cellIdx,1}=[xPoly; yPoly];
        end

        fprintf('please open the next file.\n')
        fprintf('Ok, lets go on.\n')
    end
end

%%

figure,
hist(cellMatY(:,2))
figure,
hist(cellMatT(:,2))
figure,
hist(cellMatE(:,2))

figure,
hist(cellMatY(:,1))
figure,
hist(cellMatT(:,1))
figure,
hist(cellMatE(:,1))


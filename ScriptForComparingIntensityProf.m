% ----------------------------------------------------------------------%
% Script to compare intensity profile of the cells with different length
%-----------------------------------------------------------------------%

% M is the matrix with lengths and corresponding intensities as coulums for
% each cell

for vectIdx = 1 : round(size(M,2)/2)-1
    i=vectIdx;
    vectLeng(i, 1) = i*2;
    vectLeng(i, 2) = nnz(M(:, i*2));
end

orderedVectLeng = sortrows(vectLeng, 2);

maxIdx = length(orderedVectLeng);
maxLeng = orderedVectLeng(maxIdx, 2);
intMatrix = zeros(maxIdx, maxLeng);

for rowIdx =1: maxIdx
    
    lengthCurrVect = orderedVectLeng(rowIdx, 2);
    midPointMaxLengIdx = round(maxLeng/2);
    midPointCurrentVect = round(orderedVectLeng(rowIdx, 2)/2);
    startingVectPosition = midPointMaxLengIdx - midPointCurrentVect;
    endingVectPosition = startingVectPosition + orderedVectLeng(rowIdx, 2);
    
    columnPosition = orderedVectLeng(rowIdx, 1);
    intMatrix(rowIdx, startingVectPosition+1 : endingVectPosition ) = ...
        M(1:lengthCurrVect, columnPosition);
    set(gca,'XTickLabel',[0 125 250 375 500 625 750 875 1000 1250 1375 1500 1625 1750 1875 2000 2125 2250 2375 2500 2625 2750 2875 3000 3125 3250 3375 3500 3625 3750 3875 4000 4125 4250 4325 4500 4625 4750 4875 5000 5125 5250 ] );
end

figure, imshow(intMatrix)

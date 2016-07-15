% Find the furthest away points(the pole).
            numberOfCoords = size(cCell,1);
            x = cCell(:,2);
            y = cCell(:,1);
            for k = 1 : numberOfCoords
              distances = sqrt((x-x(k)).^2 + (y-y(k)).^2);
              [maxDistance(k), indexOfMax(k)] = max(distances);
            end
            pole1PosIdx = find(maxDistance == max(maxDistance));
            pole2PosIdx = indexOfMax(pole1PosIdx);
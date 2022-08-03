function [idx, TWT_crop, data_crop] = adjustTWT(TWT, data)
    %Use this function when determining the depth of reflectors from the
    %TWT or just to adjust the TWT to begin 0 at the surface.

    % Author: Annika Horlings
    % University of Washington
    % Last updated: 16 June 2022
    
    [idx] = find(data(:, 1), 1, 'first'); %find index of first nonzero
    TWT_crop = TWT(idx:end) - TWT(idx); %crop and adjust to start at 0
    data_crop = data(idx:end, :); %crop
end
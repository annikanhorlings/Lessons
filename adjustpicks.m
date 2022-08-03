function [idx, picks_adj] = adjustpicks(picks, TWT, data)
    %Use this function if picked before adjusting the TWT to 0 at the
    %surface

    % Author: Annika Horlings
    % University of Washington
    % Last updated: 16 June 2022
    
    [idx] = find(data(:, 1), 1, 'first'); %find index of first nonzero
    picks_adj = picks - idx; %crop and adjust to start at 0;
        %note that the picks are the index of the two-way travel time
end
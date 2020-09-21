function [CELL_Activity, VIDEO_Activity] = PerCell(Table)
Table = load_burn_excel();
% Isolate Image ID and Cell ID
Image = Table.Image;
CELL_ID = Table.CELL_ID;
Soma_ROI = Table.Final_calcium_soma;
% Eliminate repeat Image ID's
[u_Image, idx] = unique(Image, 'Rows');
ImageID = num2str(u_Image);

% % For each Image
    % % For each Cell IID
        % % Look at suite2p Soma ROI
            % % Say if suite2p soma had event
        % % END
        % % Frequency of events, etc. for each cell
    % % END
    % % Percent of cells w/ active events
    % % Num. of cells per image
% % END

for ii = 1:length(u_Image)
    %%% Extract rows corresponding to a single video
    if ii < length(u_Image)
        subTable = Table(idx(ii):idx(ii+1)-1,:);
    else
        subTable = Table(size(Table,1)-idx(ii),:);
    end
    tempCELL_ID = subTable.CELL_ID;
    unique_cells = unique(tempCELL_ID);
    unique_cells = unique_cells(~isnan(unique_cells));
    tempACTIVE = subTable.Final_soma;
    tempACTIVE = tempACTIVE(~cellfun(@isempty,tempACTIVE));
    total_cells = numel(unique_cells);
    active_cells = numel(tempACTIVE);
    current_video = strcat('T',ImageID(ii,:));
    if total_cells ~= 0
        VIDEO_Activity.(current_video(:)).percent_active_cells = (active_cells/total_cells)*100;
    end
    for jj = 1:size(unique_cells,1)
        % Immediately check for NaN's in CELL_ID or if roi =/= soma
        if isnan(tempCELL_ID) == 1 | cellfun(@isempty,subTable.Final_soma)==1
            continue
        end
        
    end
end
        
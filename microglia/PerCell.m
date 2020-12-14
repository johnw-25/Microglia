function [CELL_Activity, VIDEO_Activity] = PerCell(ProcessedEvents)
Table = load_burn_excel();
Table = sortrows(Table, 'Image');
% Isolate Image ID and Cell ID
Image = Table.Image;
CELL_ID = Table.CELL_ID;
% Eliminate repeat Image ID's
[u_Image, idx] = unique(Image, 'Rows');
idx = sort(idx);
u_Image = sort(u_Image);
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
eventTracks = fieldnames(ProcessedEvents);
for ii = 1:length(u_Image)
    %%% Extract rows corresponding to a single video
    if ii < length(u_Image)
        subTable = Table(idx(ii):idx(ii+1)-1,:);
    else
        subTable = Table(size(Table,1)-idx(ii),:);
    end
    tempCELL_ID = subTable.CELL_ID;
    tempCellCell = tempCELL_ID(~isnan(tempCELL_ID));
    tempTracks = subTable.TrackID;
    tempTrackCells = tempTracks(~isnan(tempCELL_ID));
    currentVideo = strcat('T',ImageID(ii,:));
    for i = 1:length(tempTrackCells)
        cellField = strcat('C',num2str(tempCellCell(i)));
        for j = 1:length(eventTracks)
            if strcmp(tempTracks(i),eventTracks(j))
                tempEvents = ProcessedEvents.(eventTracks{j}).PeakLocs;
                tempEvents = tempEvents(~isnan(tempEvents));
                try
                    if isfield(CELL_Activity.(currentVideo), cellField)
                        tempEvents = numel(tempEvents) +  CELL_Activity.(currentVideo).(cellField).eventCount;
                        CELL_Activity.(currentVideo).(cellField).eventCount = tempEvents;
                    else
                        CELL_Activity.(currentVideo).(cellField).eventCount = numel(tempEvents);
                    end
                catch
                    CELL_Activity.(currentVideo).(cellField).eventCount = numel(tempEvents);
                end
            else
                try
                    if ~any(~cellfun('isempty',strfind(eventTracks,tempTracks(i)))) && ~isfield(CELL_Activity.(currentVideo), cellField)
                        CELL_Activity.(currentVideo).(cellField).eventCount = 0;
                    end
                catch

                end
            end
        end
    end
    unique_cells = unique(tempCELL_ID);
    unique_cells = unique_cells(~isnan(unique_cells));
    tempACTIVE = subTable.Final_soma;
    tempACTIVE = tempACTIVE(~cellfun(@isempty,tempACTIVE));
    total_cells = numel(unique_cells);
    active_cells = numel(tempACTIVE);
    
    if total_cells ~= 0
        VIDEO_Activity.(currentVideo(:)).percent_active_cells = (active_cells/total_cells)*100;
    end
    for jj = 1:size(unique_cells,1)
        % Immediately check for NaN's in CELL_ID or if roi =/= soma
        if any(isnan(tempCELL_ID)) || any(cellfun(@isempty,subTable.Final_soma)==1)
            continue
        end
        
    end
end

end
        
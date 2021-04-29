function [CELL_Activity] = PerCell(Table, ProcessedEvents)
if ~nargin
    error('No valid input arguments given.')
end
CELL_Activity = struct();
% Isolate Image ID and Cell ID
Image = Table.Image;
%CELL_ID = Table.CELL_ID;
% Eliminate repeat Image ID's
[u_Image, idx] = unique(Image, 'Rows');
idx = sort(idx);
u_Image = sort(u_Image);
ImageID = num2str(u_Image);

eventTracks = fieldnames(ProcessedEvents);
for ii = 1:length(u_Image)
    %%% Extract rows corresponding to a single video
    if ii < length(u_Image)
        subTable = Table(idx(ii):idx(ii+1)-1,:);
    else
        subTable = Table(size(Table,1)-idx(ii),:);
    end
    tempCELL_ID = subTable.CELL_ID;
    tempMice = subTable.Mouse;
    tempMice = tempMice(~isnan(tempMice));
    tempCellCell = tempCELL_ID(~isnan(tempCELL_ID));
    tempTracks = subTable.TrackID;
    tempTrackCells = tempTracks(~isnan(tempCELL_ID));
    currentVideo = strcat('T',ImageID(ii,:));
    timehitR = subTable.TimeHitRadius1;
    
    for i = 1:length(tempTrackCells)
        cellField = strcat('C',num2str(tempCellCell(i)));
        for j = 1:length(eventTracks)
            if strcmp(tempTracks(i),eventTracks(j))
                tempEvents = ProcessedEvents.(eventTracks{j}).PeakLocs;
                tempTimeHitR = ProcessedEvents.(eventTracks{j}).timehitR;
                tempEvents = tempEvents(~isnan(tempEvents));
                tempP1Events = tempEvents(tempEvents <= 90);
                tempP2Events = tempEvents((tempEvents > 90) & (tempEvents < tempTimeHitR.*0.97));
                tempP3Events = tempEvents(tempEvents >= round(timehitR(i).*0.97));
                try
                    if isfield(CELL_Activity.(currentVideo), cellField)
                        tempEvents = numel(tempEvents) +  CELL_Activity.(currentVideo).(cellField).eventCount;
                        tempP1Events = numel(tempP1Events) +  CELL_Activity.(currentVideo).(cellField).Phase1.eventCount;
                        tempP2Events = numel(tempP2Events) +  CELL_Activity.(currentVideo).(cellField).Phase2.eventCount;
                        tempP3Events = numel(tempP3Events) +  CELL_Activity.(currentVideo).(cellField).Phase3.eventCount;
                        CELL_Activity.(currentVideo).(cellField).eventCount = tempEvents;
                        CELL_Activity.(currentVideo).(cellField).MouseID = tempMice(1);
                        CELL_Activity.(currentVideo).(cellField).Phase1.eventCount = tempP1Events;
                        CELL_Activity.(currentVideo).(cellField).Phase2.eventCount = tempP2Events;
                        CELL_Activity.(currentVideo).(cellField).Phase3.eventCount = tempP3Events;
                        CELL_Activity.(currentVideo).(cellField).numProcesses = CELL_Activity.(currentVideo).(cellField).numProcesses + 1;
                    else
                        CELL_Activity.(currentVideo).(cellField).eventCount = numel(tempEvents);
                        CELL_Activity.(currentVideo).(cellField).Phase1.eventCount = numel(tempP1Events);
                        CELL_Activity.(currentVideo).(cellField).Phase2.eventCount = numel(tempP2Events);
                        CELL_Activity.(currentVideo).(cellField).Phase3.eventCount = numel(tempP3Events);
                        CELL_Activity.(currentVideo).(cellField).numProcesses = 1;
                        CELL_Activity.(currentVideo).(cellField).MouseID = tempMice(1);
                    end
                catch
                    CELL_Activity.(currentVideo).(cellField).eventCount = numel(tempEvents);
                    CELL_Activity.(currentVideo).(cellField).Phase1.eventCount = numel(tempP1Events);
                    CELL_Activity.(currentVideo).(cellField).Phase2.eventCount = numel(tempP2Events);
                    CELL_Activity.(currentVideo).(cellField).Phase3.eventCount = numel(tempP3Events);
                    CELL_Activity.(currentVideo).(cellField).numProcesses = 1;
                    CELL_Activity.(currentVideo).(cellField).MouseID = tempMice(1);
                end
            else
                try
                    if ~contains(eventTracks,tempTracks(i)) && ~isfield(CELL_Activity.(currentVideo), cellField)
                        CELL_Activity.(currentVideo).(cellField).eventCount = 0;
                        CELL_Activity.(currentVideo).(cellField).Phase1.eventCount = 0;
                        CELL_Activity.(currentVideo).(cellField).Phase2.eventCount = 0;
                        CELL_Activity.(currentVideo).(cellField).Phase3.eventCount = 0;
                        CELL_Activity.(currentVideo).(cellField).numProcesses = 0;
                        CELL_Activity.(currentVideo).(cellField).MouseID = tempMice(1);
                    end
                catch

                end
            end
        end
    end
end

end
        
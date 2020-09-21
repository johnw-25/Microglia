function sortedProcessedData = SortData(ProcessedEvents)
% % % This function takes the traces processed in processedevents and
% % % categorized them by their pathology and DPI

fn = fieldnames(ProcessedEvents);
for k = 1:numel(fn)
    if strcmp(ProcessedEvents.(fn{k}).Pathology, 'PBS')
        sortedProcessedData.PBS.(fn{k}) = ProcessedEvents.(fn{k});
    elseif strcmp(ProcessedEvents.(fn{k}).Pathology, 'TMEV')
        if ProcessedEvents.(fn{k}).DPI == 2
            sortedProcessedData.TMEV.DPI_2.(fn{k}) = ProcessedEvents.(fn{k});
        elseif 2 < ProcessedEvents.(fn{k}).DPI < 7
            sortedProcessedData.TMEV.DPI_5.(fn{k}) = ProcessedEvents.(fn{k});
        elseif ProcessedEvents.(fn{k}).DPI >= 14
            sortedProcessedData.TMEV.DPI_14.(fn{k}) = ProcessedEvents.(fn{k});
        else
            disp('Unidentified DPI, check TMEV group for abnormal DPI  times');
        end
    else
        disp('Unidentifiable pathology.');
    end
end
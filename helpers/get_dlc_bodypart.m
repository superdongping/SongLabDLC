function bp = get_dlc_bodypart(dlc, candidateNames)
%GET_DLC_BODYPART Return a DLC body part struct by trying several names.

if ischar(candidateNames) || isstring(candidateNames)
    candidateNames = cellstr(candidateNames);
end

available = string(dlc.bodyparts);
for i = 1:numel(candidateNames)
    query = string(candidateNames{i});
    idx = find(strcmpi(available, query), 1);
    if ~isempty(idx)
        fieldName = dlc.field_names{idx};
        bp = dlc.(fieldName);
        return;
    end
end

error('Required body part not found. Tried: %s. Available: %s', ...
    strjoin(candidateNames, ', '), strjoin(cellstr(available), ', '));
end

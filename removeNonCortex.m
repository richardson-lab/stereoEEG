function reducedSession = removeNonCortex(session, pathToLoc)
% Diego G. Davila
% Proekt Lab
% University of Pennsylvania School of Medicine
%
% This function removes white matter contacts and contacts not in brain
% tissue from an seeg session object. 
%
% INPUT(s):
% 1. session: an seeg session object
% 2. pathToLoc: path to the native and T1 RTK-SNAP output file containing
% contact labels
%
% OUTPUT(s):
% 1. reducedSession: session object with white matter and non-brain
% channels removed. The session object now also contains a field with a
% table matching contact labels with regions. It will also contain a vector
% describind the channel indexes that were removed. 
%
% Create a table that contains the anatomical location and contact name
locations = readtable(pathToLoc);
locations = locations(:,[1 2]);
locations.Properties.VariableNames = {'contact', 'region'};


% First, to match the formatting between session objects and location
% information given by ITK-SNAP output, we'll add a 0 before contact numbers in the location table.
for i = 1:length(locations.contact)
    if (length(locations.contact{i}) == 3)
        locations.contact{i} = [locations.contact{i}(:,1:2), '0', locations.contact{i}(:,3)];
    end
end

% create a copy to work with and return, and add the locations table
reducedSession = session;
reducedSession.locations = locations;

% loop through all the contacts and identify any that are either in white
% matter or not identified in brain
toRemove = zeros(size(session.data, 2), 1); % this will store the channels to be removed
for i = 1:length(session.channel_labels) 
    for j = 1:length(locations.contact) 
        if (strcmp(session.channel_labels{i}, locations.contact{j})) % if the IDs match
            if (contains(locations.region{j}, 'White Matter') | isempty(locations.region{j})) % if it has been labeled as white matter OR has no anatomical label
                toRemove(i) = 1; % set removal to 1
                %reducedSession.data(:,i) = []; % remove the non-cortex channel
            end
        end
    end
end
reducedSession.RemovedNonCortex = toRemove; % add the vector describing which channels were removed

% Delete the channels
reducedSession.data(:,find(toRemove)) = [];

end

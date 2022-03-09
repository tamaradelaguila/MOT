%% AUXILIAR/Y TEMPORAL CODE: SELECT VISUALLY REJECTED TRIALS


for nfish = 1:6
clearvars -except path nfish

%load data
VSDI = MOT('load', nfish);
VSDI.ref

% % sharks column to zeros 
% for ii = 1:length(VSDI.list)
%     VSDI.list(ii).shark = 0;
% end

% MANUALLY INSERT SHARKS AS '1'
% return

% aquí
sharks  = [VSDI.list(:).shark];

idx = find(sharks); 

VSDI.reject.visual = idx;

MOT('save', VSDI)

% MANUALLY DELETE COLUMN 
end 
% Created: 09/02/22
% Last Update: 
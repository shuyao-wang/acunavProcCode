function dispArfi(resname,parFile)

clear clc

if ispc
    addpath C:\users\vrk4\Documents\GitHub\SC2000\arfiProcCode\
    addpath(genpath('C:\users\vrk4\Documents\GitHub\SC2000\acunavProcCode\'))
else
    addpath /emfd/vrk4/GitHub/SC2000/arfiProcCode
    addpath(genpath('/emfd/vrk4/GitHub/SC2000/acunavProcCode'))
end

% Check inputs and set default parameters
if nargin<1
    resname = dir('res*.mat');
    if length(resname)==0, error('res files does not exist, run procArfi_acunav on binary data'); end    
    resname = resname(end).name;
end
% Pull out time stamp from filename
if ischar(resname)
    if ~exist(resname), error('res file does not exist, run procArfi_acunav on binary data'); end 
    [basePath,resname] = fileparts(resname);
elseif ~ischar(resname)
    temp = dir('res*.mat');
    if length(temp)==0, error('res files do not exist, run procArfi_acunav on binary data'); end    
    if resname == -1
        temp_name = temp(end).name;resname=temp_name;clear temp temp_name
        [basePath,resname] = fileparts(resname);
    else
        temp_name = temp(resname).name;resname=temp_name;clear temp temp_name
        [basePath,resname] = fileparts(resname);
    end
end

timeStamp = resname(5:end);
resname = fullfile(basePath, strcat(resname, '.mat'));
if nargin<2 || isempty(parFile)
    parFile = fullfile(basePath, sprintf('par_%s.mat',timeStamp));
    % Check to see if there is a time stamped parameters file, if it is in current directory, or one directory level higher
    if ~exist(parFile, 'file'), parFile = fullfile(pwd, 'parameters.mat');end
    if ~exist(parFile, 'file'), parFile = fullfile(pwd, '..', 'parameters.mat');end
end

% Load res and par files
res = load(resname);
par = load(parFile);

% Display baseName and setName
fprintf(1, 'basename/setName: %s/%s\n',par.baseName,par.setName)

%% Motion Filter
timeRange = [-inf -0.15 3.3 3.5];
order = 1;

tmask = false(size(res.t));
for i = 1:2:length(timeRange)
    tmask(res.t>timeRange(i) & res.t<timeRange(i+1)) = true;
end
tmask(par.nref+(1:par.npush*length(par.pushFocalDepth))) = false;
res.arfidata = linearmotionfilter(res.arfidata,res.t,find(tmask),order);


%% Cubic interpolate push and reverb
tidx1 = [par.nref+[-1 0] par.nref+par.npush+[2:3]];
tidx2 = [par.nref+[1:par.npush+1]];
[residtmp motion1] = linearmotionfilter(res.arfidata,res.t,tidx1,3);
res.arfidata(:,:,tidx2) = motion1(:,:,tidx2);
clear residtmp motion1 tidx1 tidx2;

%% Scan Convert
sc = struct(...
    'min_phi', res.angular(1),...
    'span_phi',res.angular(end)-res.angular(1),...
    'apex',0.1*res.apex,...
    'fs',par.fs*par.interpFactor*1e6...
    );


if strcmpi(par.mode,'ARFI')
[res.scarfidata,ax,lat] = scan_convert('sector',res.arfidata,sc.min_phi,sc.span_phi,sc.apex,2,sc.fs);
end

%% Display 
disp_lim = [-1 10];

figure

subplot(2,2,[1 3])
set(gca,'Position',[0.1 0.2 0.4 0.4])
for i=1:size(res.arfidata,3)
    if strcmpi(par.mode,'ARFI')
        imagesc(lat,ax,res.scarfidata(:,:,i),disp_lim)
        xlabel('Lateral (cm)')
        ylabel('Axial (cm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
        axis image
    elseif strcmpi(par.mode,'MMODE')
        imagesc([],res.radial,res.arfidata(:,:,i),disp_lim)
        xlabel('Push #')
        ylabel('Axial (cm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
    end
    colorbar
    pause
end
axgate = par.pushFocalDepth + [-1 1];
idxgate = [find(res.radial>axgate(1),1) find(res.radial>axgate(2),1)];
for i=1:size(res.arfidata,2)
    for j=1:size(res.arfidata,3)
        disp(i,j) = mean(res.arfidata(idxgate(1):idxgate(2),i,j));
    end
end

subplot(222)
if strcmpi(par.mode,'ARFI')
    imagesc(res.angular,res.t,disp',disp_lim)
    xlabel('Angle')
    ylabel('Track Time (ms)')
elseif strcmpi(par.mode,'MMODE')
    imagesc([],res.t,disp',disp_lim)
    xlabel('Push #')
    ylabel('Track Time (ms)')
end
subplot(224)
plot(res.t,disp')
xlabel('Track Time (ms)')
ylabel('Diplacement (\mum)')
title(sprintf('[%d-%d]\n nsamp = %d', axgate(1),axgate(2),idxgate(2)-idxgate(1)))
grid on
ylim(disp_lim)
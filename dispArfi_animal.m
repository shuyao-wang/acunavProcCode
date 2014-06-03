function dispArfi(resname,parFile)

clear clc
close all

if ispc
    addpath C:\users\vrk4\Documents\GitHub\SC2000\arfiProcCode\
    addpath(genpath('C:\users\vrk4\Documents\GitHub\acunavProcCode\'))
else
    addpath /emfd/vrk4/GitHub/SC2000/arfiProcCode
    addpath(genpath('/emfd/vrk4/GitHub/acunavProcCode'))
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
if par.ref_idx == -1
    tidx1 = [par.nref+[-1 0] par.nref+par.npush+[3:4]];
    tidx2 = [par.nref+[1:par.npush+2]];
else
    tidx1 = [par.nref+[-1 0] par.nref+par.npush+[2:3]];
    tidx2 = [par.nref+[1:par.npush+1]];
end
[residtmp motion1] = linearmotionfilter(res.arfidata,res.t,tidx1,3);
res.arfidata(:,:,tidx2) = motion1(:,:,tidx2);
clear residtmp motion1 tidx1 tidx2;

%% Get Bmode

res.bimg = abs(complex(res.I,res.Q));
res.bimg = db(res.bimg/max(res.bimg(:)));


%% Scan Convert
sc = struct(...
    'min_phi', res.angular(1),...
    'span_phi',res.angular(end)-res.angular(1),...
    'apex',0.1*res.apex,...
    'fs',par.fs*par.interpFactor*1e6...
    );


if strcmpi(par.mode,'ARFI')
    [res.scarfidata,axa,lata] = scan_convert('sector',res.arfidata,sc.min_phi,sc.span_phi,sc.apex,2,sc.fs);
    [res.scbimg,axb,latb] = scan_convert('sector',res.bimg,sc.min_phi,sc.span_phi,sc.apex,2,sc.fs);
end

%% Display
disp_lim = [-1 8];
blim = [-60 0];

figure(1)
set(1,'Position',[1950 350 560 420])
figure(2)
set(gcf,'Position',[2530 350 560 420])

if strcmpi(par.mode,'ARFI')
    figure(1)
    imagesc(latb,axb,res.scbimg(:,:,1),blim)
    xlabel('Lateral (cm)')
    ylabel('Axial (cm)')
    title(sprintf('probe: %s\n Bmode \n',par.probe),'interpreter','none')
    colormap(gray)
    axis image
    colorbar
elseif strcmpi(par.mode,'MMODE')
    tt = 1000*[0:1/par.pushPRF:(par.numBeamGroups-1)/par.pushPRF];
    figure(1)
    imagesc(tt,res.radial,res.bimg(:,:,1),blim)
    xlabel('Time (ms)')
    ylabel('Axial (cm)')
    title(sprintf('probe: %s\n Bmode \n',par.probe),'interpreter','none')
    colormap(gray)
    colorbar
end

for i=1:size(res.arfidata,3)
    if strcmpi(par.mode,'ARFI')
        figure(2)
        imagesc(lata,axa,res.scarfidata(:,:,i),disp_lim)
        xlabel('Lateral (cm)')
        ylabel('Axial (cm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
        axis image
        colorbar
    elseif strcmpi(par.mode,'MMODE')
        figure(2)
        imagesc(tt,res.radial,res.arfidata(:,:,i),disp_lim)
        xlabel('Time (ms)')
        ylabel('Axial (cm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
    end
    pause
end

axgate = par.pushFocalDepth + [-1 1];
idxgate = [find(res.radial>axgate(1),1) find(res.radial>axgate(2),1)];
for i=1:size(res.arfidata,2)
    for j=1:size(res.arfidata,3)
        disp(i,j) = mean(res.arfidata(idxgate(1):idxgate(2),i,j));
    end
end

figure(3)
set(gcf,'Position',[1950 17 1140 242])
subplot(121)
if strcmpi(par.mode,'ARFI')
    imagesc(res.angular,res.t,disp',disp_lim)
    xlabel('Angle')
    ylabel('Track Time (ms)')
    colorbar
elseif strcmpi(par.mode,'MMODE')
    imagesc(tt,res.t,disp',disp_lim)
    xlabel('Time (ms)')
    ylabel('Track Time (ms)')
    colorbar
end
subplot(122)
for i=1:par.numBeamGroups
    cla
    plot(res.t,disp(1+(i-1)*par.nBeams:i*par.nBeams,:)','.-')
    xlabel('Track Time (ms)')
    ylabel('Diplacement (\mum)')
    if strcmpi(par.mode,'ARFI')
        title(sprintf('[%d-%d]\n nsamp = %d\n Angle = %d', axgate(1),axgate(2),idxgate(2)-idxgate(1),round(res.angular(2+(i-1)*par.nBeams))))
    elseif strcmpi(par.mode,'MMODE')
        title(sprintf('[%d-%d]\n nsamp = %d\n Time = %d', axgate(1),axgate(2),idxgate(2)-idxgate(1),tt(i)))
    end
    grid on
    ylim(disp_lim)
    pause
end

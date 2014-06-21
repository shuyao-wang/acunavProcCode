function dispArfi(resname,parFile)

clear clc
close all

t = 0.5;
timeRange = [-inf -0.15 3.3 3.4];
disp_lim = [-1 8];
blim = [-40 0];
ax_lim = [1 30];


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
if strcmpi(par.probeType,'phased')
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
end

%% Display


figure(1)
if ispc
    set(1,'Position',[5 400 450 375])
else
    set(1,'Position',[1950 350 560 420])
end
figure(2)
if ispc
    set(gcf,'Position',[475 400 450 375])
else
    set(gcf,'Position',[2530 350 560 420])
end

if (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'phased'))
    figure(1)
    imagesc(10*latb,10*axb,res.scbimg(:,:,1),blim)
    xlabel('Lateral (mm)')
    ylabel('Axial (mm)')
    title(sprintf('probe: %s\n Bmode \n',par.probe),'interpreter','none')
    colormap(gray)
    axis image
    ylim(ax_lim)
    colorbar
elseif (strcmpi(par.mode,'MMODE') && strcmpi(par.probeType,'phased'))
    tt = 1000*[0:1/par.pushPRF:(par.numBeamGroups-1)/par.pushPRF];
    figure(1)
    imagesc(tt,res.radial,res.bimg(:,:,1),blim)
    xlabel('Time (ms)')
    ylabel('Axial (mm)')
    title(sprintf('probe: %s\n Bmode \n',par.probe),'interpreter','none')
    colormap(gray)
    ylim(ax_lim)
    colorbar
elseif (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'linear'))
    figure(1)
    imagesc(res.lat,res.axial,res.bimg(:,:,1),blim)
    xlabel('Time (ms)')
    ylabel('Axial (mm)')
    title(sprintf('probe: %s\n Bmode \n',par.probe),'interpreter','none')
    axis image
    colormap(gray)
    ylim(ax_lim)
    colorbar
end

figure(2)
colormap(gray)
for i=1:size(res.arfidata,3)
    if (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'phased'))
        figure(2)
        imagesc(10*lata,10*axa,res.scarfidata(:,:,i),disp_lim)
        xlabel('Lateral (mm)')
        ylabel('Axial (mm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
        axis image
        ylim(ax_lim)
        colorbar
    elseif (strcmpi(par.mode,'MMODE')&& strcmpi(par.probeType,'phased'))
        figure(2)
        imagesc(tt,res.radial,res.arfidata(:,:,i),disp_lim)
        xlabel('Time (ms)')
        ylabel('Axial (mm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
        ylim(ax_lim)
    elseif (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'linear'))
        figure(2)
        imagesc(res.lat,res.axial,res.arfidata(:,:,i),disp_lim)
        xlabel('Lat (mm)')
        ylabel('Axial (mm)')
        title(sprintf('probe: %s\n baseName/setName: %s/%s\n',par.probe,par.baseName,par.setName),'interpreter','none')
        axis image
        ylim(ax_lim)
    end
    pause(0.03)
end

idx = find(res.t>t,1);
if (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'phased'))
    figure(2)
    imagesc(10*lata,10*axa,res.scarfidata(:,:,idx),disp_lim)
    xlabel('Lateral (mm)')
    ylabel('Axial (mm)')
    title(sprintf('probe: %s\n baseName/setName: %s/%s\n t=%1.1f ms',par.probe,par.baseName,par.setName,res.t(idx)),'interpreter','none')
    axis image
    ylim(ax_lim)
    colorbar
elseif (strcmpi(par.mode,'MMODE')&& strcmpi(par.probeType,'phased'))
    figure(2)
    imagesc(tt,res.radial,res.arfidata(:,:,idx),disp_lim)
    xlabel('Time (ms)')
    ylabel('Axial (mm)')
    title(sprintf('probe: %s\n baseName/setName: %s/%s\n t=%1.1f ms',par.probe,par.baseName,par.setName,res.t(idx)),'interpreter','none')
    ylim(ax_lim)
elseif (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'linear'))
    figure(2)
    imagesc(res.lat,res.axial,res.arfidata(:,:,idx),disp_lim)
    xlabel('Lat (mm)')
    ylabel('Axial (mm)')
    title(sprintf('probe: %s\n baseName/setName: %s/%s\n t=%1.1f ms',par.probe,par.baseName,par.setName,res.t(idx)),'interpreter','none')
    axis image
    ylim(ax_lim)    
end

if isfield(res,'axial'); res.radial = res.axial;end

axgate = par.pushFocalDepth +2 + [-1 1];
idxgate = [find(res.radial>axgate(1),1) find(res.radial>axgate(2),1)];
for i=1:size(res.arfidata,2)
    for j=1:size(res.arfidata,3)
        disp(i,j) = mean(res.arfidata(idxgate(1):idxgate(2),i,j));
    end
end

figure(3)
if ispc
    set(gcf,'Position',[5 50 920 250])
else
    set(gcf,'Position',[1950 17 1140 242])
end
subplot(121)
if (strcmpi(par.mode,'ARFI') && strcmpi(par.probeType,'phased'))
    imagesc((abs(res.apex)+par.pushFocalDepth)*tand(res.angular),res.t,disp',disp_lim)
    xlabel('Lat (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('Displacement at %d mm',par.pushFocalDepth))
    colorbar
elseif (strcmpi(par.mode,'MMODE')&& strcmpi(par.probeType,'phased'))
    imagesc(tt,res.t,disp',disp_lim)
    xlabel('Time (ms)')
    ylabel('Track Time (ms)')
    title(sprintf('Displacement at %d mm',par.pushFocalDepth))
    colorbar
elseif (strcmpi(par.mode,'ARFI')&& strcmpi(par.probeType,'linear'))
    imagesc(res.lat,res.t,disp',disp_lim)
    xlabel('Lat (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('Displacement at %d mm',par.pushFocalDepth))
    colorbar
end
subplot(122)
for i=1:par.numBeamGroups
    hold on
    plot(res.t,disp(1+(i-1)*par.nBeams:i*par.nBeams,:)','b.-')
    xlabel('Track Time (ms)')
    ylabel('Diplacement (\mum)')
    if strcmpi(par.mode,'ARFI')
        title(sprintf('[%d-%d]\n nsamp = %d\n', axgate(1),axgate(2),idxgate(2)-idxgate(1)))
    elseif strcmpi(par.mode,'MMODE')
        title(sprintf('[%d-%d]\n nsamp = %d\n', axgate(1),axgate(2),idxgate(2)-idxgate(1)))
    end
    grid on
    ylim(disp_lim)
    pause(0.03)
end

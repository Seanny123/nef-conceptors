
set(0,'DefaultFigureWindowStyle','normal');

d = mcread('13_04_sitStandup.c3d');

if 0
[mf, mm, mgrid] = mcmissing(d); 
figure(1);
set(gcf,'Position',[40 200 560 420]); 
subplot(3,1,1);
bar(mf); xlabel('Marker'); ylabel('Num. of Missing frames');
subplot(3,1,2); 
bar(mm); xlabel('Frame'); ylabel('Num. of Missing markers');
subplot(3,1,3);
imagesc(-mgrid');
colormap gray; xlabel('Frame'); ylabel('Marker');
end

%% delete markers which are almost never visible

keepmarkers = 1:41;
d.markerName = d.markerName(keepmarkers,1);
nEffectiveFrames = size(d.data,1);
newData = zeros(nEffectiveFrames,0);
for i = 1:length(keepmarkers)
    newData = [newData, ...
        d.data(:, 3*(keepmarkers(i)-1)+1:3*keepmarkers(i))];
end
d.data = newData;
d.nMarkers = length(keepmarkers);
d.nFrames = size(d.data,1);
%%
if 0
mcplottimeseries(d, [1 2 5 8 ], 'dim', 3, 'timetype', 'frame');
end 

%% prune time
keeptimes = 600:800;
d.data = d.data(keeptimes,:);
d.nFrames = length(keeptimes);



%%

d = mcsmoothen(d);
if 0
mcplottimeseries(d, [1 2 5 7 ], 'dim', 3, 'timetype', 'frame');
end 

%% transform to Dempster model

% a helper function to get marker names printed
if 0
    for i = 1:d.nMarkers
        disp(sprintf('%g %s',i,d.markerName{i,1}));
    end
end
%%
m2jpar.type = 'm2jpar';
m2jpar.nMarkers = 20;

m2jpar.markerNum{1,12} = [30 29 33 34]; % head LBHD RBHD LFHD RFHD
m2jpar.markerNum{1,11} = [23 27]; % neck CLAV C7
m2jpar.markerNum{1,13} = [26]; % right shoulder RSHO
m2jpar.markerNum{1,17} = [25]; % left shoulder LSHO
m2jpar.markerNum{1,14} = [16]; % right elbow RELB
m2jpar.markerNum{1,18} = [14]; % left elbow LELB
m2jpar.markerNum{1,15} = [15 7]; % right wrist RWRB RWRA
m2jpar.markerNum{1,19} = [9 12]; % left wrist LWRA LWRB
m2jpar.markerNum{1,16} = [11]; % right hand RFIN
m2jpar.markerNum{1,20} = [13]; % left hand LFIN
m2jpar.markerNum{1,10} = [19 20]; % back center T10 STRN
m2jpar.markerNum{1,1} = [4 3 1 2 ]; % root (back bottom) 
                                       % RFWT LFWT RBWT LBWT
m2jpar.markerNum{1,2} = [4 1]; % right hip RFWT RBWT
m2jpar.markerNum{1,6} = [3 2]; % left hip LFWT LBWT
m2jpar.markerNum{1,3} = [18]; % right knee RKNE
m2jpar.markerNum{1,7} = [17]; % left knee LKNE
m2jpar.markerNum{1,4} = [36 32]; % right ankle  RHEE RANK
m2jpar.markerNum{1,8} = [35 31]; % left ankle LHEE LANK
m2jpar.markerNum{1,5} = [40 39]; % right foot RMT5 RTOE
m2jpar.markerNum{1,9} = [37 38]; % left foot LMT5 LTOE

m2jpar.markerName{1,12} = 'head';
m2jpar.markerName{1,11} = 'neck';
m2jpar.markerName{1,13} = 'rsho';
m2jpar.markerName{1,17} = 'lsho';
m2jpar.markerName{1,14} = 'relb';
m2jpar.markerName{1,18} = 'lelb';
m2jpar.markerName{1,15} = 'rwri';
m2jpar.markerName{1,19} = 'lwri';
m2jpar.markerName{1,16} = 'rhan';
m2jpar.markerName{1,20} = 'lhan';
m2jpar.markerName{1,10} = 'back';
m2jpar.markerName{1,1} = 'root';
m2jpar.markerName{1,2} = 'rhip';
m2jpar.markerName{1,6} = 'lhip';
m2jpar.markerName{1,3} = 'rkne';
m2jpar.markerName{1,7} = 'lkne';
m2jpar.markerName{1,4} = 'rank';
m2jpar.markerName{1,8} = 'lank';
m2jpar.markerName{1,5} = 'rfoo';
m2jpar.markerName{1,9} = 'lfoo';

dj = mcm2j(d, m2jpar);



%% to segment model
j2spar.type = 'j2spar';
j2spar.rootMarker = 1;
j2spar.frontalPlane = [6 10 2 ];  % be careful about order
j2spar.parent = [0 1 2 3 4 1 6 7 8 1 10 11 11 13 14 15 11 17 18 19];
j2spar.segmentName = cell(1,19);




%% create NN raw data

[nnRawDataSitting, hmeanSitting, segLengthsSitting] = ...
    j2nnRawHJ_absH(dj,j2spar);

%% backtransform to joint data for test and show video
if 0
freq = 120; 
djBack =  nnRaw2jHJ_absH(nnRawDataSitting, hmeanSitting, ...
     segLengthsSitting, j2spar, freq);
 
 
    japarVideo.showmnum = 0;
    japarVideo.animate = 1;
    japarVideo.colors = 'wkkkk';
    japarVideo.trl = 0;
    japarVideo.numbers = [];
    japarVideo.cwidth = 5;
    japarVideo.twidth = 1;
    japarVideo.az = 0;
    japarVideo.el = 20;
    japarVideo.fps = 30;
    japarVideo.limits = [];
    japarVideo.scrsize = [300 350];
    %japarVideo.scrsize = [400 300];
    japarVideo.showfnum = 0;
    japarVideo.conn2 = [];
    japarVideo.conn = [11 12; 13 11; 11 17; 13 14; 14 15; 15 16;...
        17 18; 18 19; 19 20; 10 11; 1 10; 1 2; 1 6; 2 3; 3 4; 4 5;...
        6 7; 7 8; 8 9];
    japarVideo.msize = 8;
    japarVideo.fontsize = 12;
    japarVideo.folder = '0'; % set to string '0' if no save is wanted
    japarVideo.center = 1;
    japarVideo.pers.c = [0 -4000 0];
    japarVideo.pers.th=[0 0 0];
    japarVideo.pers.e=[0 -2000 0];
     japarVideo.nGridlines = 5;
     japarVideo.botWidth = 2000;    
    japarVideo.fignr = 11;
    
    djBackResampled = mcresampleHJ(djBack, japarVideo.fps);
    xRootShifts = djBackResampled.data(2:end,1) - ...
        djBackResampled.data(1:end-1,1);
    yRootShifts = djBackResampled.data(2:end,2) - ...
        djBackResampled.data(1:end-1,2);
    
    djrecoveredCentered = mccenterxyFrames(djBack);
    djrecoveredCentered = mcresampleHJ(djrecoveredCentered, japarVideo.fps);
    
mcanimateHJ(djrecoveredCentered, japarVideo, 1, xRootShifts, yRootShifts);
end 


%%

if 1
    save nnRawSitting nnRawDataSitting hmeanSitting ...
         segLengthsSitting;
end




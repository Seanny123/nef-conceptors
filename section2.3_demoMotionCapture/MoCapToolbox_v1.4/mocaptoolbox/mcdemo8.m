function mcdemo8
% This example shows how you can color your plots and animations
% with the MoCap toolbox. 
%
% You can only produce the animation frames with the
% toolbox. These have to be compiled into a movie
% using some other software, such as QuickTime Pro
% (for Mac) or ??? (for Windows)

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
clc
% Let us create a colored animation from the variable dance2.
% The animation parameters look like this:

load mcdemodata
mapar

% Let us change the frames-per-second value to 15 
mapar.fps = 15;

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

% Let us set individual colors for six markers (head front left, head back right, shoulder left,
% hip left back, finger right, knee left, knee right, heel left)
mapar.markercolors='bwwgwrwwwwwywwwwwwmwcbwg';

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

% and let us have a look at the new colors:
mcplotframe(dance2, 150, mapar);

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
clc
% Looks good. Now let us set the markers that we want to trace and the trace length (in seconds)
mapar.trm=[1 6 12 19 21 24];
mapar.trl=3;

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

% Let us set individual colors for the traces
mapar.tracecolors='grymcb';

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

% We rotate the figure to be frontal on average

dance2=mc2frontal(dance2,9,10);

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
clc
% In the next dialog window, choose the directory 
% where you want the frames to be stored:

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
path = uigetdir([], 'Pick a Directory')

olddir = cd; cd(path)

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we make the animation:

newpar = mcanimate(dance2, mapar);

% The animation frames are created and saved in the folder

disp([path '/' mapar.folder])

% The next thing to do is to launch a video editing program (such
% as QuickTime Pro) and read the frames into a video clip


pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

cd(olddir)

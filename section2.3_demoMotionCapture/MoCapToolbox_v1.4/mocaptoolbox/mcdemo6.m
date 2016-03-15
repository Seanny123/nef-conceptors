function mcdemo6
% This example shows how you can create animations
% with the MoCap toolbox. 
%
% You can only produce the animation frames with the
% toolbox. These have to be compiled into a movie
% using some other software, such as QuickTime Pro
% (for Mac) or ??? (for Windows)

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
clc
% Let us create an animation from the variable walk2.
% The animpar structure mapar contains the connector
% information for this variable:

load mcdemodata
mapar

% Let us change the frames-per-second value to 15 
mapar.fps = 15;

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
clc
% In the next dialog window, choose the directory 
% where you want the frames to be stored:

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
path = uigetdir([], 'Pick a Directory')

olddir = cd; cd(path)

pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the animation:

newpar = mcanimate(walk2, mapar);

% The animation frames are created and saved in the folder

disp([path '/' mapar.folder])

% The next thing to do is to launch a video editing program (such
% as QuickTime Pro) and read the frames into a video clip


pause %%%%%%%%%%% hit a key to continue %%%%%%%%%%%%%%%%%%%%%%%%%%

cd(olddir)

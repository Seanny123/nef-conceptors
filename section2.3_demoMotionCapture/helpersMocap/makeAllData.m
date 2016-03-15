addpath('./rawMocapData');

disp('generating data for ExaggeratedStride');
makeNNData_ExaggeratedStride;
disp('generating data for SlowWalk');
makeNNData_SlowWalk;
disp('generating data for Walk');
makeNNData_Walk;
disp('generating data for Jog');
makeNNData_RunJog;
disp('generating data for Cartwheel');
makeNNData_Cartwheel;
disp('generating data for Waltz');
makeNNData_Waltz;
disp('generating data for Crawl');
makeNNData_Crawl;
disp('generating data for Standup');
makeNNData_Standup;
disp('generating data for Getdown');
makeNNData_Getdown;
disp('generating data for sitting');
makeNNData_sitting;
disp('generating data for standupFromStool');
makeNNData_standupFromStool;
disp('generating data for SitDownOnStool');
makeNNData_SitDownOnStool;
disp('generating data for Box1');
makeNNData_Box1;
disp('generating data for Box2');
makeNNData_Box2;
disp('generating data for Box3');
makeNNData_Box3;

% code for studies entitled:
% Extended SEIR model for death and cure population of COVID-19 in China 
% under public intervention and clinical progress
% Su Feng, Yuan Peijiang, Li Jianmin
% modified in 2020/03/29
%% 
clc
clear
load('TimeInd.mat')  % Modeling start from 2020/1/21

global ParaSet
ParaSet=[];

N_Days = 120; % days of modeling
% para initialization
ParaSet.Init = [1*1e9 0 0 0 0]; % Initialized N S E I R
ParaSet.model=[0 0 0 0];        % transmission rate by I and E; infection rate; recovery rate
ParaSet.Death=[0 0];            % parameter of death | death lag death rate
ParaSet.T = N_Days; 
ParaSet.Control = [0 0; ...     % phase-adjusted public intervention
                   0 0];        

ParaSet.Init(3)=500;            % *********** para
ParaSet.Init(4)=150; 
ParaSet.Init(5)=0;
ParaSet.Init(2)=ParaSet.Init(1)-sum(ParaSet.Init(3:5));
ParaSet.Init_1=[0 0];

ParaSet.model(1:2)=[0.09, 0.57];        % *********** para
ParaSet.model(3:4)=[1/7, 1/35];         % *********** para (half)
ParaSet.Control = [6 0.10; 16 0.05];    % *********** para

ParaSet.DeathMedical = [11 0.95 1.05];  % clinical progress para(8:9)
ParaSet.DeathMedicalMature = 62;        % mature of clinical progress
ParaSet.Death = [3 0.001];  % *********** para (10:11)
ParaSet.Cured = 5;          % *********** cure lag  para (12)
ParaSet.Suspect=[2 0.01];   
ModRes = SEIR_DC_COVID_19_V1(ParaSet);  % SEIR-DC model
TR_ind = 1:ParaSet.T;
TimeSeq = TimeInd(TR_ind,2);
ModRes.TimeSeq=TimeSeq;                 % **** time index

ErrOutput=Plot_ModRes_Plot_V1(ModRes,ParaSet); % visualization




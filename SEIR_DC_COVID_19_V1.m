function ModRes = SEIR_DC_COVID_19_V1(ParaSet)
% code for studies entitled:
% "Extended SEIR model for death and cure population of COVID-19 in China 
% under public intervention and clinical progress"
% Su Feng, Yuan Peijiang, Li Jianmin
% modified in 2020/03/29

% % ----------- 
rB=ParaSet.model(1);             % transmission rate by I
r2B2=ParaSet.model(2);           % transmission rate by  E
a=ParaSet.model(3);              % infection rate
y=ParaSet.model(4);              % recovery rate

ModRes=[];
ModRes.Label={'S','E','I','R'};
ModRes.Data=zeros(ParaSet.T,4);  % length*4(SEIR)
N=ParaSet.Init(1);
ModRes.Data(1,:)=ParaSet.Init(2:5);

S = ModRes.Data(:,1);
E = ModRes.Data(:,2);
I = ModRes.Data(:,3);
R = ModRes.Data(:,4);

rB_T=rB*ones(ParaSet.T,1);  % phase-adjusted public intervention
r2B2_T=r2B2*ones(ParaSet.T,1);
for ci=1:size(ParaSet.Control,1)
    rB_T(ParaSet.Control(ci,1):end)=rB_T(ParaSet.Control(ci,1):end).*ParaSet.Control(ci,2);
    r2B2_T(ParaSet.Control(ci,1):end)=r2B2_T(ParaSet.Control(ci,1):end).*ParaSet.Control(ci,2);
end

TRRatio=ones(ParaSet.T,1);
TR_N=ParaSet.DeathMedical(1:2);
TRRatio(TR_N(1):end)=ParaSet.DeathMedical(3);
xx = (1:ParaSet.T)';
TRRatio=TRRatio.^xx;
if length(TRRatio)>= ParaSet.DeathMedicalMature
    TRRatio(ParaSet.DeathMedicalMature:end) = TRRatio(ParaSet.DeathMedicalMature-1);
end

y_T=y.*TRRatio;
DeltaI=0*I;
ConfirmInject=zeros(ParaSet.T,1);
if isfield(ParaSet,'Injection')
    ConfirmInject(ParaSet.Injection(:,1)-1)=ParaSet.Injection(:,2);
end
for ci = 1:ParaSet.T-1
    rB=rB_T(ci);
    r2B2=r2B2_T(ci);
    y=y_T(ci);  % clinical progress
    S(ci+1) = S(ci) - rB*S(ci)*I(ci)/N(1) - r2B2*S(ci)*E(ci)/N;
    E(ci+1) = E(ci) + rB*S(ci)*I(ci)/N(1)-a*E(ci) + r2B2*S(ci)*E(ci)/N;
    I(ci+1) = I(ci) + a*E(ci) - y*I(ci) + ConfirmInject(ci);
    R(ci+1) = R(ci) + y*I(ci);
    
    DeltaI(ci+1) = a*E(ci) + ConfirmInject(ci);
end
ModRes.Data = [S E I R];
ModRes.DeltaI = DeltaI;
ModRes.CumsumI=cumsum(DeltaI);

% %% death
TRRatio=ones(ParaSet.T,1); % *ParaSet.Death(2)
TR_N=ParaSet.DeathMedical(1:2);
TRRatio(TR_N(1):end)=TRRatio(TR_N(1):end)*ParaSet.DeathMedical(2);

xx = (1:ParaSet.T)';
TRRatio=TRRatio.^xx;
if length(TRRatio)>= ParaSet.DeathMedicalMature
    TRRatio(ParaSet.DeathMedicalMature:end) = TRRatio(ParaSet.DeathMedicalMature-1);
end

TR2=ModRes.Data(:,3).*ParaSet.Death(2).*TRRatio;
TR2=cumsum(TR2);           % convert newly death to accumulative death
TR21=[nan(ParaSet.Death(1)-1,1); TR2];
TR21(end-ParaSet.Death(1)+2:end)=[];
ModRes.Death=TR21;

% %% cure
TR1=ModRes.Data(:,4)-TR21;  % death + cure = R
TR2=zeros(ParaSet.T,1);
TR2(ParaSet.Cured+1:ParaSet.T)=TR1(1:ParaSet.T-ParaSet.Cured);
ModRes.Cure=TR2;

end


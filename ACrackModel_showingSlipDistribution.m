%% calculate C for a crack model
clear all
close all
clc

%% Set up and do one time test
ct = 3e3; %m/s; S
cl = sqrt(3)*ct; %m/s; P
pho = 3e3 ; %density, kg/m^3

mu = pho*ct^2;
lambda= pho*cl^2-2*mu;

nu = (cl^2-2*ct^2)/(2*(cl^2-ct^2)); % Poisson ratio

BurialDepth = 0;

dip = 15; %deg
fault_strike_width = 2e3; %m
fault_dip_width = 1e3; %m
CenterPoint_Depth = fault_dip_width/2*sind(dip)+BurialDepth*fault_dip_width; %m
N_strike = 100;
N_dip = 50;

% dip thrust tress drop setup
DipStressChange = -1e6; %1MPa
StrikeStressChange = 0; %

% strike stress drop setup
% DipStressChange = 0; 
% StrikeStressChange = 1e6; %1MPa


[SubFault, GlobalParameters] = DiscFault(lambda,mu,CenterPoint_Depth,dip,fault_strike_width,fault_dip_width,N_strike,N_dip);

%% Calculate stress kernel dip-dip
%  DISL1 : Dislocation in strike-slip component
%  DISL2 : Dislocation in dip-slip component
%  DISL3 : Dislocation in tensile component
%  If 0<DIP<90, positive DISL1 gives left-lateral movement
%  and positive DISL2 gives reverse-fault movement
%  while positive DISL3 gives open-type movement
dislocation_dip = [0 1 0]; % reverse
dislocation_strike = [1 0 0]; % left-lateral

% direction vector
alongstrike_vec=[1 0 0];
alongdip_vec=[0 cosd(dip) sind(dip)]; % z is positive, dip~0-90;
normal_vec=cross(alongstrike_vec,alongdip_vec);
normal_vec=normal_vec'./norm(normal_vec,2); %transverse to make the normal unit vector 3*1, a colomn vector
alongstrike_vec = alongstrike_vec';%transverse to make the alongstrike unit vector 3*1, a colomn vector
alongdip_vec = alongdip_vec';%transverse to make the alongdip unit vector 3*1, a colomn vector

%test GetTractionKernal with a single value
%observation_point = [0 0 -2000];
%[ShearStrike_traction,ShearDip_traction,Normal_traction] = GetTractionKernal(lambda, mu,observation_point ,CenterPoint_Depth, dip, strike_width, dip_width, dislocation_dip);

%% Calculate K_DipDip (dip slip causing dip stress) and K_DipStrike (dip slip causing strike stress)
% K_DipNormal (although not used at the moment, still calculated, may use in the future, 2020.12.03)
%Allocate Kernal Array
K_DipDip = zeros(length(SubFault));
K_DipStrike = zeros(length(SubFault));
K_DipNormal = zeros(length(SubFault));
%Global parameter
CenterPoint_Depth = GlobalParameters.CenterPoint_Depth;
dip = GlobalParameters.dip_deg;
mu = GlobalParameters.mu;
lambda = GlobalParameters.lambda;
dislocation = dislocation_dip;

Total_time = 0;
ElaspTime = 0;
tic
for ind_source = 1:length(SubFault)
    Total_time = Total_time + ElaspTime;
    if mod(ind_source,floor(length(SubFault)/20))==0
        disp(['K_Dip Current progress: ',num2str(ind_source/length(SubFault)*100,'%.1f'),...
            '%, ElaspTime = ',num2str(ElaspTime,'%.1f'),...
            's, TotalTime = ',num2str(Total_time,'%.1f'),'s\n'])
        tic
    end
    for ind_obs = 1:length(SubFault)
        
        % obs point info
        observation_point_now = [SubFault(ind_obs).ThisSegCenterPoint_x SubFault(ind_obs).ThisSegCenterPoint_y -SubFault(ind_obs).ThisSegCenterPoint_Depth];
        % source point info
        strike_width_now = SubFault(ind_source).strike_width;
        dip_width_now = SubFault(ind_source).dip_width;
        
        [ShearStrike_traction,ShearDip_traction,Normal_traction] = GetTractionKernal(lambda, mu,observation_point_now ,CenterPoint_Depth, dip, strike_width_now, dip_width_now, dislocation);
        K_DipStrike(ind_source,ind_obs) = ShearStrike_traction' * alongstrike_vec;
        K_DipDip(ind_source,ind_obs) = ShearDip_traction' * alongdip_vec;
        K_DipNormal(ind_source,ind_obs) = Normal_traction' * normal_vec;
        
    end  
    ElaspTime = toc;
end

%% Calculate K_StrikeDip (Strike slip causing dip stress) and K_StrikeStrike (strike slip causing strike stress)
% K_StrikeNormal (although not used at the moment, still calculated, may use in the future, 2020.12.03)
%Allocate Kernal Array
K_StrikeDip = zeros(length(SubFault));
K_StrikeStrike = zeros(length(SubFault));
K_StrikeNormal = zeros(length(SubFault));
%Global parameter
CenterPoint_Depth = GlobalParameters.CenterPoint_Depth;
dip = GlobalParameters.dip_deg;
mu = GlobalParameters.mu;
lambda = GlobalParameters.lambda;
dislocation = dislocation_strike;

Total_time = 0;
ElaspTime = 0;
tic
for ind_source = 1:length(SubFault)
    Total_time = Total_time + ElaspTime;
    if mod(ind_source,floor(length(SubFault)/20))==0
        disp(['K_Strike Current progress: ',num2str(ind_source/length(SubFault)*100,'%.1f'),'%, ElaspTime = ',num2str(ElaspTime,'%.1f'),'s, TotalTime = ',num2str(Total_time,'%.1f'),'s\n'])
    end
    tic
    for ind_obs = 1:length(SubFault)
        
        % obs point info
        observation_point_now = [SubFault(ind_obs).ThisSegCenterPoint_x SubFault(ind_obs).ThisSegCenterPoint_y -SubFault(ind_obs).ThisSegCenterPoint_Depth];
        % source point info
        strike_width_now = SubFault(ind_source).strike_width;
        dip_width_now = SubFault(ind_source).dip_width;
        
        [ShearStrike_traction,ShearDip_traction,Normal_traction] = GetTractionKernal(lambda, mu,observation_point_now ,CenterPoint_Depth, dip, strike_width_now, dip_width_now, dislocation);
        K_StrikeStrike(ind_source,ind_obs) = ShearStrike_traction' * alongstrike_vec;
        K_StrikeDip(ind_source,ind_obs) = ShearDip_traction' * alongdip_vec;
        K_StrikeNormal(ind_source,ind_obs) = Normal_traction' * normal_vec;
        
    end  
    ElaspTime = toc;
end

    
%% Inverse to obtain a crack model, dip stress drop ---> slip
Kernal_stress = cat(1,cat(2,K_DipDip,K_StrikeDip),cat(2,K_StrikeDip,K_StrikeStrike));

StressVector = cat(1,DipStressChange.*ones(length(SubFault),1),StrikeStressChange.*ones(length(SubFault),1));
SlipVector = Kernal_stress\StressVector;
DipSlip = SlipVector(1:length(SubFault));
StrikeSlip = SlipVector(length(SubFault)+1:end);

%using absolute amplitude
Slip_Amp = sqrt(StrikeSlip.^2+DipSlip.^2);
%using average's amplitude
Average_Slip_avg = sqrt(mean(StrikeSlip).^2+mean(DipSlip).^2);

% compute normal stress change
Kernel_NormalTractionChange = cat(2,K_DipNormal,K_StrikeNormal);
NormalTractionChange = Kernel_NormalTractionChange*SlipVector; % Positive, decrease normal stress

% assign along fault coordinate
AlongStrike_coor = zeros(length(SubFault),1);
AlongDip_coor = zeros(length(SubFault),1);
for ind=1:length(SubFault)
    AlongStrike_coor(ind) = SubFault(ind).ThisSegCenterPoint_strike;
    AlongDip_coor(ind) = SubFault(ind).ThisSegCenterPoint_dip;
end 


%% Calculate Effective k and h
Average_Slip = mean(Slip_Amp);
% using half fault width for characteristic length
C_k = max(abs(DipStressChange),abs(StrikeStressChange))/mu*(fault_dip_width)/Average_Slip; %k
C_k_avg = max(abs(DipStressChange),abs(StrikeStressChange))/mu*(fault_dip_width)/Average_Slip_avg; %k

Average_NormalTractionChange = mean(NormalTractionChange);
C_h = Average_NormalTractionChange/mu*(fault_dip_width)/Average_Slip; %h


%% plot slip
close all

% choose the distribution that was calculated before

%load('SlipDtrb_Ns100_Nd50_dip90_BrD0_StrikeSlip.mat')
%load('SlipDtrb_Ns100_Nd50_dip90_BrD0.02_StrikeSlip.mat')
%load('SlipDtrb_Ns100_Nd50_dip90_BrD4.5_StrikeSlip.mat')
%load('SlipDtrb_Ns100_Nd50_dip90_BrD100_StrikeSlip.mat')

%load('SlipDtrb_Ns100_Nd50_dip10_BrD0_DipSlip.mat');
%load('SlipDtrb_Ns100_Nd30_dip35_BrD0_DipSlip.mat');

dir = '/Users/baoning/Dropbox/2DDynaRup/CFactorPaper/InternalReview2/SlipDistribution/';
load([dir,'SlipDtrb_Ns100_Nd50_dip45_BrD0_DipSlip.mat']);
%load([dir,'SlipDtrb_Ns100_Nd50_dip90_BrD0.02_StrikeSlip.mat']);
%load('SlipDtrb_Ns100_Nd50_dip90_BrD100_StrikeSlip.mat')




% figure(21)
% scatter(AlongStrike_coor,AlongDip_coor,20,DipSlip,'filled')
% figure(22)
% scatter(AlongStrike_coor,AlongDip_coor,20,StrikeSlip,'filled')

% interpolated plot

strike_min = min(AlongStrike_coor);
strike_max = max(AlongStrike_coor);
dip_min = min(AlongDip_coor);
dip_max = max(AlongDip_coor);
delta_strike_plot = fault_strike_width/N_strike;
delta_dip_plot = fault_dip_width/N_dip;

[Strike_grid, Dip_grid] = meshgrid([strike_min:delta_strike_plot:strike_max],[dip_min:delta_dip_plot:dip_max]);

F_DipSlip_interp = scatteredInterpolant(AlongStrike_coor,AlongDip_coor,DipSlip);
F_StrikeSlip_interp = scatteredInterpolant(AlongStrike_coor,AlongDip_coor,StrikeSlip);

DipSlip_Interp = F_DipSlip_interp(Strike_grid,Dip_grid);
StrikeSlip_Interp = F_StrikeSlip_interp(Strike_grid,Dip_grid);

%plot
% figure(23)
% imagesc([strike_min strike_max],[dip_min dip_max],DipSlip_Interp)
% set(gca,'YDir','normal')
% figure(24)
% imagesc([strike_min strike_max],[dip_min dip_max],StrikeSlip_Interp)
% set(gca,'YDir','normal')

normalized_factor = max(abs(DipStressChange),abs(StrikeStressChange))/mu*(fault_dip_width);

figure(25)
imagesc([strike_min:delta_strike_plot:strike_max],[dip_min dip_max],sqrt(StrikeSlip_Interp.^2+DipSlip_Interp.^2)./normalized_factor)
hold on
%quiver(AlongStrike_coor,AlongDip_coor,StrikeSlip,DipSlip,'linewidth',1,'Color','k')
contour(Strike_grid,Dip_grid,sqrt(StrikeSlip_Interp.^2+DipSlip_Interp.^2)./normalized_factor,'Color','k','ShowText','on')
set(gca,'YDir','normal')
hold off
set(gca,'Fontsize',18,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear
axis equal


xlabel('Along strike element index')
ylabel('Along dip element index')
colorbar
%caxis([0 1.2])
caxis([0 2.2])
colormap('jet')

%save the data for future plotting, so don't need to calculate again
if abs(DipStressChange)>abs(StrikeStressChange)
    save(['SlipDtrb_Ns',num2str(N_strike),'_Nd',num2str(N_dip),...
    '_dip',num2str(dip),'_BrD',num2str(BurialDepth),'_DipSlip.mat'],...
    'AlongStrike_coor','AlongDip_coor','fault_strike_width','N_strike',...
    'fault_dip_width','N_dip','DipSlip','StrikeSlip',...
    'DipStressChange','StrikeStressChange','mu','dip','BurialDepth');
else
    save(['SlipDtrb_Ns',num2str(N_strike),'_Nd',num2str(N_dip),...
    '_dip',num2str(dip),'_BrD',num2str(BurialDepth),'_StrikeSlip.mat'],...
    'AlongStrike_coor','AlongDip_coor','fault_strike_width','N_strike',...
    'fault_dip_width','N_dip','DipSlip','StrikeSlip',...
    'DipStressChange','StrikeStressChange','mu','dip','BurialDepth');
end

Average_Slip_avg = sqrt(mean(StrikeSlip).^2+mean(DipSlip).^2);
C_k_avg = max(abs(DipStressChange),abs(StrikeStressChange))/mu*(fault_dip_width)/Average_Slip_avg; %k    
    
% Calculated the slip-averge center

H_average_Relative = sum(sum(sqrt(StrikeSlip_Interp.^2+DipSlip_Interp.^2).*(Dip_grid/fault_dip_width)))...
                    ./ sum(sum(sqrt(StrikeSlip_Interp.^2+DipSlip_Interp.^2)));
   
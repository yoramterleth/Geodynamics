%% yoram terleth nov 2021
%% this function considers the stress ont the plate and resulting fracturing and thinning
%% can be implemented as a function using: 
function [PLATE] = f_plate_2D(C, H, theta)

close all

%% for testing
% C.pm = 3000;
% C.pi = 3600; 
% theta =  63 ; 
% C.plate_bottom = 660000 ; 
% C.E = 70e9 ; % Young's modulus ; 
% C.poisson =  .25 ; % poisson ratio
% H = 30000; 
% C.g = 9.81;  % 
% C.xmax =  1000000 ; 
% C.xstep = 1000 ; 
% C.Hstep = 1000 ; 
% 
% % constants 
% C.pw = 1000 ; % density of water in kg m3
% C.fs = 0.6 ; % coefficient of static friction 
% C.dtdh = 25./1000 ; % thermal gradient in the mantle (K m-1)
% C.p_plate = 2700 ; % plithosphere density (kg m3)

%% read in parameters 
% show full zone or only area of tensile stress of interest?
zoom = C.zoom ; 

%% call the w 
flx = flexure(C,H,theta) ; 
%close all
%% initialise a 2 dimensional (along y) plate
x = [0:C.xstep:C.xmax]; 
h = [-H/2 : C.Hstep : H/2] ;
[plate_x,plate_h] = meshgrid(x,h);

%% intialise the 2nd derivate of w over x
% wD2 = (-flx.P0 .* (flx.lambda.^3))./((mean(C.pm(:))-C.pi)* C.g).*exp(-flx.lambda.*plate_x).*(cos(flx.lambda.*plate_x)-sin(flx.lambda.*plate_x)) ; 
wD2 = 2 .* (flx.lambda.^2).*exp(-flx.lambda.*plate_x).*sin(flx.lambda .* plate_x).*((2.* flx.Pb)./((mean(C.pm(:))-C.pi).*C.g)); 
%% calculate the strain in the plate 
Exx = (-plate_h) .* wD2 ; 

%% calculate the stress in the plate 
Sxx = (C.E ./(1-(C.poisson^2))).* Exx ; 

%% identify strenght envelope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plate setup for yield stress
h_pos = [-H : C.Hstep : 0] ; 
if strcmp(C.Anderson, 'full_depth') 
[~,plate_t] = meshgrid(x, h); 
[lith_depth,~] = meshgrid(flx.w, h); 
lith_depth(lith_depth>0)=0 ; % correct potential negative depths
plate_depth = abs(lith_depth + plate_t) ; 

elseif strcmp(C.Anderson, 'linear_660')
plate_d = [-C.plate_bottom : C.plate_bottom/length(flx.w):0];
plate_d = plate_d(1:length(flx.w));
[~,plate_t] = meshgrid(x, h);
[lith_depth,~] = meshgrid(plate_d, h);
%[~,plate_t] = meshgrid(x, h); 
lith_depth(lith_depth>0)=0 ; % correct potential negative depths
plate_depth = abs(lith_depth + plate_t) ; 
%% 
else 
[lith_depth,lith_depth2] = meshgrid(flx.w, h_pos); 
lith_depth2(lith_depth>0)=0 ; % correct potential positive depths
plate_depth = flipud(abs(lith_depth2));% + plate_t) ;
%% 
end 

% hydrostatic fluid pressure 
Pw = C.pw * C.g .* plate_depth ;

% distributed strength envelope 
DeltaSxx = zeros(size(Sxx)) ; 
env_normal = (2.*C.fs.*((C.p_plate.*C.g.*plate_depth)-Pw))./(((1+(C.fs^2)).^(1/2))+C.fs) ; 
DeltaSxx(Sxx<=0) = env_normal(Sxx<=0) ; 
env_thrust = (2.*C.fs.*((C.p_plate.*C.g.*plate_depth)-Pw))./(((1+(C.fs^2)).^(1/2))-C.fs) ; 
DeltaSxx(Sxx>0) = env_thrust(Sxx>0) ; 


%% find the location along x with the max tensile strength 

w_zoom = flx.w ;
w_zoom(Exx(end,:)>0)=nan ; 

% strain: select only upper side of plate, and area that undergoes
% extension
Exx_zoom = Exx; 
Exx_zoom(plate_h<=0)= nan ; 
Exx_zoom(Exx_zoom>0)=nan ; 

% stress: select only upper side of plate, and area that undergoes tensile stress
Sxx_zoom = Sxx; 
Sxx_zoom(plate_h<=0)= nan ; 
Sxx_zoom(Sxx_zoom>0)=nan ; 

% now identify max tensional stress within zone of interest 
% [Sxx_max Sxx_max_ind] = min(Sxx_zoom(:));
[row,col]=find(Sxx_zoom==min(Sxx_zoom(:)));

% strength envelope 
DeltaSxx_zoom = DeltaSxx ; 
DeltaSxx_zoom(plate_h<=0)= nan ; 
DeltaSxx_zoom(Sxx>0)=nan ; 

% indentify areas where strenght envelope is exceeded (by not naning them)
Faulting = ones(size(Sxx)); 
Faulting(abs(DeltaSxx) > abs(Sxx)) = nan; 

Faulting_net = abs(Sxx) - abs(DeltaSxx); 

% consider only area of interest
Faulting_zoom = Faulting ; 
Faulting_zoom(plate_h<=0)= nan ; 
Faulting_zoom(Sxx>0)=nan ; 

% identify the in plate deepest point at which faulting occurs
[~,pd] = meshgrid(x, h_pos); % make array of depth within plate
Faulting_zoom = Faulting_zoom .* pd ; % mutliply logical array with in plate thickness
F_depth = abs(min(Faulting_zoom(:))); % find value of depth: i.e the ammount of lith. material that faulted away 
[F_row, F_col] =find(Faulting_zoom==min(Faulting_zoom(:)));

if strcmp(C.limit_faulting,'fraction')
    if F_depth > F_depth / C.F_limit  
    F_depth = F_depth / C.F_limit; 
    disp(['Limited faulting depth at H/' num2str(C.F_limit) '!'])
    end 
elseif strcmp(C.limit_faulting,'value')
    if F_depth > C.F_limit  
    F_depth = C.F_limit; 
    disp(['Limited faulting depth at ' num2str(C.F_limit) ' m!'])
    end 
end 

% if the built up stresses are too low, there is no faulting
% in this case, set indices and fault penetration depth to zero
if isempty(F_row) 
    F_row = 1 ; 
    disp('Yield stress never exceeded at iteration: no thickness loss.')
else 
    F_row = F_row(1) ; % find location for visual verification 
end
if isempty(F_col)
    F_col= 1 ;
else 
    F_col = F_col(1) ; 
end 
% F_depth(F_depth>H./10) = H./10 ;% 
F_depth(isnan(F_depth)) = 0 ; 


% zoom for plotting  
Faulting_net_zoom = Faulting_net ; 
Faulting_net_zoom(plate_h<=0)= nan ; 
Faulting_net_zoom(Sxx>0)=nan ; 


%% Save variables of interest

PLATE.Exx.zoom = Exx_zoom ;
PLATE.Exx.full = Exx ;
PLATE.Sxx.zoom = Sxx_zoom ;
PLATE.Sxx.full = Sxx ;
PLATE.DeltaSxx.zoom = DeltaSxx_zoom ; 
PLATE.Faulting_net.zoom = Faulting_net_zoom ; 
PLATE.Faulting_depth = F_depth ; 
PLATE.wD2 = wD2 ; 
PLATE.x = x ; 
PLATE.h = h ; 
PLATE.F_col = F_col ; 
PLATE.F_row = F_row ; 


%% visualization 
if zoom 
    Sxx_plot = Sxx_zoom ; 
    Exx_plot = Exx_zoom ;
    w_plot = w_zoom ; 
    DeltaSxx_plot = DeltaSxx_zoom ; 
    Faulting_plot = Faulting_net_zoom ; 
    
else 
    Sxx_plot = Sxx ; 
    Exx_plot = Exx ; 
    w_plot = flx.w ; 
    DeltaSxx_plot = DeltaSxx ; 
    Faulting_plot = Faulting_net ; 
end 

PLATE.DeltaSxx_plot = DeltaSxx_plot ; 
PLATE.Sxx_plot = Sxx_plot ; 
PLATE.Faulting_plot = Faulting_plot ; 
PLATE.Exx_plot = Exx_plot ; 
PLATE.w_plot = w_plot ; 

% figure
% set(gcf,'Visible', 'off');
% % re-plot profile for reference
% ax(1) = subplot(3,1,1); 
% plot(flx.x./1000, w_plot./1000,'linewidth',3,'Color',[0.6784    0.5059    0.3137])
% grid on 
% xlabel('Distance along plate [km]'), ylabel('w [km]')
% %xlim([0 950])
% % 
% % % distance from mid-plate
% % ax(2) = subplot(4,1,2);
% % pcolor(x./1000,h./1000,abs(plate_h)./1000), shading flat 
% % colormap(ax(2),flipud(viridis(264*2)))
% % xlabel('Distance along plate [km]')
% % ylabel('Depth in plate [km]')
% % c = colorbar;
% % ylabel(c,'Distance from mid-plate [km]');
% 
% % strain in plate
% ax(3) = subplot(3,1,2);
% pcolor(x./1000,h./1000,Exx_plot), shading flat, hold on 
% contour(x./1000,h./1000,Exx,[0 0],'--k','linewidth',1)
% colormap(ax(3),flipud(viridis(264*2)))
% xlabel('Distance along plate [km]')
% ylabel('Distance from mid-plate [km]')
% c = colorbar;
% ylabel(c,'Longitudinal Strain');
% %caxis([-1000 1000])
% %xlim([0 910])
% box on 
% 
% % stress in plate 
% ax(4) = subplot(3,1,3); 
% pcolor(x./1000,h./1000,Sxx_plot), shading flat ,hold on 
% contour(x./1000,h./1000,Sxx,[0 0],'--k','linewidth',1)
% colormap(ax(4),flipud(viridis(264*2)))
% xlabel('Distance along plate [km]')
% ylabel('Distance from mid-plate [km]')
% c = colorbar;
% ylabel(c,'Longitudinal stress');
% %xlim([0 910])
% box on 
% % scatter(x(col)./1000,h(row)./1000,100,'sr','LineWidth',2)
% 
% figure 
% set(gcf,'Visible', 'off');
% ax(1) = subplot(2,1,1); 
% pcolor(x./1000,h./1000,DeltaSxx_plot), shading flat ,hold on 
% contour(x./1000,h./1000,Sxx,[0 0],'--k','linewidth',1)
% colormap(ax(1),flipud(viridis(264*2)))
% xlabel('Distance along plate [km]')
% ylabel('Distance from mid-plate [km]')
% c = colorbar;
% ylabel(c,'Anderson Yield Strength');
% %xlim([0 910])
% box on 
% ax(2) = subplot(2,1,2); 
% pcolor(x./1000,h./1000,Faulting_plot), shading flat ,hold on 
% contour(x./1000,h./1000,Sxx,[0 0],'--k','linewidth',1)
% colormap(ax(2),plasma(264*2))
% xlabel('Distance along plate [km]')
% ylabel('Distance from mid-plate [km]')
% c = colorbar;
% ylabel(c,'|\sigma_{xx}| - |\sigma_{Yield}|');
% %xlim([0 910])
% box on
% contour(x./1000,h./1000,Faulting_net,[0 0],'--g','linewidth',1.6)
% scatter(x(F_col)./1000,h(F_row)./1000,100,'sr','LineWidth',2)


 end 


%% SB main 
% Terleth nov. 2021 for geodynamics 
%% INITIALISE

close all 
clearvars
clc

% read int constants 
[C]=parameters(); 

% set intitial plate thickness
H = 70000 ; 

% adjust plate length 
C.xmax = C.plate_bottom ./ sin(C.t_b) ; % the length of the subducting plate really depends on theta... from buffet 2006 

% initialise arrays to save the moments 
T_g = zeros(length(C.t_b),length(C.depth));
T_res =  zeros(length(C.t_b),length(C.depth));

% initialise array storing efective plate thickness over time
H_store = zeros(1,length(C.time)); 

% initialiase array storing subduction angle over time
Theta_store = zeros(1,length(C.time)); 

%% "TIME" ITERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1 : length(C.time)


% intiate x for length iterations 
% = [0:C.xstep:C.xmax]; 
    
    
for i = 1:length(C.depth) -1
    
% adjust mu 
mu = C.mu(i) ; 

% adjust local mantle density 
pm = C.pm(i) ; 

% calculate density difference mantle/oc. lithosphere
delta_p = pm - C.p_plate; 

% adjust plate length 
l_full = C.depth(end)./ sin(C.t_b) ; % the length of the subducting plate really depends on theta... from Buffet (2006) 
l = l_full - (C.depth(i)./sin(C.t_b));  % l is lthe length of lith below any point
x = (C.depth(i)./sin(C.t_b)); % x is the position of the plate we are at

% get the position along z for each considered angle at the current x 
flx = flexure(C,H,rad2deg(C.t_b)) ; 

%% TORQUE 

% torque due to flow 
T_f = 2 .* mu .* C.v .* l .* (((sin(C.t_b)./((pi -C.t_b)+ sin(C.t_b)))) + ((sin(C.t_b).^2)./(C.t_b.^2-(sin(C.t_b).^2)))); 

% torque due to gravity 
T_g(:,i) = (1/2) .* delta_p .* C.g .* H .* (l.^2) .* cos(C.t_b) ; 

% torque due to elastic plate bending 
D = (C.E * (H ^3))/(12*(1-(C.poisson^2))) ; % flexural rigidity
% 2nd derivative over x of w  

% wD2 = (-flx.P0 .* (flx.lambda.^3))/((pm-C.pi).* C.g).*exp(-flx.lambda.*x).*(cos(flx.lambda.*x)-sin(flx.lambda.*x)) ; 
wD2 = 2 .* (flx.lambda.^2).*exp(-flx.lambda.*x).*sin(flx.lambda .*x).*((2.* flx.Pb)./((mean(C.pm(:))-C.pi).*C.g));
% and the moment at this location % the minus is adjustment for the change
% of reference frame
T_b = D .* wD2 ; 

% add flexure to the flow side, as it likely reduces angle of subduction 
T_res(:,i)= T_f + T_b;

end % end of looping over the depths

%% now sum over the depths 
T_g_total = sum(T_g,2); 
T_res_total = sum(T_res,2); 

%% solve graphically and find min angle 
compare = rmmissing([T_g_total T_res_total]); 
[mindif, min_ind] = min(abs((compare(:,1)-compare(:,2))));
t_b_min = rad2deg(C.t_b(min_ind)) ;
disp(['Theta b is ' num2str(t_b_min) ' degrees.'])

 
%% ADJUSTMENT OF H DUE TO FAULTING IN THE PLATE 

% initialise PLATE 

% call 2d plate consideration  
C.xmax = l_full ; % consider only plate lenght eq to t_b_min
[PLATE] = f_plate_2D(C, H,t_b_min) ; 
C.xmax = C.plate_bottom ./ sin(C.t_b) ; % reset to original

% Store the effective pate thickness used in the finished iteration
H_store(t) = H ;

% Adjust new plate thickness 
H = H - PLATE.Faulting_depth ; 

%% save & visdualise change in theta and in H 
%% plot the profile 
% 
figure 
set(gcf,'Visible', 'off');
plot(PLATE.x./1000, PLATE.w_plot./1000,'linewidth',2,'Color',[0.6784    0.5059    0.3137])
grid on 
xlabel('Distance from w_{min} along x [km]'), ylabel('Displacement along w of mid-plate [km]')
title(['Effective H = ' num2str(H./1000),' km; dip angle = ' num2str(t_b_min) ' degrees.'])
F = getframe(gcf) ;
Frame_w(t)= F;


% final torque balance 
figure
set(gcf,'Visible', 'off');
plot(rad2deg(C.t_b),T_g_total,'linewidth',1.2); hold on ;
plot(rad2deg(C.t_b),T_res_total,'linewidth',1.2);
xlabel('\theta_{b} [^{\circ}]')
ylabel('Torque [N m^{-1}]')
legend('T_{gravity}','T_{resistive}')
grid on 
title(['Effective H = ' num2str(H./1000),' km; dip angle = ' num2str(t_b_min) ' degrees.'])
ylim([0,10e20])
F = getframe(gcf) ;
Frame_torque(t)= F;

% store the value of t_b_min 
Theta_store(t) = t_b_min ;

%% visualise the colorplots

figure
set(gcf,'Visible', 'off');
% re-plot profile for reference
ax(1) = subplot(3,1,1); 
plot(PLATE.x./1000, PLATE.w_plot.*1000,'linewidth',3,'Color',[0.6784    0.5059    0.3137])
grid on 
xlabel('Distance along plate [km]'), ylabel('w [km]')
xlim([0 1500])
% 
% % distance from mid-plate
% ax(2) = subplot(4,1,2);
% pcolor(x./1000,h./1000,abs(plate_h)./1000), shading flat 
% colormap(ax(2),flipud(viridis(264*2)))
% xlabel('Distance along plate [km]')
% ylabel('Depth in plate [km]')
% c = colorbar;
% ylabel(c,'Distance from mid-plate [km]');

% strain in plate
ax(3) = subplot(3,1,3);
% pcolor(PLATE.x./1000,PLATE.h./1000,PLATE.Exx_plot), shading flat, hold on 
% contour(PLATE.x./1000,PLATE.h./1000,PLATE.Exx_plot,[0 0],'--k','linewidth',1)
% colormap(ax(3),flipud(viridis(264*2)))
% xlabel('Distance along plate [km]')
% ylabel('h [km]')
% c = colorbar;
% ylabel(c,'Longitudinal Strain');
% %caxis([-1000 1000])
% xlim([0 1500])
% box on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcolor(PLATE.x./1000,PLATE.h./1000,PLATE.Faulting_plot), shading flat ,hold on 
contour(PLATE.x./1000,PLATE.h./1000,PLATE.Sxx_plot,[0 0],'--k','linewidth',1)
colormap(ax(3),plasma(264*2))
xlabel('Distance along plate [km]')
ylabel('h[km]')
c = colorbar;
ylabel(c,'|\sigma_{xx}| - |\sigma_{Yield}|');
xlim([0 1500])
box on
contour(PLATE.x./1000,PLATE.h./1000,PLATE.Faulting_plot,[0 0],'--g','linewidth',1.6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stress in plate 
ax(4) = subplot(3,1,2); 
pcolor(PLATE.x./1000,PLATE.h./1000,PLATE.Sxx_plot), shading flat ,hold on 
contour(PLATE.x./1000,PLATE.h./1000,PLATE.Sxx_plot,[0 0],'--k','linewidth',1)
colormap(ax(4),flipud(viridis(264*2)))
xlabel('Distance along plate [km]')
ylabel('h [km]')
c = colorbar;
ylabel(c,'Longitudinal stress');
xlim([0 1500])
box on 
% scatter(x(col)./1000,h(row)./1000,100,'sr','LineWidth',2)
F = getframe(gcf) ;
Frame_stress(t)= F;

figure 
ax(3) = subplot(3,1,1); 
plot(PLATE.x./1000, PLATE.w_plot./1000,'linewidth',3,'Color',[0.6784    0.5059    0.3137])
grid on 
xlabel('Distance along plate [km]'), ylabel('w [km]')
xlim([0 1500])
set(gcf,'Visible', 'off');
ax(1) = subplot(3,1,2); 
pcolor(PLATE.x./1000,PLATE.h./1000,PLATE.DeltaSxx_plot), shading flat ,hold on 
contour(PLATE.x./1000,PLATE.h./1000,PLATE.Sxx_plot,[0 0],'--k','linewidth',1)
colormap(ax(1),flipud(viridis(264*2)))
xlabel('Distance along plate [km]')
ylabel('h [km]')
c = colorbar;
ylabel(c,'Anderson Yield Strength');
xlim([0 1500])
box on 
ax(2) = subplot(3,1,3); 
pcolor(PLATE.x./1000,PLATE.h./1000,PLATE.Faulting_plot), shading flat ,hold on 
contour(PLATE.x./1000,PLATE.h./1000,PLATE.Sxx_plot,[0 0],'--k','linewidth',1)
colormap(ax(2),plasma(264*2))
xlabel('Distance along plate [km]')
ylabel('h[km]')
c = colorbar;
ylabel(c,'|\sigma_{xx}| - |\sigma_{Yield}|');
xlim([0 1500])
box on
contour(PLATE.x./1000,PLATE.h./1000,PLATE.Faulting_plot,[0 0],'--g','linewidth',1.6)
%scatter(PLATE.x(PLATE.F_col)./1000,PLATE.h(PLATE.F_row)./1000,100,'sr','LineWidth',2)

F = getframe(gcf) ;
Frame_yield(t)= F;


end 

%% Visualize iteration evoltuion 
figure 
yyaxis left 
plot(C.time,Theta_store,'linewidth',1.6)
xlabel('Model iterations')
ylabel('Subd. Angle \theta [^{\circ}]')
grid on 
ylim([min(Theta_store)-1, max(Theta_store)+1]);

yyaxis right
plot(C.time,H_store./1000,'linewidth',1.6)
ylabel('Effective lith. thickness [km]')
ylim([min(H_store./1000)-2, max(H_store./1000)+2]);

%% figure of the mantle properties 
figure 
subplot(1,2,1)
plot(C.mu, C.depth./1000,'linewidth',1.5,'color','#77AC30')
ylabel('depth [km]'),xlabel('mantle viscosity [Pa s]')
axis ij 
set(gca,'XAxisLocation','top')
grid on 
yline(410,'--r','linewidth',1.6)
xlim([min(C.mu)-1e19,max(C.mu)+1e19])
subplot(1,2,2)
plot(C.pm, C.depth./1000,'linewidth',1.5,'color','#7E2F8E')
xlabel('Mantle density [kg m^{-3}]')
axis ij 
grid on 
yline(410,'--r','410 boundary','linewidth',1.6)
xlim([min(C.pm)-100,max(C.pm)+100])
set(gca,'XAxisLocation','top')

%% make videos 
videos(Frame_torque,'v_torque_balance'); 
videos(Frame_stress, 'V_stress'); 
videos(Frame_yield,'V_yieldandH');
videos(Frame_w,'V_w'); 


% function computing the behavior of an elastic plate subject to a certain
% dip angle
%yoram terleth - nov 2021

%% expects: 
% - C.pm [kg m3] the density of the mantle and 
% - C.pi [kg m3] the density of the stuff above the plate
% - C.g grav constant 
% - C.xmax [m] model doain in x direction 
% - C.xstep [m] model resolution in x direction
% - C.poisson Poisson ration
% - C.E Young shear modulus
% - theta [deg] the dip angle of the subducting elasticv plate
% - C.plate_bottom [ m] depth of the lower model limit or lower plate limit 
% - H [m] ingoing plate thikness 

function [flx] = flexure(C,H,theta)

%% for testing
% C.pm = 3000;
% C.pi = 3600; 
% theta = 63 ; 
% C.plate_bottom = 660000 ; 
% C.E = 70e9 ; % Young's modulus ; 
% C.poisson =  .25 ; % poisson ratio
% H = 70000; 
% C.g = 9.81;  % 
% C.xmax = 1000000 ; 
% C.xstep = 1000 ; 
%% initialise x 
x = [0:C.xstep: C.xmax];

%% calculate lambda the distance from surface to lowest point 
lambda = (pi * tan(deg2rad(theta))) ./ (2 * C.plate_bottom) ; 



%% calculate w over x (the actual position of the plate)

% felxural rigidity 
D = (C.E * (H^3))/(12*(1-(C.poisson^2))) ; 

% P0 from notes eric 
P0 = ((2*(mean(C.pm) - C.pi)*C.g * D)./lambda ); %/1e10 ;
%Pb = P0 ./ 4 ; 
% Pb from notes and rearranging anf ymax at x = 0 boundary condition
Pb = (C.plate_bottom .* (mean(C.pm(:))-C.pi) .* C.g)./2 ; 


% vertical plate displacement 
%w = -((P0 .* lambda) ./ (2 .* ( mean(C.pm(:)) - C.pi) .* C.g)) .* exp(-lambda .* x) .* (cos(lambda.*x) + sin(lambda.*x)); 
w = -((2 .* Pb .* lambda) ./ ((mean(C.pm(:))-C.pi).*C.g)) .* exp(-lambda .* x) .* cos(lambda.*x);

%% plot the profile 
% 
% figure(1) 
% plot(x./1000, w./1000,'linewidth',.4)%,'Color',[0.6784    0.5059    0.3137]),hold on
% hold on
% grid on 
% xlabel('Distance from w_{min} along x [km]'), ylabel('Displacement along w of mid-plate [km]')
% %ylim([-660,10])

%% store the relevant output
flx.lambda = lambda ; 
flx.Pb = Pb ; 
flx.w = w; 
flx.D = D ; 
flx.P0 = P0 ; 
flx.x = x ; 


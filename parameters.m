function [C] = parameters()

%% gives all the parameters & constants

%% model parameters: large loop  

% time steps to consider
C.time =  [1:5] ; 

% the range of theta b angles we'd like to try 
C.t_b = deg2rad([10:1:80])';

C.depth = [0:50000:660000] ;                                                    

% do we want to consider the strength envelope from the actual depth in the
% mantle, or just the depth within the plate. If we use the full mantle,
% then the depth profile dominates... 
C.Anderson = 'linear_660' ; % or 'full_depth' or empty for just h 
%% model parameters: plate internal 

% model domain in x direction 
C.xmax =  660000 ; 
C.xstep = 1000 ;

% internal to plate: step in thicness considered
C.Hstep = 1000 ; 

% zoom in plots? 
C.zoom = 0 ; 

% limit to the faulting extent 
C.limit_faulting = 'value' ;    % 'fraction' to limit faulting by a fraction of the current plate thickness
                                % 'value' to limit faulting to a set value
                                % in m 
                                % 'no' to not limit faulting 
C.F_limit = 1000 ;                % limit value to faulting: 
                                % - if 'fraction' should be a positive integer to divide plate thickness by
                                % - if 'value' should be the max fault
                                % length in m. 

%% lithosphere parameters

% plate velocity [m s-1]
C.v = 7.9e-11 ; % 1.26839e-9  


%C.h = 20000 ; % 2e4 ;                                                         
C.Rmin = 200000 ; 

% infill density: continental lithosphere
C.pi = 2600; 

C.plate_bottom = 660000 ; 
C.E = .5e11 ;% ; .6e11 ;% .5e3 ; % Young's modulus ; 
C.poisson =  .25 ; % poisson ratio
%H = 30000; 

C.xmax =  1000000 ; % plate length 
C.xstep = 1000 ; % horizontal model step
C.Hstep = 1000 ; % vertical model step

C.fs = 0.6 ; % coefficient of static friction 
C.p_plate = 2700 ; % plithosphere density (kg m3)


%% mantle parameters 

% mantle viscosity - depth dependent 
C.mu = zeros(length(C.depth),1);
% C.mu((C.plate_bottom .* cos(deg2rad(90-rad2deg(C.t_b)))) <410000) = 10^22 ;
% C.mu((C.plate_bottom .* cos(deg2rad(90-rad2deg(C.t_b))))>410000) = 5*10^(19); %20.5); 
C.mu(C.depth <= 410000) = 10^19 ; 
C.mu(C.depth>410000) = 10^23 ; 

% mantle density 
C.pm = zeros(length(C.depth),1) ; %3000
%C.pm((C.plate_bottom .* cos(deg2rad(90-rad2deg(C.t_b)))) <410000) = 2800 ; % [kg/m3] olivine dens
%C.pm((C.plate_bottom .* cos(deg2rad(90-rad2deg(C.t_b)))) >410000) = 3840 ; % [kg/m3] wadsleyite dens
C.pm(C.depth <=410000) = 2800 ; 
C.pm(C.depth>410000) = 3850 ; 
% constants 
C.pw = 1000 ; % density of water in kg m3

%% general constants
C.g = 9.81;  % grav acceleration [m s-2]

% coordinate vector length
C.r = 1 ; 








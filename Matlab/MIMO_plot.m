%% Transfer Function Approach (using matlab tf and bode commands)
clear; close all; clc;
format short;

%% Define Variables
% DOF numbering: 1 = swing, 2 = flywheel 1, 3 = flywheel 2 (2-pitch only)

J1 = 22.45; % kg*m^2, swing moment of inertia
K1 = 517;   % Nm/rad, swing restoring stiffness
C1 = 1.65;  % Nm*s/rad, swing friction
J2 = 4.19;  % kg*m^2, flywheel moment of inertia
K2 = 40.75; % Nm/rad, magnetic spring stiffness
C2 = 0.85;  % Nm*s/rad, effective system damping
J3 = J2;    % second flywheel identical to first
K3 = K2;
C3 = C2;
Na = 3/1;   % Gearbox transmission ratio
Nb = 28/40; % Belt drive transmission ratio
Cgen = 0;    % Generator damping (Nm*s/rad)
Cpto = Na^2*Cgen;
texc_amp = 65;  % Excitation torque amplitude (Nm)

%% Single resonator
% G11 = Excitation torque to swing velocity
% G12 = Gearbox torque to swing velocity
% G21 = Excitation torque to flywheel velocity
% G22 = Gearbox torque to flywheel velocity
% G21 = G12... Gearbox torque to swing velocity = Excitation torque to flywheel velocity

% Define laplace variable
s = tf('s') ; 

% Define impedances
z1 = J1*s + C1 + (K1 /s) ;
z2 = J2*s ;
z3 = C2 + ((Nb^2 * K2)/s) ;

% Define inputs
G11 = (z2 + z3) / (z1*z2 + z1*z3 + z2*z3) ;
G12 = -z2 / (z1*z2 + z1*z3 + z2*z3) ;
G21 = G12;  % Coupled
G22 = (z1 + z2) / (z1*z2 + z1*z3 + z2*z3) ;

G_mimo = [G11, G12; 
          G21, G22];

% Plot the Bode diagram for all input-output pairs
figure(1);
bode(G_mimo);
grid on;

%% Double resonator
% G11 = Excitation torque to swing velocity
% G12 = Gearbox torque 2 to swing velocity
% G13 = Gearbox torque 3 to swing velocity

% G21 = Excitation torque to flywheel 2 velocity
% G22 = Gearbox torque 2 to flywheel 2 velocity
% G23 = Gearbox torque 3 to flywheel 2 velocity

% G31 = Excitation torque to flywheel 3 velocity
% G32 = Gearbox torque 2 to flywheel 3 velocity
% G33 = Gearbox torque 3 to flywheel 3 velocity

% G21 = G12... Excitation torque to flywheel 2 velocity = Gearbox torque 2 to swing velocity
% G23 = G32... Gearbox torque 3 to flywheel 2 velocity = Gearbox torque 2 to flywheel 3 velocity
% G13 = G31... Gearbox torque 3 to swing velocity =  Excitation torque to flywheel 3 velocity


% Define laplace variable
s = tf('s') ; 

% Define impedances
z1 = J1*s + C1 + (K1/s) ;
z2 = J2*s ;
z3 = J3*s ;
z4 = C2 + ((Nb^2 * K2)/s) ;
z5 = C3 + ((Nb^2 * K3)/s) ;

% Define inputs
denom = (z2*z4*z5 + z3*z4*z5 + z1*(z2 + z4)*(z3 + z5) + z2*z3*(z4 + z5)) ;
G11 = ((z2 + z4)*(z3 + z5)) / denom ;
G12 = (-z2*(z3+z5)) / denom ;
G13 = (z3*(z2 + z4)) / denom ;
G21 = G12 ;
G22 = (z3*z5 + z1*(z3 + z5) + z2*(z3 + z5)) / denom;
G23 = (z2*z3) / denom ;
G31 = G13;
G32 = G23;
G33 = (z3*z4 + z1*(z2 + z4) + z2*(z3 + z4)) / denom;

G_mimo = [G11, G12, G13; 
          G21, G22, G23;
          G31, G32, G33];

% Plot the Bode diagram for all input-output pairs
figure(2);
bode(G_mimo);
grid on;

%% Plot comparison
%
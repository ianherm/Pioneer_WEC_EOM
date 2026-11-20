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

% Frequency vector shared by all plots (Hz)
f = logspace(log10(0.05), log10(20), 3000);

red  = [1 0 0];   
blue = [0 0 1];   

%% Single resonator (use G for MIMO)
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

% Bode diagram (Hz) for all input-output pairs
figure(1);
plotBodeGrid(G_mimo, 'G', f, blue);

%% Double resonator (use H for MIMO)
% H11 = Excitation torque to swing velocity
% H12 = Gearbox torque 2 to swing velocity
% H13 = Gearbox torque 3 to swing velocity

% H21 = Excitation torque to flywheel 2 velocity
% H22 = Gearbox torque 2 to flywheel 2 velocity
% H23 = Gearbox torque 3 to flywheel 2 velocity

% H31 = Excitation torque to flywheel 3 velocity
% H32 = Gearbox torque 2 to flywheel 3 velocity
% H33 = Gearbox torque 3 to flywheel 3 velocity

% H12 = H21, H31 = H13, identical flywheels -> H12 = H21 = H31 = H13
% H22 = H33
% H23 = H32 (unique to double resonator)

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
H11 = ((z2 + z4)*(z3 + z5)) / denom ;
H12 = (-z2*(z3+z5)) / denom ;
H13 = (-z3*(z2 + z4)) / denom ;
H21 = H12 ;
H22 = (z3*z5 + z1*(z3 + z5) + z2*(z3 + z5)) / denom;
H23 = (z2*z3) / denom ;
H31 = H13;
H32 = H23;
H33 = (z3*z4 + z1*(z2 + z4) + z2*(z3 + z4)) / denom;

H_mimo = [H11, H12, H13; 
          H21, H22, H23;
          H31, H32, H33];

% Plot Bode
figure(2);
plotBodeGrid(H_mimo, 'H', f, red);

%% Plot comparison (H = red, G = blue)

% Figure 3: H11 v G11
figure(3);
plotCompare({H11, G11}, {red, blue}, {'H11', 'G11'}, ...
    'H11 v G11: Excitation Torque to Swing Velocity', f);

% Figure 4: H21/H12/H31/H13 v G12/G21
% (all identical within each group, so one curve per group is plotted)
figure(4);
plotCompare({H12, G12}, {red, blue}, {'H21/H12/H31/H13', 'G12/G21'}, ...
    ['H21/H12/H31/H13 v G12/G21: Curves represent Excitation Torque to Flywheel(s) Velocity & ' ...
     'Gear Box(s) Torque to Swing Velocity'], f);

% Figure 5: H23/H32 (single profile)
figure(5);
plotCompare({H23}, {red}, {'H23/H32'}, ...
    'H23/H32: Gearbox Torque 3 to Flywheel Velocity 2, Gearbox Torque 3 to Flywheel Velocity 2', f);

% Figure 6: G22 v H22/H33
figure(6);
plotCompare({H22, G22}, {red, blue}, {'H22/H33', 'G22'}, ...
    'H22/H33 v G22: Gearbox Torque to Flywheel Velocity, Same resonator', f);

%% Local functions 
function [magdB, phdeg] = bodeHz(sys, f)
    % Frequency response evaluated at f (Hz)
    resp  = squeeze(freqresp(sys, 2*pi*f));
    magdB = 20*log10(abs(resp));
    phdeg = unwrap(angle(resp)) * 180/pi;
end

function plotBodeGrid(M, name, f, col)
    if nargin < 4, col = [0 0.447 0.741]; end   % default MATLAB blue
    n = size(M,1);  m = size(M,2);
    tiledlayout(2*n, m, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:n
        for j = 1:m
            [magdB, phdeg] = bodeHz(M(i,j), f);

            nexttile((2*i-2)*m + j);               % magnitude tile
            plot(f, magdB, 'Color', col, 'LineWidth', 1.2); grid on;
            title(sprintf('%s%d%d', name, i, j));
            ylabel('Magnitude (dB)'); xlim([0.10 1]); ylim([-80 20]);

            nexttile((2*i-1)*m + j);               % phase tile
            plot(f, phdeg, 'Color', col, 'LineWidth', 1.2); grid on;
            ylabel('Phase (deg)'); xlabel('Frequency (Hz)');
            xlim([0.10 1]); ylim([-360 360]); yticks(-360:180:360);
        end
    end
end

function plotCompare(sysList, colors, labels, ttl, f)
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ttl, 'FontWeight', 'bold');
    axM = nexttile; hold on; grid on;
    axP = nexttile; hold on; grid on;
    for k = 1:numel(sysList)
        [magdB, phdeg] = bodeHz(sysList{k}, f);
        plot(axM, f, magdB, 'Color', colors{k}, 'LineWidth', 1.4, ...
            'DisplayName', labels{k});
        plot(axP, f, phdeg, 'Color', colors{k}, 'LineWidth', 1.4, ...
            'DisplayName', labels{k});
    end
    set([axM axP], 'XLim', [0.1 1]);
    set(axM, 'YLim', [-80 20]);
    set(axP, 'YLim', [-360 360], 'YTick', -360:180:360);
    ylabel(axM, 'Magnitude (dB)');
    ylabel(axP, 'Phase (deg)'); xlabel(axP, 'Frequency (Hz)');
    legend(axM, 'Location', 'best');
end
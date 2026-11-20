clear; close all; clc;
format short;

%% Compare response of one and two pitch-resonator swing-sets
% DOF numbering: 1 = swing, 2 = flywheel 1, 3 = flywheel 2 (2-pitch only)

%% Define Variables
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
Cgen = 0.20;    % Generator damping (Nm*s/rad)
Cpto = Na^2*Cgen;
texc_amp = 65;  % Excitation torque amplitude (Nm)

%% 1-Pitch Resonator: 2 DOF

J_1pitch = [J1+J2, J2;
            J2,    J2];

C_1pitch = [C1, 0;
            0, C2+Cpto];

K_1pitch = [K1, 0;
            0, K2*Nb^2];

F_1pitch = [texc_amp; 0];

%% 2-Pitch Resonator: 3 DOF

J_2pitch = [J1+J2+J3, J2, J3;
            J2,       J2, 0;
            J3,       0,  J3];

C_2pitch = [C1,      0,       0;
            0, C2+Cpto,       0;
            0,      0, C3+Cpto];

K_2pitch = [K1, 0,       0;
            0, K2*Nb^2, 0;
            0, 0, K3*Nb^2];

F_2pitch = [texc_amp; 0; 0];
%% Calculate Frequency Responses
% (suffix = which system; row index of each result = DOF 1, 2, 3)

w = linspace(0.1, 8, 4000); % rad/s
[amp_1pitch,ph_1pitch] = frequencyResponse(J_1pitch,C_1pitch,K_1pitch,F_1pitch,w);
[amp_2pitch,ph_2pitch] = frequencyResponse(J_2pitch,C_2pitch,K_2pitch,F_2pitch,w);
%% Plot Amplitude Comparison

figure(1)
% Swing angles (DOF 1)
subplot(2,1,1)
semilogy(w,amp_1pitch(1,:),'b-','LineWidth',2)
hold on
semilogy(w,amp_2pitch(1,:),'r--','LineWidth',2)
grid on
ylabel('Swing Angle Amplitude (deg)')
title('Frequency Response: 1-Pitch vs 2-Pitch')
legend('1-Pitch Swing (\theta_1)','2-Pitch Swing (\theta_1)','10 deg limit')
set(gca,'FontSize',12,'LineWidth',1.5)
% Relative flywheel angles (DOF 2, 3)
subplot(2,1,2)
semilogy(w,amp_1pitch(2,:),'b-','LineWidth',2)
hold on
semilogy(w,amp_2pitch(2,:),'r--','LineWidth',2)
semilogy(w,amp_2pitch(3,:),'Color',[1 0.871 0.282], ...
    'LineStyle',':','LineWidth',2)
grid on
xlabel('\omega (rad/s)')
ylabel('Relative Flywheel Angle Amplitude (deg)')
legend('1-Pitch Flywheel (\theta_2)','2-Pitch Flywheel 1 (\theta_2)', ...
    '2-Pitch Flywheel 2 (\theta_3)')
set(gca,'FontSize',12,'LineWidth',1.5)
%% Plot Phase Comparison

figure(2)
% Swing angles (DOF 1)
subplot(2,1,1)
plot(w,ph_1pitch(1,:),'b-','LineWidth',2)
hold on
plot(w,ph_2pitch(1,:),'r--','LineWidth',2)
grid on
ylabel('Swing Phase (deg)')
title('Frequency Response Phase')
legend('1-Pitch Swing (\theta_1)','2-Pitch Swing (\theta_1)')
set(gca,'FontSize',12,'LineWidth',1.5)
% Relative flywheel angles (DOF 2, 3)
subplot(2,1,2)
plot(w,ph_1pitch(2,:),'b-','LineWidth',2)
hold on
plot(w,ph_2pitch(2,:),'r--','LineWidth',2)
plot(w,ph_2pitch(3,:),'Color',[1 0.871 0.282], ...
    'LineStyle',':','LineWidth',2)
grid on
xlabel('\omega (rad/s)')
ylabel('Relative Flywheel Phase (deg)')
legend('1-Pitch Flywheel (\theta_2)','2-Pitch Flywheel 1 (\theta_2)', ...
    '2-Pitch Flywheel 2 (\theta_3)')
set(gca,'FontSize',12,'LineWidth',1.5)
%% Local frequency response function

function [amp,phase,X] = frequencyResponse(J,C,K,F,w)
    % Calculate frequency response
    X = zeros(size(J,1),length(w));
    for k = 1:length(w)
        X(:,k) = (K-w(k)^2*J+1i*w(k)*C)\F;
    end
    % Convert amplitude and phase to degrees
    amp = rad2deg(abs(X));
    phase = rad2deg(angle(X));
end
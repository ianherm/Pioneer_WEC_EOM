
clear; close all; clc;
format short;

%% Compare response of one and two pitch-resonator swing-sets

%% Define Variables

Js = 22.45; % kg*m^2, swing moment of inertia
Ks = 517;   % Nm/rad, swing restoring stiffness
Cs = 1.65;  % Nm*s/rad, swing friction

J1 = 4.19;  % kg*m^2, flywheel moment of inertia
K1 = 40.75; % Nm/rad, magnetic spring stiffness
C1 = 0.85;  % Nm*s/rad, effective system damping

J2 = J1;
K2 = K1;
C2 = C1;

Na = 3/1;   % Gearbox transmission ratio
Nb = 28/40; % Belt drive transmission ratio

Cgen = 0.20;    % Generator damping (Nm*s/rad)
Cpto = Na^2*Cgen;

texc_amp = 65;  % Excitation torque amplitude (Nm)

%% 1-Pitch Resonator: 2 DOF

J_1pitch = [Js+J1, J1;
            J1,    J1];

C_1pitch = [Cs, 0;
            0, C1+Cpto];

K_1pitch = [Ks, 0;
            0, K1*Nb^2];

F_1pitch = [texc_amp; 0];

%% 2-Pitch Resonator: 3 DOF

J_2pitch = [Js+J1+J2, J1, J2;
            J1,       J1, 0;
            J2,       0,  J2];

C_2pitch = [Cs,      0,       0;
            0, C1+Cpto,       0;
            0,      0, C2+Cpto];

K_2pitch = [Ks, 0,       0;
            0, K1*Nb^2, 0;
            0, 0, K2*Nb^2];

F_2pitch = [texc_amp; 0; 0];

%% Calculate Frequency Responses

w = linspace(0.1, 8, 4000); % rad/s

[amp_1,ph_1] = frequencyResponse(J_1pitch,C_1pitch,K_1pitch,F_1pitch,w);
[amp_2,ph_2] = frequencyResponse(J_2pitch,C_2pitch,K_2pitch,F_2pitch,w);

%% Plot Amplitude Comparison

figure(1)

% Swing angles
subplot(2,1,1)
semilogy(w,amp_1(1,:),'b-','LineWidth',2)
hold on
semilogy(w,amp_2(1,:),'r--','LineWidth',2)
yline(10,'--k','10 deg limit')
grid on
ylabel('Swing Angle Amplitude (deg)')
title('Frequency Response: 1-Pitch vs 2-Pitch')
legend('1-Pitch Swing','2-Pitch Swing','10 deg limit')
set(gca,'FontSize',12,'LineWidth',1.5)

% Relative flywheel angles
subplot(2,1,2)
semilogy(w,amp_1(2,:),'b-','LineWidth',2)
hold on
semilogy(w,amp_2(2,:),'r--','LineWidth',2)
semilogy(w,amp_2(3,:),'Color',[0.2 0.6 0.2], ...
    'LineStyle',':','LineWidth',2)
grid on
xlabel('\omega (rad/s)')
ylabel('Relative Flywheel Angle Amplitude (deg)')
legend('1-Pitch Flywheel 1','2-Pitch Flywheel 1', ...
    '2-Pitch Flywheel 2')
set(gca,'FontSize',12,'LineWidth',1.5)

%% Plot Phase Comparison

figure(2)

% Swing angles
subplot(2,1,1)
plot(w,ph_1(1,:),'b-','LineWidth',2)
hold on
plot(w,ph_2(1,:),'r--','LineWidth',2)
grid on
ylabel('Swing Phase (deg)')
title('Frequency Response Phase')
legend('1-Pitch Swing','2-Pitch Swing')
set(gca,'FontSize',12,'LineWidth',1.5)

% Relative flywheel angles
subplot(2,1,2)
plot(w,ph_1(2,:),'b-','LineWidth',2)
hold on
plot(w,ph_2(2,:),'r--','LineWidth',2)
plot(w,ph_2(3,:),'Color',[0.2 0.6 0.2], ...
    'LineStyle',':','LineWidth',2)
grid on
xlabel('\omega (rad/s)')
ylabel('Relative Flywheel Phase (deg)')
legend('1-Pitch Flywheel 1','2-Pitch Flywheel 1', ...
    '2-Pitch Flywheel 2')
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

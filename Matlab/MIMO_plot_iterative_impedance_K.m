%% Transfer Function Approach: sweep J2/J3 asymmetry
clear; close all; clc;
format short;

% H11 = Excitation torque to swing velocity
% H12 = Gearbox torque 2 to swing velocity
% H13 = Gearbox torque 3 to swing velocity
% H21 = Excitation torque to flywheel 2 velocity
% H22 = Gearbox torque 2 to flywheel 2 velocity
% H23 = Gearbox torque 3 to flywheel 2 velocity
% H31 = Excitation torque to flywheel 3 velocity
% H32 = Gearbox torque 2 to flywheel 3 velocity
% H33 = Gearbox torque 3 to flywheel 3 velocity
% H12 = H21, H13 = H31, H23 = H32 -> six unique plots

%% Base parameters (struct, so any sweep just overrides one field)
p.J1 = 22.45;   % swing inertia (kg*m^2)
p.K1 = 517;     % swing stiffness (Nm/rad)
p.C1 = 1.65;    % swing friction (Nm*s/rad)

J_flywheel    = 4.19;    % Lee et al. 2024
K_magspring   = 40.75;
C_sysfriction = 0.85;

p.J2 = J_flywheel;  p.K2 = K_magspring;  p.C2 = C_sysfriction;   % flywheel 2
p.J3 = J_flywheel;  p.K3 = K_magspring;  p.C3 = C_sysfriction;   % flywheel 3

p.Na = 3/1;   % Gearbox transmission ratio
p.Cgen = 0;    % Generator damping (Nm*s/rad)
p.Cpto = p.Na^2*p.Cgen;

p.Nb = 28/40;       % belt ratio

f = logspace(log10(0.05), log10(20), 3000);   % Hz

%% Sweep: K2 and K3 move apart, K_diff 
n  = 5;
K3 = linspace(K_magspring, K_magspring + 30, n);
K2 = linspace(K_magspring, K_magspring - 30, n);
K_diff = K3 - K2;

P = repmat(p, 1, n);
for k = 1:n
    P(k).K2 = K2(k);
    P(k).K3 = K3(k);
end

H = cell(1, n);
for k = 1:n
    H{k} = buildMIMO(P(k));
end

colors = summer(n); 
%% Label data
L.diff = K_diff;   fmt.diff = 'K_{diff} = %.2f';
L.K2   = K2;       fmt.K2   = 'K_2 = %.2f';
L.K3   = K3;       fmt.K3   = 'K_3 = %.2f';

lbl = makeLabels({'diff','K2','K3'}, L, fmt);   % same legend for every plot

%% Plots
figure(1);  plotSweep(H, 1, 2, f, colors, lbl, ...
    'H12 = H21: Gearbox torque 2 to swing velocity = Excitation torque to flywheel 2 velocity');

figure(2);  plotSweep(H, 1, 3, f, colors, lbl, ...
    'H13 = H31: Gearbox torque 3 to swing velocity = Excitation torque to flywheel 3 velocity');

figure(3);  plotSweep(H, 2, 3, f, colors, lbl, ...
    'H23 = H32: Gearbox torque 3 to flywheel 2 velocity = Gearbox torque 2 to flywheel 3 velocity');

figure(4);  plotSweep(H, 1, 1, f, colors, lbl, ...
    'H11: Excitation torque to swing velocity');

figure(5);  plotSweep(H, 2, 2, f, colors, lbl, ...
    'H22: Gearbox torque 2 to flywheel 2 velocity');

figure(6);  plotSweep(H, 3, 3, f, colors, lbl, ...
    'H33: Gearbox torque 3 to flywheel 3 velocity');

%% Save all figures (file name = figure title)
saveDir = 'C:\Users\ianhe\Projects\pioneerWEC\Figures\mimo_vary_K';   % <-- your folder path

figs = findall(groot, 'Type', 'figure');
for k = 1:numel(figs)
    tl = findobj(figs(k), 'Type', 'tiledlayout');
    ttl = tl.Title.String;
    fname = regexprep(ttl, '[\\/:*?"<>|\s]+', '_');      % strip characters illegal in filenames
    exportgraphics(figs(k), fullfile(saveDir, [fname '.png']), 'Resolution', 300);
end


%% Local functions
function H_mimo = buildMIMO(p)
    s  = tf('s');
    z1 = p.J1*s + p.C1 + p.K1/s;
    z2 = p.J2*s;
    z3 = p.J3*s;
    z4 = p.C2 + (p.Nb^2 * p.K2)/s;
    z5 = p.C3 + (p.Nb^2 * p.K3)/s;

    denom = z2*z4*z5 + z3*z4*z5 + z1*(z2+z4)*(z3+z5) + z2*z3*(z4+z5);

    H11 = ((z2+z4)*(z3+z5)) / denom;
    H12 = (-z2*(z3+z5)) / denom;
    H13 = (-z3*(z2+z4)) / denom;
    H21 = H12;
    H22 = (z3*z5 + z1*(z3+z5) + z2*(z3+z5)) / denom;
    H23 = (z2*z3) / denom;
    H31 = H13;
    H32 = H23;
    H33 = (z3*z4 + z1*(z2+z4) + z2*(z3+z4)) / denom;

    H_mimo = [H11 H12 H13; H21 H22 H23; H31 H32 H33];
end

function [magdB, phdeg] = bodeHz(sys, f)
    resp  = squeeze(freqresp(sys, 2*pi*f));
    magdB = 20*log10(abs(resp));
    phdeg = unwrap(angle(resp)) * 180/pi;
end

function lbl = makeLabels(names, L, fmt)
    % names: cell array of fields to include, e.g. {'diff','J2'}
    n   = numel(L.(names{1}));
    lbl = cell(1, n);
    for k = 1:n
        parts = cell(1, numel(names));
        for q = 1:numel(names)
            parts{q} = sprintf(fmt.(names{q}), L.(names{q})(k));
        end
        lbl{k} = strjoin(parts, ', ');
    end
end

function plotSweep(Hset, i, j, f, colors, labels, ttl)
    % Overlay element (i,j) of every case in Hset on one magnitude/phase figure
    tl  = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ttl, 'FontWeight', 'bold', 'Interpreter', 'none');
    axM = nexttile; hold(axM, 'on'); grid(axM, 'on');
    axP = nexttile; hold(axP, 'on'); grid(axP, 'on');
    for k = 1:numel(Hset)
        [magdB, phdeg] = bodeHz(Hset{k}(i,j), f);
        plot(axM, f, magdB, 'Color', colors(k,:), 'LineWidth', 1.4, 'DisplayName', labels{k});
        plot(axP, f, phdeg, 'Color', colors(k,:), 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end
    set([axM axP], 'XLim', [0.1 1]);
    set(axM, 'YLim', [-80 20]);
    set(axP, 'YLim', [-360 360], 'YTick', -360:180:360);
    ylabel(axM, 'Magnitude (dB)');
    ylabel(axP, 'Phase (deg)'); xlabel(axP, 'Frequency (Hz)');

    lg = legend(axM);
    lg.Layout.Tile = 'east';   % outside the axes, right side, vertically centered
end
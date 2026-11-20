clear; close; clc; close all;
format loose;
format short;
%% 2 Pitch-Resonator on Swing-Set: Define variables, matrices, and input torques
str = '2-pitch resonator, set-up';
fprintf('%s\n', strrep(pad([' ' str ' '], 80, 'both', '-'), str, ['<strong>' str '</strong>']));

%  Define Variables
Js = 22.45; % kg m^2, swing moment of inertia
Ks = 517; % Nm/rad, swing restoring stiffness
Cs = 1.65; % Nm*s/rad, swing friction
J1 = 4.19; % kg m^2, flywheel moment of inertia
K1 = 40.75; % Nm/rad, magnetic spring stiffness
C1 = 0.85; %Nm*s/rad, Effective system damping (friction)
J2=J1; 
K2=K1;
C2=C1;
Na = 3/1 ; % Gearbox transmission ratio between flywheel & generator shaft
Nb = 28/40; % Belt drive transmission ratio between magnetic spring and flywheel

% 65-Nm amplitude produces maximum swing-set motion +/- 10 deg
texc_amp = 65;      % N*m, amplitude of 65*sin(wt)

% t1 = theta_rel' * C_pto ; % Generator/motor acts as "pure damper"
% Cgen: 0.07 to 0.29 Nms/rad, corresponds to 0.64-2.58 Nms/rad gearbox
Cgen = 0.20; % Generator/motor damping coefficient (N*m*s/rad)
Cpto1 = Na^2 * Cgen;
Cpto2 = Cpto1 ; 

% Define inertia, damping, stiffness, force matrices 
J=[Js+J1+J2, J1, J2;
      J1,    J1, 0;
      J2,     0, J2];

C=[Cs,    0,     0;
    0, C1+Cpto1, 0;
    0,    0,   C2+Cpto2];

K=[Ks,    0,     0;
    0, K1*Nb^2,  0
    0,    0,   K2*Nb^2];

F = [texc_amp; 0; 0];         % unit direction: excitation on the swing DOF only

disp('J matrix is')
disp(J)
disp('K matrix is ')
disp(K)
disp('C matrix is ')
disp(C)
disp('F matrix is ')
disp(F)


%% Undamped/unforced evals, evecs, wn. 
str = 'Undamped/unforced. Use undamped evec to check if C proportional';
fprintf('%s\n', strrep(pad([' ' str ' '], 80, 'both', '-'), str, ['<strong>' str '</strong>']));

[evec,eval]=eig(K,J);

disp('Undampened Evals')
evals=diag(eval)

disp('Undamped Natural Frequencies')
wns=sqrt(evals)

disp('Undamped, Mass-normalized eigenvectors (MATLAB automatically mass normalizes)')
evec

disp('C matrix proportional? Ut*C*U diagonal?')
evec'*C*evec
disp('Non-proportional damping!')

disp('                           ')
disp('J matrix check: Ut*J*U = I? Mass-orthonormality')
disp(round(evec'*J*evec,10))

%% Convert to state-space (ss) A (state) and B (input). left/right evec, eval, wn. 
str = 'State-space x* = Ax + BQ. Left/right evec, eval, wn';
fprintf('%s\n', strrep(pad([' ' str ' '], 80, 'both', '-'), str, ['<strong>' str '</strong>']));

% Calc A (state) and B (input matrices)
Ass=[zeros(size(J)),eye(size(J));
    -inv(J)*K,-inv(J)*C];
Bss=[zeros(size(F));inv(J)*F];

disp('Ass matrix:')
disp(round(Ass,10))
disp('Bss matrix:')
disp(Bss)

disp('Use damp function to quickly find zetas')
damp(Ass)

[revecs,evals,levecs]=eig(Ass);

disp('Eigenvalues')
evalsvec=diag(evals)

disp('Use s=-zeta*wn+wn*sqt(1-zeta^2) to find wns and zetas')
disp('Natural Frequencies')
wnsvec=abs(evalsvec)

disp('Damping Ratio')
zetavec=real(-evalsvec)./wnsvec

disp('Right evecs before sorting')
revecs
disp('Left evecs before sorting')
levecs

%% Sort left and right evecs and normalize
str = 'Sort left/right evecs and normalize';
fprintf('%s\n', strrep(pad([' ' str ' '], 80, 'both', '-'), str, ['<strong>' str '</strong>']));

% Sort eigenvalues by |lambda| (= wn) ascending
% identical (conjugate pairs) magnitudes sorted by phase 
% (- imaginary 1, + imaginary 2)
[evalsvecsort,ind]=sort(evalsvec);

% apply same ind to left/right evecs so columns stay matched to modes.
revalsort=evalsvec(ind);
revecssort=revecs(:,ind);
levecssort=levecs(:,ind);

disp('Right evecs after sorting:')
revecssort
disp('Left evecs after sorting:')
levecssort

disp('Biorthogonality check after sorting YT * X = I?')
levecssort'*revecssort

% Scale each left/right pair by sqrt(YT * X) so YT*X = 1
for ii=1:length(revalsort),
    aind(ii)=sqrt(levecssort(:,ii)'*revecssort(:,ii));
    levecssortnorm(:,ii)=levecssort(:,ii)/aind(ii)';
    revecssortnorm(:,ii)=revecssort(:,ii)/aind(ii);
end

% disp('Use s=-zeta*wn+wn*sqt(1-zeta^2) to find wns and zetas')
wnsvec=abs(revalsort) ;
zetavec=real(-revalsort)./wnsvec ;

disp('Biorthogonality check after normalizing YT * X = I?')
disp('Show yt*x=1')
round(levecssortnorm'*revecssortnorm,1)
disp('Show yt*A*x=lam')
levecssortnorm'*Ass*revecssortnorm

disp('Right evecs after normalization/sorting')
revecssortnorm
disp('Left evecs after normalization/sorting')
levecssortnorm

%% Convert to modal coordinates q'(t) = Amodal*q + gamma
str = 'Convert to modal coordinates';
fprintf('%s\n', strrep(pad([' ' str ' '], 80, 'both', '-'), str, ['<strong>' str '</strong>']));

disp('A modal matrix (Am)')
Adecoupled=levecssortnorm'*Ass*revecssortnorm

disp('qdot in qdot=Am*q+gamma*u')
qdot=diag(levecssortnorm'*revecssortnorm)

disp('gamma in qdot=Am*q+gamma*u');
gamma=levecssortnorm'*Bss

%% Steady-state response to Q(t) = 65*sin(om*t) nt
% Complex: Q0*exp(j*om*t): real part -> cos(om*t) response, imag part -> sin(om*t) response
% Physical Space: x* = Ass*x + Bss*Q(t), Q(t)=sin(wt), Bss contains 65 amp
% Modal: 1*qdot - Am*q = gamma*Q(t), qdot = Am*q + gamma*Q(t), gamma = Y'B'
% steady-state response: qi = A1 sin(wt) + A2 cos(wt)
% A1 = -lam_i*gamma_i/(lam_i^2+w^2),  A2 = -w*gamma_i/(lam_i^2+w^2)

t = linspace(0,2,1000); % 1000 linearly spaced points from 0 to 2 sec

wn_pto_lee24_hz = 0.347 ;                         % PTO natural freq (from Lee 2024), Hz
wn_pto_lee24_rads = wn_pto_lee24_hz * 2 * pi ;    % PTO natural freq, rad/s
wn_swing_lee24_hz = 0.77 ;                        % Swing natural freq (from Lee 2024), Hz
wn_swing_lee24_rads = wn_swing_lee24_hz * 2 * pi ;% Swing natural freq, rad/s

w1 = wn_swing_lee24_rads ; % fixed w for modal/physical response, rad/s

for ii=1:size(Ass,1)
    lam=Adecoupled(ii,ii);              % modal eigenvalue lambda_i
    G=gamma(ii);                        % gamma = Y'*Bss already includes the 65 N*m
    A1=(-lam*G)/(lam^2+w1^2);           % sin coefficient
    A2=(-w1*G)/(lam^2+w1^2);            % cos coefficient
    qss(ii,:)=A1*sin(w1*t)+A2*cos(w1*t);
end

% Convert modal coordinates back to physical states
zss=revecssortnorm*qss;
    % zss(1,:) = theta_s
    % zss(2,:) = theta_rel1
    % zss(3,:) = theta_rel2
    % zss(4,:) = theta_s'
    % zss(5,:) = theta_rel1'
    % zss(6,:) = theta_rel2'

%% Plot Modal response at single specified frequency w
% Color definitions
q1_real = [0.153 0.702 0.145];   % green
q1_imag = [0.302 0.471 1.000];   % blue
q2_real = [0.922 0.851 0.000];   % gold
q2_imag = [0.902 0.000 0.000];   % red
q3_real = [0.494 0.184 0.556];   % purple
q3_imag = [0.000 0.750 0.750];   % teal

% Modal coordinates q1, q2, q3
figure(11)
plot(t, real(qss(1,:)), 'Color', q1_real, 'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, real(qss(2,:)), 'Color', q2_real, 'LineWidth', 2, 'LineStyle', '--')
plot(t, real(qss(3,:)), 'Color', q3_real, 'LineWidth', 2, 'LineStyle', '-.')
plot(t, imag(qss(1,:)), 'Color', q1_imag, 'LineWidth', 2, 'LineStyle', ':')
plot(t, imag(qss(2,:)), 'Color', q2_imag, 'LineWidth', 2, 'LineStyle', '-.')
plot(t, imag(qss(3,:)), 'Color', q3_imag, 'LineWidth', 2, 'LineStyle', '--')
hold off
title(sprintf('Modal Coordinates q1, q2, q3 at w = %.3f rad/s',w1),'FontSize',20)
xlabel('Time (s)','FontSize',20)
ylabel('Modal amplitude','FontSize',20)
legend('q1-real','q2-real','q3-real','q1-imag','q2-imag','q3-imag')
grid on
set(gca,'FontSize',12,'LineWidth',2)

% Modal coordinates q4, q5, q6 (reuse colors of the matching mode pair)
figure(12)
plot(t, real(qss(4,:)), 'Color', q1_real, 'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, real(qss(5,:)), 'Color', q2_real, 'LineWidth', 2, 'LineStyle', '--')
plot(t, real(qss(6,:)), 'Color', q3_real, 'LineWidth', 2, 'LineStyle', '-.')
plot(t, imag(qss(4,:)), 'Color', q1_imag, 'LineWidth', 2, 'LineStyle', ':')
plot(t, imag(qss(5,:)), 'Color', q2_imag, 'LineWidth', 2, 'LineStyle', '-.')
plot(t, imag(qss(6,:)), 'Color', q3_imag, 'LineWidth', 2, 'LineStyle', '--')
hold off
title(sprintf('Modal Coordinates q4, q5, q6 at w = %.3f rad/s',w1),'FontSize',20)
xlabel('Time (s)','FontSize',20)
ylabel('Modal amplitude','FontSize',20)
legend('q4-real','q5-real','q6-real','q4-imag','q5-imag','q6-imag')
grid on
set(gca,'FontSize',12,'LineWidth',2)

%% Colors / styles for swing and flywheels
c_swing = [0.000 0.447 0.741];   % blue
c_fw1   = [0.788 0.161 0.00];    % red
c_fw2   = [1 0.871 0.282];       % yellow

%% Plot physical response at single specified frequency w
% Physical angles (deg)
figure(13)
plot(t, rad2deg(real(zss(1,:))), 'Color', c_swing, 'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, rad2deg(real(zss(2,:))), 'Color', c_fw1, 'LineWidth', 2, 'LineStyle', '-')
plot(t, rad2deg(real(zss(3,:))), 'Color', c_fw2, 'LineWidth', 2, 'LineStyle', '--')
hold off
title(sprintf('Angles at w = %.3f rad/s',w1),'FontSize',20)
xlabel('Time (s)','FontSize',20), ylabel('Angle (deg)','FontSize',20)
legend('swing \theta_s','flywheel 1 rel. \theta_{rel,1}','flywheel 2 rel. \theta_{rel,2}')
grid on
set(gca,'FontSize',12,'LineWidth',2)

% Physical angular velocities (rad/s)
figure(14)
plot(t, real(zss(4,:)), 'Color', c_swing, 'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, real(zss(5,:)), 'Color', c_fw1, 'LineWidth', 2, 'LineStyle', '-')
plot(t, real(zss(6,:)), 'Color', c_fw2, 'LineWidth', 2, 'LineStyle', '--')
hold off
title(sprintf('Angular Velocity at w = %.3f rad/s',w1),'FontSize',20)
xlabel('Time (s)','FontSize',20), ylabel('Velocity (rad/s)','FontSize',20)
legend('swing \omega_s','flywheel 1 rel. \omega_{rel,1}','flywheel 2 rel. \omega_{rel,2}')
grid on
set(gca,'FontSize',12,'LineWidth',2)

%% Plot frequency amplitude/phase response
f_hz = linspace(0.1,1,4000) ;   % Exc frequency range, Hz
w_rads = f_hz * 2 * pi ;        % Angular exc frequency, rad/s
X = zeros(3, numel(w_rads));

for k = 1:numel(w_rads)
    X(:,k) = (K - w_rads(k)^2*J + 1j*w_rads(k)*C) \ F;
end

amp_s    = rad2deg(abs(X(1,:)));
amp_rel1 = rad2deg(abs(X(2,:)));
amp_rel2 = rad2deg(abs(X(3,:)));
ph_s     = rad2deg(angle(X(1,:)));
ph_rel1  = rad2deg(angle(X(2,:)));
ph_rel2  = rad2deg(angle(X(3,:)));

figure(15)
subplot(2,1,1)
semilogy(f_hz, amp_s,    'Color', c_swing, 'LineWidth', 2, 'LineStyle', '-')
hold on
semilogy(f_hz, amp_rel1, 'Color', c_fw1,   'LineWidth', 2, 'LineStyle', '-')
semilogy(f_hz, amp_rel2, 'Color', c_fw2,   'LineWidth', 2, 'LineStyle', '--')
yline(10,'--k','10 deg limit')
xline(wns/(2*pi),':','Color',[0.4 0.4 0.4])   % natural frequencies in Hz
hold off
xlim([f_hz(1) f_hz(end)])
grid on
ylabel('Amplitude (deg)')
legend('swing','flywheel 1 rel.','flywheel 2 rel.','10 deg limit')
title('2-Pitch Resonator Frequency Response')
set(gca,'FontSize',12,'LineWidth',2)

subplot(2,1,2)
plot(f_hz, ph_s,    'Color', c_swing, 'LineWidth', 2, 'LineStyle', '-')
hold on
plot(f_hz, ph_rel1, 'Color', c_fw1,   'LineWidth', 2, 'LineStyle', '-')
plot(f_hz, ph_rel2, 'Color', c_fw2,   'LineWidth', 2, 'LineStyle', '--')
hold off
xlim([f_hz(1) f_hz(end)])
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
legend('swing','flywheel 1 rel.','flywheel 2 rel.')
set(gca,'FontSize',12,'LineWidth',2)

%% Frequency Response: Cgen Sweep, Plot response amplitude
Cgen_range = 0.07:0.01:0.29;
amp_rel1_sweep = zeros(length(Cgen_range),length(w_rads));

for i=1:length(Cgen_range)
    Cpto = Na^2*Cgen_range(i);
    C_sweep = diag([Cs,C1+Cpto,C2+Cpto]);
    for k=1:length(w_rads)
        X = (K-w_rads(k)^2*J+1j*w_rads(k)*C_sweep)\F;
        amp_rel1_sweep(i,k) = rad2deg(abs(X(2)));
    end
end

figure(17)
hold on
colors = parula(length(Cgen_range));
for i=1:length(Cgen_range)
    plot(f_hz,amp_rel1_sweep(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
hold off
colormap(parula)
cb = colorbar;
clim([min(Cgen_range) max(Cgen_range)])
ylabel(cb,'C_{gen} (Nms/rad)')
xlim([f_hz(1) f_hz(end)])
grid on
xlabel('Frequency (Hz)')
ylabel('Relative Flywheel 1 Amplitude (deg)')
title('2-Pitch Resonator: Varied Generator Damping')
set(gca,'FontSize',12,'LineWidth',2)
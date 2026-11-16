clear; close; clc; close all;
format loose;
format short;
%% 1 Pitch-Resonator on Swing-Set: Define variables, matrices, and input torques
str = '1-pitch resonator, set-up';
fprintf('%s\n', strrep(pad([' ' str ' '], 80, 'both', '-'), str, ['<strong>' str '</strong>']));

%  Define Variables
% DOF 1 = swing, DOF 2 = flywheel (relative to swing)
J1 = 22.45; % kg m^2, swing moment of inertia
K1 = 517; % Nm/rad, swing restoring stiffness
C1 = 1.65; % Nm*s/rad, swing friction
J2 = 4.19; % kg m^2, flywheel moment of inertia
K2 = 40.75; % Nm/rad, magnetic spring stiffness
C2 = 0.85; %Nm*s/rad, Effective system damping (friction)
Na = 3/1 ; % Gearbox transmission ratio between flywheel & generator shaft
Nb = 28/40; % Belt drive transmission ratio between magnetic spring and flywheel

% 65 Nm amplitude produces maximum swing-set motion +/- 10 deg
texc_amp = 65;      % N*m, amplitude of 65*sin(wt)

% t1 = theta_rel' * C_pto ; % Generator/motor acts as "pure damper"
% Cgen: 0.07 to 0.29 Nms/rad, corresponds to 0.64-2.58 Nms/rad gearbox
Cgen = 0 ; % Generator/motor damping coefficient (N*m*s/rad)
Cpto = Na^2 * Cgen ;

% Define inertia, damping, stiffness, force matrices 
J=[J1+J2, J2;
    J2, J2];

C=[C1, 0;
    0, C2+Cpto];

K=[K1, 0;
    0, K2*Nb^2];

F = [texc_amp; 0];         % unit direction: excitation on the swing DOF only

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

t = linspace(0,2,1000); % 1000 linearly space points from 0 set to 10 sec

wn_pto_lee24_hz = 0.347 ; % PTO natural freq (from Lee 2024), Hz
wn_pto_lee24_rads = wn_pto_lee24_hz * 2 * pi ; % PTO natural freq (from Lee 2024), m/s
wn_swing_lee24_hz = 0.77 ; % Swing natural freq (from Lee 2024), Hz
wn_swing_lee24_rads = wn_swing_lee24_hz * 2 * pi ; % Swing natural freq (from Lee 2024), Hz

w1 = wn_swing_lee24_rads ; % fixed w for modal/physical response, m/s

for ii=1:size(Ass,1)
    lam=Adecoupled(ii,ii);              % modal eigenvalue lambda_i
    G=gamma(ii);                        % gamma = Y'*Bss already includes the 65 N*m
    A1=(-lam*G)/(lam^2+w1^2);           % sin coefficient
    A2=(-w1*G)/(lam^2+w1^2);            % cos coefficient
    qss(ii,:)=A1*sin(w1*t)+A2*cos(w1*t);
end

% Convert modal coordinates back to physical states
zss=revecssortnorm*qss;
    % zss(1,:) = theta_1   (swing)
    % zss(2,:) = theta_2   (flywheel, relative to swing)
    % zss(3,:) = theta_1'
    % zss(4,:) = theta_2'


%% Plot Modal response at single specified frequency w

% Modal coordinates
figure(1)
plot(t, real(qss(1,:)),  'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, real(qss(2,:)),  'LineWidth', 2, 'LineStyle', '--')
plot(t, imag(qss(1,:)),  'LineWidth', 2, 'LineStyle', ':')
plot(t, imag(qss(2,:)),  'LineWidth', 2, 'LineStyle', '-.')
hold off
title('Modal Coordinates q1, q2','FontSize',20)
xlabel('Time (s)','FontSize',20)
ylabel('Modal amplitude','FontSize',20)
legend('q1-real','q2-real','q1-imag','q2-imag')
grid on
set(gca,'FontSize',12,'LineWidth',2)

figure(2)
plot(t, real(qss(3,:)),  'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, real(qss(4,:)),  'LineWidth', 2, 'LineStyle', '--')
plot(t, imag(qss(3,:)),  'LineWidth', 2, 'LineStyle', ':')
plot(t, imag(qss(4,:)),  'LineWidth', 2, 'LineStyle', '-.')
hold off
title('Modal Coordinates q3, q4','FontSize',20)
xlabel('Time (s)','FontSize',20)
ylabel('Modal amplitude','FontSize',20)
legend('q3-real','q4-real','q3-imag','q4-imag')
grid on
set(gca,'FontSize',12,'LineWidth',2)


%% Plot physical response at single specified frequency w
% Physical angles (deg)
figure(3)
plot(t,rad2deg(real(zss(1:2,:))),'LineWidth',2)
title(sprintf('Angles at w = %.3f rad/s',w1),'FontSize',20)
xlabel('Time (s)','FontSize',20), ylabel('Angle (deg)','FontSize',20)
legend('swing \theta_1','flywheel rel. \theta_2'), grid on
set(gca,'FontSize',12,'LineWidth',2)

% Physical angular velocities (rad/s)
figure(4)
plot(t,real(zss(3:4,:)),'LineWidth',2)
title('Angular Velocities','FontSize',20)
xlabel('Time (s)','FontSize',20), ylabel('Velocity (rad/s)','FontSize',20)
legend('swing \omega_1','flywheel rel. \omega_2'), grid on
set(gca,'FontSize',12,'LineWidth',2)

%% Plot frequency amplitude/phase response
f_hz = linspace(0.1,1,4000) ; % Exc frequency range hz, based on lee et. al 2024 fig 8
w_rads = f_hz * 2 * pi ; % Angular exc frequency range, ms 
X = zeros(2, numel(w_rads));

% EOM: J*theta'' + C*theta'+K*theta = F*sin(w,t), theta=(theta_1,theta_2)
for k = 1:numel(w_rads)
    X(:,k) = (K - w_rads(k)^2*J + 1j*w_rads(k)*C) \ F;   % [theta_1; theta_2] amplitudes
end

amp_1 = rad2deg(abs(X(1,:)));    % swing amplitude, deg
amp_2 = rad2deg(abs(X(2,:)));    % flywheel relative amplitude, deg
ph_1  = rad2deg(angle(X(1,:)));  % swing phase 
ph_2  = rad2deg(angle(X(2,:)));  % flywheel-relative phase 

figure(5)

subplot(2,1,1)
plot(f_hz, amp_1, f_hz, amp_2, 'LineWidth', 2)
xlim([f_hz(1) f_hz(end)])
grid on
ylabel('Amplitude (deg)')
legend('swing','flywheel rel.')
title('Frequency response')

subplot(2,1,2)
plot(f_hz, ph_1, f_hz, ph_2, 'LineWidth', 2)
xlim([f_hz(1) f_hz(end)])
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
legend('swing','flywheel rel.')

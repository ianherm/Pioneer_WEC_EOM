clear; close; clc;
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

t = linspace(0,10,1000); % 1000 linearly space points from 0 set to 10 sec
w = 0.25:0.05:0.45 ; % Excitation frequency range

% Complex: Q0*exp(j*om*t): real part -> cos(om*t) response, imag part -> sin(om*t) response
% Physical Space: x* = Ass*x + Bss*Q(t), Q(t)=sin(wt), Bss contains 65 amp
% Modal: 1*qdot - Am*q = gamma*Q(t), qdot = Am*q + gamma*Q(t), gamma = Y'B'
% steady-state response: qi = A1 sin(wt) + A2 cos(wt)
% A1 = -lam_i*gamma_i/(lam_i^2+w^2),  A2 = -w*gamma_i/(lam_i^2+w^2)

w1 = 4.9;                              % fixed excitation frequency 
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
    % zss(2,:) = theta_rel
    % zss(3,:) = theta_s'
    % zss(4,:) = theta_rel'


    %% Plot

    % Modal coordinates
    figure(11)
    plot(t,real(qss(1,:)),'b-','LineWidth',2)
    hold on
    plot(t,real(qss(2,:)),'r--','LineWidth',2)
    plot(t,real(qss(3,:)),'y:','LineWidth',2)
    title('Modal Coordinates','FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Modal Amplitude','FontSize',20)
    legend('q1','q2','q3')
    grid on
    set(gca,'FontSize',12,'LineWidth',2)

    figure(12)
    plot(t,real(qss(4,:)),'b-','LineWidth',2)
    hold on
    plot(t,real(qss(5,:)),'r--','LineWidth',2)
    plot(t,real(qss(6,:)),'y:','LineWidth',2)
    title('Modal Velocities','FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Modal Amplitude','FontSize',20)
    legend('q4','q5','q6')
    grid on
    set(gca,'FontSize',12,'LineWidth',2)

    % Physical angles
    figure(13)
    plot(t,rad2deg(real(zss(1,:))),'b-','LineWidth',2)
    hold on
    plot(t,rad2deg(real(zss(2,:))),'r--','LineWidth',2)
    plot(t,rad2deg(real(zss(3,:))),'y:','LineWidth',2)
    title(sprintf('Angles at w = %.3f rad/s',w1),'FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Angle (deg)','FontSize',20)
    legend('Swing \theta_s','Flywheel 1 \theta_{rel,1}', ...
        'Flywheel 2 \theta_{rel,2}')
    grid on
    set(gca,'FontSize',12,'LineWidth',2)

    % Physical angular velocities
    figure(14)
    plot(t,real(zss(4,:)),'b-','LineWidth',2)
    hold on
    plot(t,real(zss(5,:)),'r--','LineWidth',2)
    plot(t,real(zss(6,:)),'y:','LineWidth',2)
    title('Angular Velocities','FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Velocity (rad/s)','FontSize',20)
    legend('Swing \omega_s','Flywheel 1 \omega_{rel,1}', ...
        'Flywheel 2 \omega_{rel,2}')
    grid on
    set(gca,'FontSize',12,'LineWidth',2)

    %% Frequency Response

    w = linspace(0.1,8,4000);
    X = zeros(3,numel(w));

    for k=1:numel(w)
        X(:,k)=(K-w(k)^2*J+1j*w(k)*C)\F;
    end

    amp_s=rad2deg(abs(X(1,:)));
    amp_rel1=rad2deg(abs(X(2,:)));
    amp_rel2=rad2deg(abs(X(3,:)));

    ph_s=rad2deg(angle(X(1,:)));
    ph_rel1=rad2deg(angle(X(2,:)));
    ph_rel2=rad2deg(angle(X(3,:)));

    figure(15)
    semilogy(w,amp_s,'b-','LineWidth',2)
    hold on
    semilogy(w,amp_rel1,'r--','LineWidth',2)
    semilogy(w,amp_rel2,'y:','LineWidth',2)
    yline(10,'--k','10 deg limit')
    xline(wns,':','Color',[0.4 0.4 0.4])
    grid on
    xlabel('\omega (rad/s)')
    ylabel('Amplitude (deg)')
    title('2-Pitch Resonator Frequency Response')
    legend('Swing','Flywheel 1 relative','Flywheel 2 relative', ...
        '10 deg limit')
    set(gca,'FontSize',12,'LineWidth',2)

    figure(16)
    plot(w,ph_s,'b-','LineWidth',2)
    hold on
    plot(w,ph_rel1,'r--','LineWidth',2)
    plot(w,ph_rel2,'y:','LineWidth',2)
    grid on
    xlabel('\omega (rad/s)')
    ylabel('Phase (deg)')
    title('2-Pitch Resonator Phase Response')
    legend('Swing','Flywheel 1 relative','Flywheel 2 relative')
    set(gca,'FontSize',12,'LineWidth',2)

    %% Frequency Response: Cgen Sweep, Plot response amplitude

    Cgen_range = 0.07:0.01:0.29;
    w = linspace(0.1,8,4000);

    amp_rel1_sweep = zeros(length(Cgen_range),length(w));

    for i=1:length(Cgen_range)

        Cpto = Na^2*Cgen_range(i);
        C_sweep = diag([Cs,C1+Cpto,C2+Cpto]);

        for k=1:length(w)
            X = (K-w(k)^2*J+1j*w(k)*C_sweep)\F;

            amp_rel1_sweep(i,k) = rad2deg(abs(X(2)));
        end

    end

    figure(17)
    hold on

    % Create sequential color gradient
    colors = parula(length(Cgen_range));

    for i=1:length(Cgen_range)
        plot(w,amp_rel1_sweep(i,:),'Color',colors(i,:),'LineWidth',1.5)
    end

    colormap(parula)
    cb = colorbar;
    clim([min(Cgen_range) max(Cgen_range)])
    ylabel(cb,'C_{gen} (Nms/rad)')

    grid on
    xlabel('\omega (rad/s)')
    ylabel('Relative Flywheel 1 Amplitude (deg)')
    title('2-Pitch Resonator: Varied Generator Damping')
    set(gca,'FontSize',12,'LineWidth',2)
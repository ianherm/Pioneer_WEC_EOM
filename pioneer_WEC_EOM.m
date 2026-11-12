clear; close; clc;
format loose;
format short;
%% 1 resonator

%  Define Variables
Js = 22.45; % kg m^2, swing moment of inertia
Ks = 517; % Nm/rad, swing restoring stiffness
Bs = 1.65; % Nm*s/rad, swing friction
J1 = 4.19; % kg m^2, flywheel moment of inertia
K1 = 40.75; % Nm/rad, magnetic spring stiffness
B1 = 0.85; %Nm*s/rad, Effective system damping (friction)
Na = 3/1 ; % Gearbox transmission ratio between flywheel & generator shaft
Nb = 28/40; % Belt drive transmission ratio between magnetic spring and flywheel

texc_amp = 65;      % N*m, amplitude of 65*sin(wt)

% t1 = theta_rel' * B_pto ; % Generator/motor acts as "pure damper"
% Bgen: 0.07 to 0.29 Nms/rad, corresponds to 0.64-2.58 Nms/rad gearbox
% B_gen = 0.07:0.29:0.01 ;
Bgen = 0.20; % Generator/motor damping coefficient (N*m*s/rad)
Bpto = Na^2 * Bgen;


% Define inertia, damping, stiffness, force matrices 
J=[Js+J1, J1;
    J1, J1];

B=[Bs, 0;
    0, B1+Bpto];

K=[Ks, 0;
    0, K1*Nb^2];

F = [texc_amp; 0];         % unit direction: excitation on the swing DOF only

disp('J matrix is ')
disp(J)
disp('K matrix is ')
disp(K)
disp('B matrix is ')
disp(B)
disp('F matrix is ')
disp(F)

[evec,eval]=eig(K,J);

disp('Undampened Evals')
evals=diag(eval)

disp('Undamped Natural Frequencies')
wns=sqrt(evals)

disp('Mass normalized eigenvectors (MATLAB automatically mass normalizes')
evec

disp('B matrix proportional? Ut*C*U')
evec'*B*evec
disp('Non-proportional damping!')

disp('                           ')
disp('J matrix check: Ut*J*U = I?')
disp(round(evec'*J*evec,10))


Ass=[zeros(size(J)),eye(size(J));
    -inv(J)*K,-inv(J)*B];
Bss=[zeros(size(F));inv(J)*F];

% 
disp('A matrix is ')
disp(round(Ass,10))
disp('B matrix is Fo=1')
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

[evalsvecsort,ind]=sort(evalsvec)
revalsort=evalsvec(ind);
revecssort=revecs(:,ind);
levecssort=levecs(:,ind);

disp('Right evecs after sorting')
revecssort
disp('Left evecs after sorting')
levecssort

levecssort'*revecssort

for ii=1:4,
    aind(ii)=sqrt(levecssort(:,ii)'*revecssort(:,ii));
    levecssortnorm(:,ii)=levecssort(:,ii)/aind(ii)';
    revecssortnorm(:,ii)=revecssort(:,ii)/aind(ii);
end

disp('Use s=-zeta*wn+wn*sqt(1-zeta^2) to find wns and zetas')
wnsvec=abs(revalsort)
zetavec=real(-revalsort)./wnsvec

disp('Eigenvalues of A after sorting')
evalsvecsort
disp('Right evecs after normalization/sorting')
revecssortnorm
disp('Left evecs after normalization/sorting')
levecssortnorm

levecssort'*revecssort
disp('Show yt*x=1')
((levecssortnorm'*revecssortnorm))
disp('Show yt*A*x=lam')
((levecssortnorm'*Ass*revecssortnorm))
% 
disp('A modal matrix')
Adecoupled=levecssortnorm'*Ass*revecssortnorm
disp('qdot in qdot=Am*q+gamma*u')
qdot=diag(levecssortnorm'*revecssortnorm)
disp('gamma in qdot=Am*q+gamma*u');
gamma=levecssortnorm'*Bss


tvec=linspace(0,10,1000);


disp('a*zdot+b*z=G')

om=wnsvec(1);
mu=0.2;
ru=0.25;

for ii=1:size(Ass,1)
    a=qdot(ii);
    b=Adecoupled(ii,ii);
    G=mu*ru*om^2*gamma(ii);
    A1=(a*G*om)/(b^2+a^2*om^2);
    A2=(b*G)/(b^2+a^2*om^2);
    qss(ii,:)=A1*sin(om*tvec)+A2*cos(om*tvec);
end

zss=revecssortnorm*qss;

figure(11)
plot(tvec,real(qss(1:2,:)),tvec,imag(qss(1:2,:)),'LineWidth',2)
title('Modal Displacements','FontSize',20);
xlabel('Time (s)','FontSize',20)
ylabel('Disp (m)','FontSize',20)
legend('q1-real','q2-real','q1-imag','q2-imag')
h=gca;
set(h,'FontSize',12,'LineWidth',2)
grid on


figure(12)
plot(tvec,real(qss(3:4,:)),tvec,imag(qss(3:4,:)),'LineWidth',2)
title('Modal Displacements','FontSize',20);
xlabel('Time (s)','FontSize',20)
ylabel('Disp (m)','FontSize',20)
legend('q3-real','q4-real','q3-imag','q4-imag')
h=gca;
set(h,'FontSize',12,'LineWidth',2)
grid on

figure(13)
plot(tvec,real(zss(1:4,:)),'LineWidth',2)
title('Displacements','FontSize',20);
xlabel('Time (s)','FontSize',20)
ylabel('Disp (m)','FontSize',20)
legend('y1','y2','y3','y4')
h=gca;
set(h,'FontSize',12,'LineWidth',2)
grid on
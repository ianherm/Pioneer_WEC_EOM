clear all;
format short e;

% 

m=2;
k=1000;
c=5;


M=[m,0;
    0,m];

K=[2*k,-k;
    -k,k];

C=[2*c,-c;
    -c,3*c];

F=[1;
    0];

disp('M matrix is ')
disp(M)
disp('K matrix is ')
disp(K)
disp('C matrix is ')
disp(C)
disp('F matrix is ')
disp(F)


[evec,eval]=eig(K,M);

disp('Undampened Evals')
evals=diag(eval)

disp('Undamped Natural Frequencies')
wns=sqrt(evals)

disp('Mass normalized eigenvectors')
evec

disp('C matrix not proportional Ut*C*U')
evec'*C*evec

Ass=[zeros(size(M)),eye(size(M));
    -K/M,-C/M];
Bss=[zeros(size(F));inv(M)*F];


disp('A matrix is ')
disp(Ass)
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

for ii=1:8,
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
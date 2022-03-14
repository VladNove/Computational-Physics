% Pendulul gravitational dublu
% Novetschi Vlad
clear; close all; clc;

%am implementat pendulul gravitational triplu, folosind ecuatile deduse din lagrangian.
%din pacate, posibil din cauza unei erori de calcul, sistemul este instabil numeric, rezultand in
%viteze si unghiuri extrem de mari dupa cateva iteratii, ce depasesc precizia matlab-ului.
%am incercat sa rezolvam aceasta problema mentinand unghiurile in intervalul 0 2pi 
%si limitand vitezele pendulului, insa sistemul nu mentine o energie constanta.
%am lucrat impreuna cu colegul meu Anghel Andrei pentru calculul si implementarea ecuatilor.

RT=0; % selector real time(1) / slow motion(0)

g=9.80665; % m/s^2; acceleratia gravitationala terestra

% Parametrii fizici ai sistemului mecanic:

L1=2.2; L2=1.3; L3 = 1.1 % m; lungimile tijelor

m1=1.2; m2=1.7; m3= 1.9 % kg; masele corpurilor

% Conditiile initiale - orice unghiuri intre -180 si +180!

theta10=100; theta20=70; theta30=20 % grade; unghiurile initiale

theta10=theta10*pi/180; theta20=theta20*pi/180; theta30=theta30*pi/180 % conversia in rad

OM10=10; OM20=-20; OM30=0; % grade/s; vitezele unghiulare initiale

OM10=OM10*pi/180; OM20=OM20*pi/180; OM30=OM30*pi/180% conversia in rad/s

% Definirea duratelor caracteristice:

omega1=sqrt(g/L1); omega2=sqrt(g/L2); omega3=sqrt(g/L3) % pulsatii proprii ale componentelor

T1=2*pi/omega1; T2=2*pi/omega2; T3=3*pi/omega3% perioade proprii ale componentelor

T=max(T1,T2);
T=max(T,T3); % "timp caracteristic" al miscarii pendulului dublu

ti=0; tf=5*T; N=10000; t=linspace(ti,tf,N); dt=t(2)-t(1); % timpul discret

% Prealocare si valori de start:

theta1=zeros(1,N); theta2=theta1; theta3=theta1; % prealocare unghiuri

OM1=zeros(1,N); OM2=OM1; OM3=OM1;% prealocare viteze unghiulare

T=zeros(1,N);
U=zeros(1,N);
H=zeros(1,N);

theta1(1)=theta10; theta2(1)=theta20; theta3(1)=theta30; % unghiuri de start pas 1

theta1(2)=theta10+OM10*dt; theta2(2)=theta20+OM20*dt; theta3(2)=theta30+OM30*dt; % unghiuri de start pas 2

OM1(1)=OM10; OM2(1)=OM20;  OM3(1)=OM30; % valori de start ale vitezelor unghiulare

% Notatii ajutatoare:

miu=1+m1/m2; % coeficient adimensional

r=L2/L1; % coeficient adimensional

%a11=miu; a22=r; % coeficienti diagonala principala (constanti)

%Coeficienti constanti de pe diagonala principala
a11 = L1^2*(m1 + m2 + m3);
a22 = L2^2*(m2+m3);
a33 = m3*L3^2;
%%%
tic;

b1=0;
b2=0;
b3=0;

for i=2:N-1 % ciclul recurentelor

    %aux=theta2(i)-theta1(i);
    %a21=cos(aux); a12=a21*r; % coeficienti diagonala secundara (variabili)

    %coeficienti variabili 
    a12 = (m2+m3)*L1*L2*cos(theta1(i)-theta2(i));
    a21 = (m3+m2)*L1*L2*cos(theta2(i)-theta1(i));
    a13 = m3*L1*L3*cos(theta1(i)-theta3(i));
    a31 = m3*L1*L3*cos(theta1(i)-theta3(i));
    a32 = m3*L2*L3*cos(theta2(i)-theta3(i));
    a23 = m3*L2*L3*cos(theta2(i)-theta3(i));
    %%
    % Pentru vitezele unghiulare curente folosim derivatele la stanga:

    OM1(i)=(theta1(i)-theta1(i-1))/dt; % viteza corpului 1 la pasul i

    OM2(i)=(theta2(i)-theta2(i-1))/dt; % viteza corpului 2 la pasul i
    
    OM3(i)=(theta3(i)-theta3(i-1))/dt; % viteza corpului 3 la pasul i
    
    U(i)=-g*((m1+m2+m3)*L1*cos(theta1(i))+(m2+m3)*L2*cos(theta2(i))+m3*L3*cos(theta3(i)));
    vabs2=OM2(i)^2+OM1(i)^2+2*OM1(i)*OM2(i)*cos(theta2(i)-theta1(i));
    unghivabs2=atan(OM2(i)*sin(theta2(i)-theta1(i))/(OM1(i)+OM2(i)*cos(theta2(i)-theta1(i))));
    T(i)=1/2*(m1*OM1(i)^2+m2*vabs2+m3*OM3(i)^2+m3*(OM3(i)^2+vabs2^2+2*vabs2*OM3(i)*cos(theta3(i)-theta2(i)+unghivabs2)));
    H(i)=U(i)+T(i);
    
    
     limitaviteza = 0.001;
    if (abs(OM1(i))>limitaviteza)
      OM1(i)=limitaviteza*sign(OM1(i));
    endif
    if (abs(OM2(i))>limitaviteza)
      OM2(i)=limitaviteza*sign(OM2(i));
    endif
    if (abs(OM3(i))>limitaviteza)
      OM3(i)=limitaviteza*sign(OM3(i));
    endif
   
    
    %b1=r*OM2(i)^2*sin(aux)-g/L1*miu*sin(theta1(i)); % termen "liber" 1

    %b2=-OM1(i)^2*sin(aux)-g/L1*sin(theta2(i)); % termen "liber" 2

    %termeni liberi
    b1 = g*L1*sin(theta1(i))*(m1+m2+m3) + m2*L1*L2*sin(theta1(i)-theta2(i))*OM1(i)*OM2(i) + m3*L1*L3*sin(theta1(i)-theta3(i)) + m3*L1*L2*sin(theta1(i)-theta2(i))*OM1(i)*OM2(i)  + (m2+m3)*L1*L2*(sin(theta2(i)-theta1(i))*(OM1(i)-OM2(i))*OM2(i)) + m3*L1*L3*(sin(theta3(i)-theta1(i))*(OM1(i)-OM3(i))*OM3(i));
    b1=-b1;
    
    b2 = g*L2*(m2*sin(theta2(i)) + m3*sin(theta2(i))) + OM1(i)*OM2(i)*L1*L2*sin(theta2(i)-theta1(i))*(m2 + m3)+ m3*L2*L3*sin(theta2(i)-theta3(i))*OM2(i)*OM3(i) + (m2 + m3)*L1*L2*sin(theta2(i)-theta1(i))*(OM1(i)-OM2(i))*OM1(i) + m3*L2*L3*sin(theta3(i) - theta2(i))*(OM2(i) - OM3(i))*OM3(i);
    b2=-b2;
    
    b3 = m3*g*L3*sin(theta3(i)) * m3*L2*L3 * sin(theta2(i) - theta3(i))*OM2(i)*OM3(i) - m3*L1*L3*sin(theta1(i) - theta3(i))*OM1(i)*OM3(i) + m3*L1*L3*sin(theta3(i) - theta1(i))*(OM1(i) - OM3(i))*OM1(i) + m3*L2*L3*sin(theta3(i) - theta2(i))*(OM2(i) - OM3(i))*OM2(i);
    b3=-b3;
    
    %%
    
    
    
    A=[a11,a12,a13;a21,a22,a23;a31,a32,a33]; B=[b1;b2;b3]; % matrice sistem si coloana termeni liberi

    E=A\B; % rezolvarea sistemului liniar in forma matriceala

    eps1=E(1); eps2=E(2); eps3=E(3);% acceleratiile unghiulare curente

    % Recurentele de ordinul II:

    theta1(i+1)=2*theta1(i)-theta1(i-1)+dt^2*eps1; % corp 1

    theta2(i+1)=2*theta2(i)-theta2(i-1)+dt^2*eps2; % corp 2
    
    theta3(i+1)=2*theta3(i)-theta3(i-1)+dt^2*eps3; % corp 3
    
    theta1(i+1) = mod (theta1(i+1)+pi, 2*pi)-pi;
    theta2(i+1) = mod (theta2(i+1)+pi, 2*pi)-pi;
    theta3(i+1) = mod (theta3(i+1)+pi, 2*pi)-pi;
   

endfor;

OM1(N)=(theta1(N)-theta1(N-1))/dt; % viteza unghiulara corp 1 la pasul N

OM2(N)=(theta2(N)-theta2(N-1))/dt; % viteza unghiulara corp 2 la pasul N

OM3(N)=(theta3(N)-theta3(N-1))/dt; % viteza unghiulara corp 3 la pasul N
toc; % afiseaza timpul de calcul al solutiei numerice

% Coordonate carteziene ale corpurilor:

x1=L1*sin(theta1); x2=x1+L2*sin(theta2); x3=x2+L3*sin(theta3);% coordonate orizontale

y1=-L1*cos(theta1); y2=y1-L2*cos(theta2); y3=y2-L3*cos(theta3);% coordonate verticale

% Energiile cinetica, potentiala, si respectiv totala:

%T=1/2*(m1*L1^2*OM1.^2+m2*(L1^2*OM1.^2+L2^2*OM2.^2+2*L1*L2*OM1.*OM2.*cos(theta2-theta1)));

%U=-g*((m1+m2)*L1*cos(theta1)+m2*L2*cos(theta2));

%H=T+U; % hamiltonianul sistemului - energia totala

figure(1);

Lmax=L1+L2+L3; % semilatura cadrului grafic

coef=30; % controleaza dimensiunile grafice ale corpurilor

rg1=coef*m1^(1/3); rg2=coef*m2^(1/3); rg3=coef*m3^(1/3);% raze "grafice"

tic; simt=0; % porneste cronometrul si initializeaza timpul simularii

while simt<=tf % ciclul grafic

  j=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t din discretizare

  plot([0 x1(j) x2(j) x3(j)],[0 y1(j) y2(j) y3(j)],'-g','LineWidth',3); hold on; % tije

  xlabel('x/m'); ylabel('y/m');

  plot(0,0,'.k','MarkerSize',10); % articulatie de suspensie

  plot(x1(j),y1(j),'.r','MarkerSize',rg1); % corpul 1

  plot(x1(j),y1(j),'.k','MarkerSize',10); % articulatie corp 1

  plot(x2(j),y2(j),'.b','MarkerSize',rg2); % corpul 2
  
  plot(x2(j),y2(j),'.k','MarkerSize',10); % articulatie corp 2
  
  plot(x3(j),y3(j),'.b','MarkerSize',rg3); % corpul 3

  axis([-Lmax Lmax -Lmax Lmax]); axis square; % cadrul grafic

  text(3/5*Lmax,3/5*Lmax,['E = ',num2str(round(H(j))),' J']);

  if RT==1 % real time(1) / slow motion(0)

    simt=toc; % actualizeaza timpul simularii cu ceasul sistemului

    text(3/5*Lmax,4/5*Lmax,['t = ',num2str(round(t(j))),' s']);

  else

    simt=simt+1e-1; % incrementeaza cu o decisecunda

    text(3/5*Lmax,4/5*Lmax,['t=',num2str(round(t(j)*10)),' ds']);

  end

  pause(1e-10); 
  hold off

end
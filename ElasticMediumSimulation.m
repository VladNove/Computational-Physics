%Novetschi Vlad 312CA

%realizand numeroase experimente am observat ca viteza de propagare a pulsului
%creste direct proportional cu constantele elastice are resorturilor si
%invers proportional cu masa oscilatorilor.

% Mediu elastic 1D
clear; close all; clc;
% Selectori:
TR=0; % dinamica in timp real - 1 / dinamica cu timp controlat - 0
TIP=0; % unda longitudinala - 1 / unda transversala - 0
% Parametrii fizici ai sistemului:
P=100; % numarul de corpuri
k1 = 35;
k2 = 15;
m1 = 1.5;
m2 = 0.5;
m=ones(1,P); % kg; masele corpurilor
m(1:P/2)=m1;
m(P/2:P)=m2;
k=ones(1,P+1); % N/m; constantele elastice ale resorturilor
k(1:P/2)=k1;
k(P/2:P)=k2;
a=0.5; % m; lungimile resorturilor nedeformate
Tmax=2*pi*sqrt(mean(m)/mean(k)); % timp caracteristic sistemului
ti=0; tf=20*Tmax; N=2000; t=linspace(ti,tf,N); dt=t(2)-t(1); % timp discret
eta=zeros(N,P); % prealocare deplasari (N momente de timp; P pozitii)
% Deplasari initiale si valori de start:
eta0=zeros(1,P); % m; mediu neperturbat initial
eta(1,:)=eta0; eta(2,:)=eta0; % pasul 1 si pasul 2 temporale
% Conditii la capete (frontiera):
A=2*a; % amplitudinea perturbatiei
OM=pi/2; % pulsatia perturbatiei
etas=zeros(1,N); etad=zeros(1,N); % functii de timp (capat dreapta fixat)
etas=A*sin(OM*t);
for i=2:N-1
  if (etas(i)<0)
    break;
  endif
end
etas(i-1:N)=0;
for i=2:N-1 % ciclul temporal
    for j=2:P-1 % ciclul spatial
    % Recurenta de ordinul II (vezi Curs 4):
    eta(i+1,j)=2*eta(i,j)-eta(i-1,j)+...
        dt^2/m(j)*(k(j)*(eta(i,j-1)-eta(i,j))+k(j+1)*(eta(i,j+1)-eta(i,j)));
    end % ciclul spatial
    % Relatie particulara la j=1 (primul oscilator)
    eta(i+1,1)=2*eta(i,1)-eta(i-1,1)+...
        dt^2/m(1)*(k(1)*(etas(i)-eta(i,1))+k(2)*(eta(i,2)-eta(i,1)));
    % Relatie particulara la j=P (ultimul oscilator)
    eta(i+1,P)=2*eta(i,P)-eta(i-1,P)+...
        dt^2/m(P)*(k(P)*(eta(i,P-1)-eta(i,P))+k(P+1)*(etad(i)-eta(i,P)));
end % ciclul temporal
x=a:a:P*a; % coordonatele de echilibru ale corpurilor
figure(1); % Simularea dinamica a undelor
set(1,'Position',[50,100,1200,300]); % redimensionare fereastra
tic; simt=0; % porneste cronometrul si initializeaza timpul simularii
while simt<=tf % ciclul grafic
    hold off; % sterge plot precedent
    index=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t de simt
    if TIP==1
      plot((x(1:P/2)+eta(index,1:P/2)')*[1,1],[-A,A],'-b','MarkerSize',20); hold on;
      plot((x(P/2:P)+eta(index,P/2:P)')*[1,1],[-A,A],'-g','MarkerSize',20); hold on;
      title('Unda longitudinala');
    else
      plot(x(1:P/2),eta(index,(1:P/2)),'.b','MarkerSize',20); hold on;
      plot(x(P/2:P),eta(index,(P/2:P)),'.g','MarkerSize',20); hold on;
      ylabel('eta / m');
      title('Unda transversala');
    end;
    xlabel('x / m');
    axis([x(1)-a x(P)+a,-2*A,2*A]);
    if TR==1 % 1 - in timp real
        simt=toc; % actualizeaza timpul simularii cu ceasul sistemului
        text(0.8*x(P),1.5*A,['t = ',num2str(round(t(index))),' s']);
    else
        simt=simt+5e-1; % incrementeaza timpul simularii
        text(0.8*x(P),1.5*A,['t = ',num2str(round(t(index)*10)),' ds']);
    end
    pause(1e-6)
end

    
    





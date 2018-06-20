%Bartlett Capon Pseudospectro

clear all;
close all;
clc;

M=2;            %número de elementos en el arreglo de antenas  
D=1;            %número de fuentes de las que se desea el AoA
sig2=.1;        %varianza del ruido
th1=110.99999999;   %ángulo de la fuente transmisora
%th2=60;
th1=th1*pi/180;     
%th2=th2*pi/180;
a1=[1]; 
%a2=[1];
a=[1];
for i=2:M
    a1=[a1 exp(1j*i*pi*sin(th1))];      %vector principal del arreglo
 %   a2=[a2 exp(1j*i*pi*sin(th2))];
end
%A=[a1' a2'];
A=[a1'];
Rss=[1];
Rxx=A*Rss*A'+sig2*eye(M);

theta=0:180;

for k=1:180;
   %th(k)=-pi/2+pi*k/180;
   th(k)=theta(k)*pi/180;
   clear a
   a=[1];
   for jj=2:M
      a = [a exp(1j*jj*pi*sin(th(k)))];
   end
   P_Bartlett(k)=real(a*Rxx*a');   
   P_Capon(k)=real(inv(a*inv(Rxx)*a'));
   
   %MUSIC
   [V,Dia] = eig(Rxx);
   [Y,Index] = sort(diag(Dia));
   EN = V(:,Index(1:M-D));
   
   P_MUSIC(k)=real(inv(abs(a*EN*EN'*a')));
end

figure;
subplot(311)
plot(th*180/pi,(P_Bartlett/max(P_Bartlett)),'-b')
grid on
xlabel('Angle(?)')
ylabel('|P(\theta)| (dB)')
%title('Pseudoespectro con Algoritmo de Bartlett')
subplot(312)
plot(th*180/pi,(P_Capon/max(P_Capon)),'-g')
grid on
xlabel('Angle(?)')
ylabel('|P(\theta)| (dB)')
title('Pseudoespectro con Algoritmo de Capon')
subplot(313)
plot(th*180/pi,(P_MUSIC/max(P_MUSIC)),'-r')
grid on
xlabel('Angle(?)')
ylabel('|P(\theta)| (dB)')
%title('Pseudoespectro con Algoritmo de MUSIC');

figure
plot(th*180/pi,(P_Bartlett/max(P_Bartlett)), ...
    th*180/pi,(P_Capon/max(P_Capon)), ...
    th*180/pi,(P_MUSIC/max(P_MUSIC)));
grid on
xlabel('Angle(?)')
ylabel('|P(\theta)| (dB)')
%title('Pseudoespectro con Algoritmo de MUSIC, Bartlett & Capon');
grid on;
legend('Bartlett','Capon','MUSIC');

%corr

clc;
clear all;
close all;

M = 2;             % Número de antenas
f = 1.2e9;          % Frecuencia de operación
d = [0 0.26];  % Vector de distancias según la posicicón de las antenas
c = 3e8;           % Velocidad de la luz
K = 5000;          % Número de muestras

%x = load("");

%Rxx = 1/K*(x*x'); % Obteniendo Rxx de acuerdo con SmartAntennas

Rxx = [1.1+0i 0.1006-0.9949i; 0.1006+0.9949i 1.1+0i];     %esta Rxx se carga manualmente obtenida de la simulación. Luego se cargará del archivo que se capture

theta = 0:180;

for k=1:180
  th(k) = theta(k)*pi/180;
  a = [];
  for jj=1:length(d)
    a = [a exp(1j*d(jj)*(2*pi*f/c)*sin(th(k)))];
  end
  
  P_Bartlett(k) = a*Rxx*a';
  P_Capon(k) = inv(a*inv(Rxx)*a');
  
  %MUSIC
  [V, Dia] = eig(Rxx);
  [Y, Index] = sort(diag(Dia));
  EN = V(:,Index(1:M-1));
  
  P_MUSIC(k) = inv(abs(a*EN*EN'*a')); 
end

plot(th*180/pi,(P_Bartlett/max(P_Bartlett)), ...
th*180/pi,(P_Capon/max(P_Capon)), ...
th*180/pi,(P_MUSIC/max(P_MUSIC)));

grid on;
legend("Bartlett", "Capon", "MUSIC");

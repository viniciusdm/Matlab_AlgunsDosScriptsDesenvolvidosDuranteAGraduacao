%Autor: Maria D. Miranda - Adaptação de Vinicius Bueno de Moraes, NUSP: 1O256432
%Teste 3 - PTC3451 - EPUSP 2021
%Item c. - AlgoritMO - Steepest Descent

%Limpeza do Workspace e Command Window
clc;
clear;

%Parâmetros estipulados
M = 2;
N = 1500;
mu = 0.01;

%Suporte
varS = 0.25;
theta = pi/6;
omega = (2 * pi) / 40;
phiV = 2 * pi * rand;
phiU = 2 * pi * rand;
n = 0:N-1;
A = 5;

%Sinais propostos
s = sqrt(varS) * randn(N, 1);                %Sinal de interesse
x = sin(omega * n + theta + phiV)';          %Interferência correlacionada com a entrada
u = A * sin(omega * n + phiV)';              %Sinal de entrada
d = s + x;                                   %Resposta desejada
 
r = xcorr(u, M-1, 'unbiased');
ru = r(M:end);
Ru = toeplitz(ru)                          %R - Matriz de autocorrelação
pdu = xcorr(d, u, M-1, 'unbiased');
p = pdu(M:end)                             %p - Vetor de correlação cruzada

varD = cov(d)                          
varU = cov(u)

wOtimo = Ru\p                              %Coeficientes de Wiener
JMin = varD - p' * wOtimo

yOtimo = filter(wOtimo, 1, u);             %Estimação ótima
eOtimo = d - yOtimo;                       %Erro ótimo

wOtimoN=kron(wOtimo,ones(1,N))';

%Implementação do Algoritmo LMS (Função no fim do Arquivo)
 
%Chamada da Função LMSWn para retorno de parâmetros de interesse - Definida no fim do arquivo
[yLMS, eLMS, wLMS, eaLMS] = LMSWn(M, mu, u, d, wOtimo);

%Figuras - LMS
figure(3);
suptitle({'Comparações - Solução Ótima X Algoritmo LMS com N = 1500', 'e Passo de Adapt. (\mu) de 0,01 - item c.', ''});

%Comparação entre o erro dos métodos ao longo das itereçoes
subplot(211);
plot([0:N-1],eLMS(1:N),'b'); 
hold on;
plot([0:N-1],s(1:N),'g'); 
hold off;
legend('eLMS(n)','s(n)'); 
grid;
axis([0 N-1 -2.5 2.5]);
xlabel({'n', ''});
title('Comportamento do Erro dos métodos frente ao número de iterações');

%Coeficientes - Algoritmo LMS convergindo para a solução ótima
subplot(212);
plot(wOtimoN(:,1),'m'); 
hold on;
plot(wLMS(1,1:N)','b');
plot(wOtimoN(1:N,2),'m');
plot(wLMS(2,1:N)','b');
hold off;
legend('wOtimo','wLMS(n)');
grid;
axis([0 N-1 -1.2 1.2]);
xlabel({'n', ''});
title({'', 'Coeficientes - Algoritmo LMS convergindo para a solução ótima'});

%Item c.2 - EQM X EQME

%Parâmetros da Simulação desejada
mu = 0.01;

%Declaração dos Vetores EQM e EQME
EQM = zeros(N,1); 
EQME = zeros(N,1);

% Número de realizações
K = 500; 

for k=1:K
%Variáveis a serem renovadas a cada Iteracao
phiV = 2 * pi * rand;
s = sqrt(varS) * randn(N,1);
x = sin(omega * n + theta + phiV)';
u = A * sin(omega * n + phiV);
d = s + x;

%Parâmetros - Funçao LMS
[yLMS,eLMS,wLMS,eaLMS] = LMSWn(M,mu,u,d,wOtimo);

%Cálculo EQMS
EQM = EQM + eLMS .* eLMS;
EQME = EQME + eaLMS .* eaLMS;
end

EQM=EQM/K; 
EQME=EQME/K;
EQMdB=10*log10(EQM);
EQMEdB=10*log10(EQME);
EQMEteo = JMin * ((trace(Ru)*mu)/(2-(trace(Ru)*mu)));
EQMEteodB = 10*log10(EQMEteo);

%Geração dos gráficos comparativos - EQM X EQME
figure(4);
suptitle({'Algoritmo LMS - Uma análise: EQM X EQME ao longo de iterações',''});

subplot(211); 
plot(n,EQM, 'b'); 
hold on;
plot(JMin*ones(N,1),'m--'); 
hold off;
grid; 
legend('EQM(n)', 'JMin - Wiener');
axis([0 N-1 0 2]); 
xlabel('n'); 
title('EQM LMS x JMin - Wiener');

subplot(212); 
plot(n,EQME,'b'); 
hold on;
plot(EQMEteo*ones(N,1),'m--'); 
hold off;
grid; 
legend('EQME(n)', 'EQMEmin(n)');
axis([0 N-1 0 2]);
xlabel('n'); 
title('EQME ao longo das iterações frente ao mínimo');

%Gráfico - EQM em função do Passo de Adaptação
figure(5);
suptitle({'EQM [dB] e EQME [dB] - LMS em função do Passo de Adaptação (\mu)', ''} );

subplot(211)
plot(n,EQMdB(:, 1),'r'); 
hold on;

%Parâmetros da Simulação desejada
mu = 0.02;

%Declaração dos Vetores EQM e EQME
EQM = zeros(N,1); 
EQME = zeros(N,1);

% Número de realizações
K = 500; 

for k=1:K
%Variáveis a serem renovadas a cada Iteracao
phiV = 2 * pi * rand;
s = sqrt(varS) * randn(N,1);
x = sin(omega * n + theta + phiV)';
u = A * sin(omega * n + phiV);
d = s + x;

%Parâmetros - Funçao LMS
[yLMS,eLMS,wLMS,eaLMS] = LMSWn(M,mu,u,d,wOtimo);

%Cálculo EQMS
EQM = EQM + eLMS .* eLMS;
EQME = EQME + eaLMS .* eaLMS;
end

EQM=EQM/K; 
EQME=EQME/K;
EQMdB=10*log10(EQM);
EQMEdB=10*log10(EQME);
EQMEteo = JMin * ((trace(Ru)*mu)/(2-(trace(Ru)*mu)));
EQMEteodB = 10*log10(EQMEteo);

plot(n,EQMdB,'b'); 

%Parâmetros da Simulação desejada
mu = 0.03;

%Declaração dos Vetores EQM e EQME
EQM = zeros(N,1); 
EQME = zeros(N,1);

% Número de realizações
K = 500; 

for k=1:K
%Variáveis a serem renovadas a cada Iteracao
phiV = 2 * pi * rand;
s = sqrt(varS) * randn(N,1);
x = sin(omega * n + theta + phiV)';
u = A * sin(omega * n + phiV);
d = s + x;

%Parâmetros - Funçao LMS
[yLMS,eLMS,wLMS,eaLMS] = LMSWn(M,mu,u,d,wOtimo);

%Cálculo EQMS
EQM = EQM + eLMS .* eLMS;
EQME = EQME + eaLMS .* eaLMS;
end

EQM=EQM/K; 
EQME=EQME/K;
EQMdB=10*log10(EQM);
EQMEdB=10*log10(EQME);
EQMEteo = JMin * ((trace(Ru)*mu)/(2-(trace(Ru)*mu)));
EQMEteodB = 10*log10(EQMEteo);

plot(n,EQMdB(:,1),'g');
plot(n,10*log10(JMin)*ones(N,1),'m--'); 
hold off;
grid; 
legend('\mu=0,01','\mu=0,02','\mu=0,03','JMin');
ylim([-6.5 -1.5]);
ylabel('EQM [dB]');
xlabel('n');
title('Simulaçao - EQM [dB] para \mu=0,01, \mu=0,02 e \mu=0,03');

subplot(212)

%Parâmetros da Simulação desejada 
mu = 0.01;

%Declaração dos Vetores EQM e EQME
EQM = zeros(N,1); 
EQME = zeros(N,1);

% Número de realizações
K = 500; 

for k=1:K
%Variáveis a serem renovadas a cada Iteracao
phiV = 2 * pi * rand;
s = sqrt(varS) * randn(N,1);
x = sin(omega * n + theta + phiV)';
u = A * sin(omega * n + phiV);
d = s + x;

%Parâmetros - Funçao LMS
[yLMS,eLMS,wLMS,eaLMS] = LMSWn(M,mu,u,d,wOtimo);

%Cálculo EQMS
EQM = EQM + eLMS .* eLMS;
EQME = EQME + eaLMS .* eaLMS;
end

EQM=EQM/K; 
EQME=EQME/K;
EQMdB=10*log10(EQM);
EQMEdB=10*log10(EQME);
EQMEteo = JMin * ((trace(Ru)*mu)/(2-(trace(Ru)*mu)));
EQMEteodB = 10*log10(EQMEteo);

plot(n,EQMEdB,'r', 'DisplayName', '\mu=0,01'); 
hold on;
plot(n,EQMEteodB*ones(N,1),'r--', 'DisplayName', '\mu=0,01 teo');

%Parâmetros da Simulação desejada
mu = 0.02;

%Declaração dos Vetores EQM e EQME
EQM = zeros(N,1); 
EQME = zeros(N,1);

% Número de realizações
K = 500; 

for k=1:K
%Variáveis a serem renovadas a cada Iteracao
phiV = 2 * pi * rand;
s = sqrt(varS) * randn(N,1);
x = sin(omega * n + theta + phiV)';
u = A * sin(omega * n + phiV);
d = s + x;

%Parâmetros - Funçao LMS
[yLMS,eLMS,wLMS,eaLMS] = LMSWn(M,mu,u,d,wOtimo);

%Cálculo EQMS
EQM = EQM + eLMS .* eLMS;
EQME = EQME + eaLMS .* eaLMS;
end

EQM=EQM/K; 
EQME=EQME/K;
EQMdB=10*log10(EQM);
EQMEdB=10*log10(EQME);
EQMEteo = JMin * ((trace(Ru)*mu)/(2-(trace(Ru)*mu)));
EQMEteodB = 10*log10(EQMEteo);

plot(n,EQMEdB,'b', 'DisplayName', '\mu=0,02');
plot(n,EQMEteodB*ones(N,1),'b--', 'DisplayName', '\mu=0,02 teo'); 

%Parâmetros da Simulação desejada 
mu = 0.03;

%Declaração dos Vetores EQM e EQME
EQM = zeros(N,1); 
EQME = zeros(N,1);

% Número de realizações
K = 500; 

for k=1:K
%Variáveis a serem renovadas a cada Iteracao
phiV = 2 * pi * rand;
s = sqrt(varS) * randn(N,1);
x = sin(omega * n + theta + phiV)';
u = A * sin(omega * n + phiV);
d = s + x;

%Parâmetros - Funçao LMS
[yLMS,eLMS,wLMS,eaLMS] = LMSWn(M,mu,u,d,wOtimo);

%Cálculo EQMS
EQM = EQM + eLMS .* eLMS;
EQME = EQME + eaLMS .* eaLMS;
end

EQM=EQM/K; 
EQME=EQME/K;
EQMdB=10*log10(EQM);
EQMEdB=10*log10(EQME);
EQMEteo = JMin * ((trace(Ru)*mu)/(2-(trace(Ru)*mu)));
EQMEteodB = 10*log10(EQMEteo);

plot(n,EQMEdB,'g', 'DisplayName', '\mu=0,03');
plot(n,EQMEteodB*ones(N,1),'g--', 'DisplayName', '\mu=0,013 teo'); 
hold off;
grid;
legend();
ylim([-15 -4]);
ylabel('EQME [dB]');
xlabel('n');
title('Simulaçao - EQME [dB] X Teórioco, para \mu = 0,01, \mu = 0,02 e \mu = 0,03');


%Função - LMS
function [yLMS, eLMS, wLMS, eaLMS]=LMSWn(M, mu, u, d, wOtimo)
N=length(u); 
wLMS=zeros(M,1); 
yLMS=zeros(N,1);
eLMS=zeros(N,1); 
uv=zeros(M,1);
eaLMS=zeros(N, 1);

for i=1:N
uv=[u(i); uv(1:M-1)];
yLMS(i)=uv'*wLMS(:,i);
eLMS(i)=d(i)-yLMS(i);
wLMS(:,i+1)=wLMS(:,i)+mu*eLMS(i)*uv;
eaLMS(i)=uv'*(wOtimo-wLMS(:,i));

end
end

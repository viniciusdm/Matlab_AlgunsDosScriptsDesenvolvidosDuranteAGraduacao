%Autor: Maria D. Miranda - Adaptação de Vinicius Bueno de Moraes, NUSP: 1O256432
%Teste 3 - PTC3451 - EPUSP 2021
%Item b. - AlgoritMO - Steepest Descent

%Limpeza do Workspace e Command Window
clc;
clear;

%Parâmetros estipulados
M = 2;
N = 1500;

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

wOtimo = Ru\p                              %Coeficientes de Wiener
JMin = varD - p' * wOtimo

yOtimo = filter(wOtimo, 1, u);             %Estimação ótima
eOtimo = d - yOtimo;                       %Erro ótimo

%Faixa de valores para convergência do Algoritmo - Em função de R - item
%b.1
maxAutValor = max(eig(Ru));
maxPasso = 2/maxAutValor

%item b.2
N=1500; 
mu=0.01;

wSD=zeros(M,N); 
ySD=zeros(N,1); 
eSD=zeros(N,1);
JSD=zeros(N,1); 
uv=zeros(M,1);

for i=1:N
wSD(:,i+1)=wSD(:,i)+mu*(p-Ru*wSD(:,i));
uv=[u(i); uv(1:M-1)];
ySD(i)=uv'*wSD(:,i);
eSD(i)=d(i)-ySD(i);
JSD(i)=varD-2*wSD(:,i)'*p+wSD(:,i)'*Ru*wSD(:,i);
end

JMinSteepest=varD-2*wOtimo'*p+wOtimo'*Ru*wOtimo
wOtimoN=kron(wOtimo,ones(1,N))';

%Figuras - steepest descent
figure(2);
suptitle({'Comparações - Solução Ótima X Algoritmo de Steepest Descent com', 'N = 1500 e Passo de Adapt. (\mu) de 0,01 - item b.', ''});

%Comparação entre o erro dos métodos ao longo das itereç?es
subplot(311);
plot([0:N-1],eSD(1:N),'b'); 
hold on;
plot([0:N-1],s(1:N),'g'); 
hold off;
legend('eSD(n)','s(n)'); 
grid;
axis([0 N-1 -2.5 2.5]);
xlabel({'n', ''});
title('Comportamento do Erro dos métodos frente ao número de iterações');

%Coeficientes de Steepest Descent convergindo para a solução ótima
subplot(312);
plot(wOtimoN(:,1),'m'); 
hold on;
plot(wSD(1,1:N)','b');
plot(wOtimoN(1:N,2),'m');
plot(wSD(2,1:N)','b');
hold off;
legend('wOtimo','wSD(n)');
grid;
axis([0 N-1 -1.2 1.2]);
xlabel({'n', ''});
title({'', 'Coeficientes de Steepest Descent convergindo para a solução ótima'});


%Comparativo - EQM dos métodos
subplot(313);
plot(JSD,'b'); 
hold on;
plot(JMinSteepest*ones(N,1),'m'); 
hold off;
legend('JSD(n)','JMin - Wiener');
grid;
axis([0 N-1 0 2]); 
xlabel({'n', ''});
title('Item b.3 - EQM das Soluções - Convergindo ao longo das iterações');

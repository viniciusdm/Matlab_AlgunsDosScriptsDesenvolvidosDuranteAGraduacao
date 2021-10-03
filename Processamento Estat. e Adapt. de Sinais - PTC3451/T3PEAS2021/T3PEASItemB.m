%Autor: Maria D. Miranda - Adapta��o de Vinicius Bueno de Moraes, NUSP: 1O256432
%Teste 3 - PTC3451 - EPUSP 2021
%Item b. - AlgoritMO - Steepest Descent

%Limpeza do Workspace e Command Window
clc;
clear;

%Par�metros estipulados
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
x = sin(omega * n + theta + phiV)';          %Interfer�ncia correlacionada com a entrada
u = A * sin(omega * n + phiV)';              %Sinal de entrada
d = s + x;                                   %Resposta desejada
 
r = xcorr(u, M-1, 'unbiased');
ru = r(M:end);
Ru = toeplitz(ru)                          %R - Matriz de autocorrela��o
pdu = xcorr(d, u, M-1, 'unbiased');
p = pdu(M:end)                             %p - Vetor de correla��o cruzada

varD = cov(d)                          

wOtimo = Ru\p                              %Coeficientes de Wiener
JMin = varD - p' * wOtimo

yOtimo = filter(wOtimo, 1, u);             %Estima��o �tima
eOtimo = d - yOtimo;                       %Erro �timo

%Faixa de valores para converg�ncia do Algoritmo - Em fun��o de R - item
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
suptitle({'Compara��es - Solu��o �tima X Algoritmo de Steepest Descent com', 'N = 1500 e Passo de Adapt. (\mu) de 0,01 - item b.', ''});

%Compara��o entre o erro dos m�todos ao longo das itere�?es
subplot(311);
plot([0:N-1],eSD(1:N),'b'); 
hold on;
plot([0:N-1],s(1:N),'g'); 
hold off;
legend('eSD(n)','s(n)'); 
grid;
axis([0 N-1 -2.5 2.5]);
xlabel({'n', ''});
title('Comportamento do Erro dos m�todos frente ao n�mero de itera��es');

%Coeficientes de Steepest Descent convergindo para a solu��o �tima
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
title({'', 'Coeficientes de Steepest Descent convergindo para a solu��o �tima'});


%Comparativo - EQM dos m�todos
subplot(313);
plot(JSD,'b'); 
hold on;
plot(JMinSteepest*ones(N,1),'m'); 
hold off;
legend('JSD(n)','JMin - Wiener');
grid;
axis([0 N-1 0 2]); 
xlabel({'n', ''});
title('Item b.3 - EQM das Solu��es - Convergindo ao longo das itera��es');

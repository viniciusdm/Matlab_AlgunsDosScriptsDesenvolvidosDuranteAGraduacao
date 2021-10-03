%Autor: Maria D. Miranda - Adaptação de Vinicius Bueno de Moraes, NUSP: 1O256432
%Teste 3 - PTC3451 - EPUSP 2021
%Item a. - Filtro ótimo linear e Solução de Wiener

%Limpeza do Workspace e Command Window
clc;
clear;

%Parâmetros estipulados
M = 2;
N = 4500;

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

%Geração da figura enunciando o comportomamento do Filtro Adaptativo
N2 = 200;

figure(1);
suptitle({'Resultado da Filtragem Adaptativa para os sinais propostos - item a.',''});

%Análise inicial dos sinais envolvidos
subplot(311);
plot([0:N2-1],u(1:N2),'g');
hold on;
plot([0:N2-1],d(1:N2),'k'); 
hold off;
grid; 
legend('u(n)','d(n)');
xlabel({'n', ''});
axis([0 N2-1 -5.1 5.1]);
title('Primeira análise - Sinal de entrada frente a Resposta do Filtro Desejada d(n)');

%Saída esperada do Filtro Adaptativo
subplot(312);
plot([0:N2-1],x(1:N2),'g'); 
hold on;
plot([0:N2-1],d(1:N2),'k');
plot([0:N2-1],yOtimo(1:N2),'m.'); 
hold off;
grid; 
legend('x(n)','d(n)','yOtimo(n)');
xlabel({'n', ''}); 
axis([0 N2-1 -2.2 2.2]);
title('Saída do Filtro se sobrepondo exatamente ao sinal correlacionado a entrada u(n) - (x(n))');

%Eerro do Filtro Adaptativo
subplot(313);
plot([0:N2-1],d(1:N2),'k'); 
hold on;
plot([0:N2-1],s(1:N2),'g'); 
plot([0:N2-1],eOtimo(1:N2),'m.'); 
hold off;
grid; 
legend('d(n)','s(n)','eOtimo(n)');
xlabel({'n', ''}); 
axis([0 N2-1 -2.2 2.2]);
title('Erro ótimo se sobrepondo exatamente ao Sinal descorrelacionado da entrada u(n) - (s(n))');

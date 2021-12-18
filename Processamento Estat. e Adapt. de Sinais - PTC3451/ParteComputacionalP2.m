%Escola Politécnica da Universidade de São Paulo
%PTC3451 - Processamento Estatístico e Adaptativo de Sinais
%Prof. Maria D. Miranda - EPUSP, 18/12/2021
%Autor: Vinicius Bueno de Moraes, Adaptado de Maria D. Miranda | NUSP: 10256432

clc; clear all;  %Limpeza do ambiente
%Segunda Prova - Parte Computacional
disp('Segunda Prova - Parte Computacional');
disp(' ');

%Questão 1 - itens ((a), (b) e (c)) - Rotinas Auxiliares
disp('Questão 1 - itens ((a), (b) e (c)) - Rotinas Auxiliares');
disp('');

RxTeorico = [21 19.39 17.59; 19.39 21 19.39; 17.59 19.39 21]
[QTeorico, DTeorico] = eig(RxTeorico);
VarianciaTeorico = diag(DTeorico)
PropVarExplicadaTeorico = 100*diag(DTeorico)/sum(diag(DTeorico))

%Questão 1 - item (d) - Parte Computacional
disp('Questão 1 - item (d) - Parte Computacional');

load dadosAP2.mat;                                                     %Carrega os dados fornecidos (x1, x2 e x3)

%Estimativas e composicão do Modelo PCA
X = [x1; x2; x3];                                                      %Composicão da Matriz X

RxEstimado = X*X'/size(X, 2)                                           %Estimativa da Matriz de Covariância
[QEstimado, DEstimado] = eig(RxEstimado);
VarianciaEstimado = diag(DEstimado)
PropVarExplicadaEstimado = 100*diag(DEstimado)/sum(diag(DEstimado))    %Variância Explicada
Cc = QEstimado*DEstimado*QEstimado';

y1 = QEstimado(:, 1)'*X;
y2 = QEstimado(:, 2)'*X;
y3 = QEstimado(:, 3)'*X;

%Plots - Geracão dos Gráficos
figure(1);
suptitle({'Questão 1 - item (d) iii.', '',''});

%Gráfico dos sinais x(n) ao longo do tempo 
subplot(321);
plot(x1, 'b');
legend('x1(n)');
xlabel('n');
grid;
axis([0 size(X,2) -15 15]);
title('Sinal x1(n) ao longo do tempo');

subplot(323);
plot(x2, 'b');
legend('x2(n)');
xlabel('n');
grid;
axis([0 size(X,2) -15 15]);
title('Sinal x2(n) ao longo do tempo');

subplot(325);
plot(x3, 'b');
legend('x3(n)');
xlabel('n');
grid;
axis([0 size(X,2) -15 15]);
title('Sinal x3(n) ao longo do tempo');

%Gráfico dos sinais y(n) ao longo do tempo 
subplot(322);
plot(y1, 'b');
legend('y1(n)');
xlabel('n');
grid;
axis([0 size(X,2) -15 15]);
title('Sinal y1(n) ao longo do tempo');

subplot(324);
plot(y2, 'b');
legend('y2(n)');
xlabel('n');
grid;
axis([0 size(X,2) -15 15]);
title('Sinal y2(n) ao longo do tempo');

subplot(326);
plot(y3, 'b');
legend('y3(n)');
xlabel('n');
grid;
axis([0 size(X,2) -15 15]);
title('Sinal y3(n) ao longo do tempo');

%Gráfico - Diagrama de Dispersão
figure(2);
suptitle({'Questão 1 - item (d) iv.', '',''});

subplot(321);
scatter(x1, x2, 1);
hold on
plot(2*QEstimado(1,3)*[-10 10], 2*QEstimado(2,3)*[-10 10],'g','Linewidth',2)
plot(2*QEstimado(1,2)*[-10 10], 2*QEstimado(2,2)*[-10 10],'m','Linewidth',2)
plot(2*QEstimado(1,1)*[-10 10], 2*QEstimado(2,1)*[-10 10],'r','Linewidth',2)
hold off;
grid;
xlabel('x1(n)');
ylabel('x2(n)');
legend('scatter plot','CP1', 'CP2', 'CP3');
axis([-20 20 -20 20]);
title('Diagrama de dispersão de x1(n) por x2(n)');

subplot(323);
scatter(x1, x3, 1);
hold on
plot(2*QEstimado(1,3)*[-10 10], 2*QEstimado(3,3)*[-10 10],'g','Linewidth',2)
plot(2*QEstimado(1,2)*[-10 10], 2*QEstimado(3,2)*[-10 10],'m','Linewidth',2)
plot(2*QEstimado(1,1)*[-10 10], 2*QEstimado(3,1)*[-10 10],'r','Linewidth',2)
hold off; 
grid;
xlabel('x1(n)');
ylabel('x3(n)');
legend('scatter plot','CP1', 'CP2', 'CP3');
axis([-20 20 -20 20]);
title('Diagrama de dispersão de x1(n) por x3(n)');

subplot(325);
scatter(x2, x3, 1);
hold on
plot(2*QEstimado(2,3)*[-10 10], 2*QEstimado(3,3)*[-10 10],'g','Linewidth',2)
plot(2*QEstimado(2,2)*[-10 10], 2*QEstimado(3,2)*[-10 10],'m','Linewidth',2)
plot(2*QEstimado(2,1)*[-10 10], 2*QEstimado(3,1)*[-10 10],'r','Linewidth',2)
hold off; 
grid;
xlabel('x2(n)');
ylabel('x3(n)');
legend('scatter plot','CP1', 'CP2', 'CP3');
axis([-20 20 -20 20]);
title('Diagrama de dispersão de x2(n) por x3(n)');

subplot(322);
scatter(y1, y2, 1);
grid;
xlabel('y1(n)');
ylabel('y2(n)');
axis([-20 20 -20 20]);
title('Diagrama de dispersão de y1(n) por y2(n)');

subplot(324);
scatter(y1, y3, 1);
grid;
xlabel('y1(n)');
ylabel('y3(n)');
axis([-20 20 -20 20]);
title('Diagrama de dispersão de y1(n) por y3(n)');

subplot(326);
scatter(y2, y3, 1);
grid;
xlabel('y2(n)');
ylabel('y3(n)');
axis([-20 20 -20 20]);
title('Diagrama de dispersão de y2(n) por y3(n)');

%Gráfico - Histogramas
figure(3);
suptitle({'Questão 1 - item (d) v.', '',''});

subplot(321);
histogram(x1);
grid;
axis([-15 15 0 800]);
title('Histograma de x1(n)');

subplot(323);
histogram(x2);
grid;
axis([-15 15 0 800]);
title('Histograma de x2(n)');

subplot(325);
histogram(x3);
grid;
axis([-15 15 0 800]);
title('Histograma de x3(n)');

subplot(322);
histogram(y1);
grid;
axis([-15 15 0 800]);
title('Histograma de y1(n)');

subplot(324);
histogram(y2);
grid;
axis([-15 15 0 800]);
title('Histograma de y2(n)');

subplot(326);
histogram(y3);
grid;
axis([-15 15 0 800]);
title('Histograma de y3(n)');

%Questão 2 - item (a) - Pedriodograma e DEP de x1(n)
N=size(x1, 2);
L=size(x1, 2);          %Comprimento da Janela - Pede-se o número de amostras do sinal
K=1;                    %Existira apenas um seguimento, das especificacoes e item anterior
Nf=N;                   
Nf2=round(Nf/2);

Ipm1=zeros(1,Nf);
janela1=rectwin(L)';
 
for r=0:K-1
 xr1(1:L)=x1(r+[1:L]).*janela1; 
 Ipm1=Ipm1+abs(fft(xr1(1:L),Nf)).^2;
end

U1=sum(janela1.^2)/L; 
Um1=sum(janela1)^2;  

%Estimativa da DSP
Ipm1=Ipm1/(L*K*U1);  
 
%Gráfico - (DEP)
figure(4)
plot([0:Nf2-1],10*log10(Ipm1(1:Nf2)),'r');
grid;
title({'Questão 2 - item (a)', '', 'Densidade Espectral de Potência de x1(n) [dB]'});
xlabel('Frequência');
legend('Janela retangular');

%Questão 2 - item (b) - Determinação do Sinal
disp('Questão 2 - itens (b) - Determinação do Sinal');
disp(' ');
varR1E = mean(Ipm1(1000:Nf/2))                %Variância - ruído a apartir de 1000
[v11,v12] = max(abs(Ipm1));                   %Intermediário 1 – freq. e amplitude
[v21,v22] = max(abs(Ipm1(Ipm1<max(Ipm1))));   %Intermediário 2 – freq. e amplitude
w1sinalEstimado1 = (2*pi*v12/Nf)*(180/3.1415) %Freq senoide 1 (graus)
w2sinalEstimado1 = (2*pi*v22/Nf)*(180/3.1415) %Freq senoide 2 (graus)
A1 = sqrt((v11-varR1E)*4*U1*L/(Um1))          %Amplitude senoide 1
A2 = sqrt((v21-varR1E)*4*U1*L/(Um1))          %Amplitude senoide 2










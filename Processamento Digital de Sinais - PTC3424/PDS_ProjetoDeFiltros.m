%Escola Polit�cnica da Universidade de S�o Paulo | PTC3424 - Processamento Digital de Sinais
%Docente: Maria D. Miranda
%Nome: Vinicius Bueno de Moraes - NUSP: 10256432

%Referente a Quet�o 2 da Prova Computacional (P2)

%%%%Script para esbo�o dos Diagramas de Toler�ncia, de Ambos os Filtros%%%%

clear all;

%Par�metros de entrada
freq = 0:0.001:1; 
passagemA = zeros(1,length(freq));
passagemB = zeros(1,length(freq));
tolsAinf = zeros(1,length(freq));
tolsAsup = zeros(1,length(freq));
tolsBinf = zeros(1,length(freq));
tolsBsup = zeros(1,length(freq));

%Composi��o do Diagrama do Filtro A
for i=1:length(freq)
    if (freq(i) >= 0.36) && (freq(i) < 0.76)
        passagemA(i) = 1;
    else
        passagemA(i) = 0;
    end
end

for i=1:length(freq)
    if (freq(i) >= 0) && (freq(i) < 0.26)
        tolsAinf(i) = 0.02;
        tolsAsup(i) = 0.02;
    elseif (freq(i) > 0.46) && (freq(i) < 0.66)
        tolsAinf(i) = 0.99;
        tolsAsup(i) = 1.01;
    elseif (freq(i) > 0.86) && (freq(i) < 1)
        tolsAinf(i) = 0.02;
        tolsAsup(i) = 0.02;
    end
end

%Composi��o do Diagrama do Filtro B
for i=1:length(freq)
    if (freq(i) >= 0.25) && (freq(i) < 0.67)
        passagemB(i) = 1;
    else
        passagemB(i) = 0;
    end
end

for i=1:length(freq)
    if (freq(i) >= 0) && (freq(i) < 0.24)
        tolsBinf(i) = 0.1;
        tolsBsup(i) = 0.1;
    elseif (freq(i) > 0.26) && (freq(i) < 0.66)
        tolsBinf(i) = 0.9;
        tolsBsup(i) = 1.1;
    elseif (freq(i) > 0.68) && (freq(i) < 1)
        tolsBinf(i) = 0.1;
        tolsBsup(i) = 0.1;
    end
end

%Esbo�o do Diagrama de Toler�ncias do Filtro A
figure(1);
plot(freq, passagemA, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsAsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsAinf, '--k', 'LineWidth', 1);
xlim([0 1]);
ylim([0 1.2]);
title('Diagrama de Toler�ncias - Filtro A');
legend('Ganho', 'Ripple - Passagem', 'Ripple - Rejei��o');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

%Esbo�o do Diagrama de Toler�ncias do Filtro B
figure(2);
plot(freq, passagemB, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsBsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsBinf, '--k', 'LineWidth', 1);
xlim([0 1]);
ylim([0 1.2]);
title('Diagrama de Toler�ncias - Filtro B');
legend('Ganho', 'Ripple - Passagem', 'Ripple - Rejei��o');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

%Compara��o entres os Filtros
figure(3);
plot(freq, passagemA, 'r', 'LineWidth', 2);
hold on;
plot(freq, passagemB, 'k', 'LineWidth', 2);
xlim([0 1]);
ylim([0 1.2]);
title('Compara��o entre os Filtros');
legend('Ganho - Filtra A', 'Ganho - Filtro B');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

%Gera��o das Janelas

%Tamanho das Janelas Obtidas no Documento Te�rico que referencia esse + 1

%Filtro A - Janela de Hamming - N = 40
NAHamming = 41;
LAHamming = (NAHamming-1)/2;

%Filtro A - Janela de Blackman - N = 60
NABlackman = 61;
LABlackman = (NABlackman-1)/2;

%Filtro B  - Janela de Retangular - N = 200
NBRetangular = 201;
LBRetangular = (NBRetangular-1)/2;

%Filtro B - Janela de Hamming - N = 400
NBHamming = 401;
LBHamming = (NBHamming-1)/2;

%Filtro B - Janela de Blackman - N = 600
NBBlackman = 601;
LBBlackman =(NBBlackman-1)/2;

%%%%

%Frequ�ncias de Corte - Filtra A
omegaCA1 = 0.36*pi;
omegaCA2 = 0.76*pi;

%Frequ�ncias de Corte - Filtro B
omegaCB1 = 0.25*pi;
omegaCB2 = 0.67*pi;

%%%%

%n's
nAH = -LAHamming : LAHamming;
nAB = -LABlackman : LABlackman;
nBR = -LBRetangular : LBRetangular;
nBH = -LBHamming : LBHamming;
nBB = -LBBlackman : LBBlackman;

%%%%

%Defini��o da respsota ao pulso unit�rio do Filtro A
hdAH = (omegaCA2/pi)*sinc((omegaCA2/pi)*nAH)-(omegaCA1/pi)*sinc((omegaCA1/pi)*nAH);
hdAB = (omegaCA2/pi)*sinc((omegaCA2/pi)*nAB)-(omegaCA1/pi)*sinc((omegaCA1/pi)*nAB);

%Defini��o da respsota ao pulso unit�rio do Filtro B
hdBR = (omegaCB2/pi)*sinc((omegaCB2/pi)*nBR)-(omegaCB1/pi)*sinc((omegaCB1/pi)*nBR);
hdBH = (omegaCB2/pi)*sinc((omegaCB2/pi)*nBH)-(omegaCB1/pi)*sinc((omegaCB1/pi)*nBH);
hdBB = (omegaCB2/pi)*sinc((omegaCB2/pi)*nBB)-(omegaCB1/pi)*sinc((omegaCB1/pi)*nBB);

%%%%

%Janela de Hamming - N = 40 - Filtro A;
janelaHA = hamming(NAHamming)';
hoAH = hdAH.*janelaHA;

%Janela de Blackman - N = 60 - Filtro A;
janelaAB = blackman(NABlackman)';
hoAB = hdAB.*janelaAB;

%Janela Retangular - N = 200 - Filtro B;
janelaBR = rectwin(NBRetangular)';
hoBR = hdBR.*janelaBR;

%Janela de Hamming - N = 400 - Filtro B;
janelaBH = hamming(NBHamming)';
hoBH = hdBH.*janelaBH;

%Janela de Blackman - N = 600 - Filtro B;
janelaBB = blackman(NBBlackman)';
hoBB = hdBB.*janelaBB;

%%%%

[HoAH,omegaAH]=freqz(hoAH,1);
[HoAB,omegaAB]=freqz(hoAB,1);

[HoBR,omegaBR]=freqz(hoBR,1);
[HoBH,omegaBH]=freqz(hoBH,1);
[HoBB,omegaBB]=freqz(hoBB,1);

figure(4);
plot(freq, passagemA, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsAsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsAinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaAH/pi,abs(HoAH), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela de Hamming - N = 40 (wc1 = 0,36pi e wc2 = 0,76pi) - Filtro A (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela de Hamming');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

figure(5);
plot(freq, passagemA, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsAsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsAinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaAB/pi,abs(HoAB), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela de Blackman - N = 60 (wc1 = 0,36pi e wc2 = 0,76pi) - Filtro A (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela de Blackman');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

figure(6);
plot(freq, passagemB, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsBsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsBinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaBR/pi,abs(HoBR), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela Retangular - N = 200 (wc1 = 0,25pi e wc2 = 0,67pi) - Filtro B (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela Retangular');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

figure(7);
plot(freq, passagemB, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsBsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsBinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaBH/pi,abs(HoBH), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela de Hamming - N = 400 (wc1 = 0,25pi e wc2 = 0,67pi) - Filtro B (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela de Hamming');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

figure(8);
plot(freq, passagemB, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsBsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsBinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaBB/pi,abs(HoBB), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela de Blackman - N = 600 (wc1 = 0,25pi e wc2 = 0,67pi) - Filtro B (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela de Blackman');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

%%%%Rotinas Referente a Janela de Kaiser

%Par�metros Referentes ao Filtro A
NA = 31;
LA = NA/2;
nAK = 0:NA-1;
betaA = 4,601;

%Par�metros Referentes ao Filtro B
NB = 85;
LB = NB/2;
nBK = 0:NB-1;
betaB = 0;

%Frequ�ncias de Corte - Filtra A
omegaCA1 = 0.36*pi;
omegaCA2 = 0.76*pi;

%Frequ�ncias de Corte - Filtro B
omegaCB1 = 0.25*pi;
omegaCB2 = 0.67*pi;

%%%%

%Defini��o da respsota ao pulso unit�rio do Filtro A
hdAK = (omegaCA2/pi)*sinc((omegaCA2/pi)*(nAK-LA))-(omegaCA1/pi)*sinc((omegaCA1/pi)*(nAK-LA));

%Defini��o da respsota ao pulso unit�rio do Filtro B
hdBK = (omegaCB2/pi)*sinc((omegaCB2/pi)*(nBK-LB))-(omegaCB1/pi)*sinc((omegaCB1/pi)*(nBK-LB));

%Janela de Kaiser - Filtro A
gKA = kaiser(NA, betaA)';
hoAK = hdAK.*gKA;

%Janela de Kaiser - Filtro A
gKB = kaiser(NB, betaB)';
hoBK = hdBK.*gKB;

[HoAK,omegaAK]=freqz(hoAK,1);

[HoBK,omegaBK]=freqz(hoBK,1);

figure(9);
plot(freq, passagemA, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsAsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsAinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaAK/pi,abs(HoAK), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela de Kaiser - N = 31 (wc1 = 0,36pi e wc2 = 0,76pi) - Filtro A (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela de Kaiser');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

figure(10);
plot(freq, passagemB, 'r', 'LineWidth', 2);
hold on;
plot(freq, tolsBsup, '--k', 'LineWidth', 1);
hold on;
plot(freq, tolsBinf, '--k', 'LineWidth', 1);
hold on;
plot(omegaBK/pi,abs(HoBK), 'b', 'LineWidth', 2);
hold on;
xlim([0 1]);
ylim([0 1.2]);
suptitle('M�dulo da Resposta em Frequ�ncia Sobreposto ao Diagrama de Toler�ncias');
title('Janela de Kaiser - N = 85 (wc1 = 0,25pi e wc2 = 0,67pi) - Filtro B (Passa-faixas)');
legend('Ganho do Filtro', 'Ripple - Passagem', 'Ripple - Rejei��o', 'Janela de Kaiser');
xlabel('Frequ�ncia - pi[rad/s]');
ylabel('M�dulo da Resposta em Frequ�ncia - |H(exp(jw))|');

%M�dulo e Fase - Filtro A - Janela de Kaiser vs Janela de Hamming
figure(11);
freqz(hoAK,1);
title('Janela de Kaiser - Filtro A');
grid on;

figure(12);
freqz(hoAH,1);
title('Janela de Hamming - Filtro A');
grid on;

%M�dulo e Fase - Filtro B - Janela de Kaiser vs Janela Retangular
figure(13);
freqz(hoBK,1);
title('Janela de Kaiser - Filtro B');
grid on;

figure(14);
freqz(hoBR,1);
title('Janela Retangular - Filtro B');
grid on;

%Fim
fprintf('Script Executado com Sucesso');
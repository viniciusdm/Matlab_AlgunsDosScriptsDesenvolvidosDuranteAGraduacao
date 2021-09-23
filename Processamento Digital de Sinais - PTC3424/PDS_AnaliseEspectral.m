%Escola Politécnica da Universidade de São Paulo | PTC3424 - Processamento Digital de Sinais
%Docente: Maria D. Miranda
%Nome: Vinicius Bueno de Moraes - NUSP: 10256432

%Referente a Quetão 1 da Prova Computacional (P2)

%%%%Script para Gráfico das FFTs e afins%%%%

clear all;

%Carrega o sinal fornecido
load dadosprova2021.mat;
N = length(x); 
n = [0:N-1];

%Define os Ms, M1 = N e M2 = 2*N
M1 = N;
M2 = 2*N;

%Frequência de amostragem utilizando periodo de amostragem dado no
%enunciado.
fa = 1/10^-3; 

%Janela Retamgular - Sinal x(n);
janelaR = rectwin(N)';
XjanelaRetangular = x.*janelaR;

%Janela de Hamming - Sinal x(n);
janelaH = hamming(N)';
XjanelaHamming = x.*janelaH;

%Janela de Blackman - Sinal x(n);
janelaB = blackman(N)';
XjanelaBlackman = x.*janelaB;

%Comando fftshift - TFD no intervalo de -pi a pi rad - Janela Retangular
%M = N
TFDNXjanelaRetangular1 = fftshift(fft(XjanelaRetangular, M1));
TFDNXjanelaHamming1 = fftshift(fft(XjanelaHamming, M1));
TFDNXjanelaBlackman1 = fftshift(fft(XjanelaBlackman, M1));
%M = 2*N
TFDNXjanelaRetangular2 = fftshift(fft(XjanelaRetangular, M2));
TFDNXjanelaHamming2 = fftshift(fft(XjanelaHamming, M2));
TFDNXjanelaBlackman2 = fftshift(fft(XjanelaBlackman, M2));

%Correções de Escala - Janela Retangular
%M = N
TFDNXjanelaRetangular1 = 2*TFDNXjanelaRetangular1/sum(rectwin(N));
TFDNXjanelaHamming1 = 2*TFDNXjanelaHamming1/sum(hamming(N));
TFDNXjanelaBlackman1 = 2*TFDNXjanelaBlackman1/sum(blackman(N));
%M = 2N
TFDNXjanelaRetangular2 = 2*TFDNXjanelaRetangular2/sum(rectwin(N));
TFDNXjanelaHamming2 = 2*TFDNXjanelaHamming2/sum(hamming(N));
TFDNXjanelaBlackman2 = 2*TFDNXjanelaBlackman2/sum(blackman(N));
%Correção de Escala
ff1 = fa*((0:M1-1)-ceil((M1-1)/2))/M1;
ff2 = fa*((0:M2-1)-ceil((M2-1)/2))/M2;

% Janela Retangular
figure(1);
subplot(221);
plot(ff1, abs(TFDNXjanelaRetangular1),'k');
title(['FFT com M = N pontos']);
ylabel('|Xjanela Retangular(f)|');
grid on;
%
subplot(223);
plot(ff1, 20*log10(abs(TFDNXjanelaRetangular1)),'k'); 
title(['FFT com M = N pontos']);
ylabel('|Xjanela Retangular(f)| dB');
grid on; 
%
subplot(222)
plot(ff2, abs(TFDNXjanelaRetangular2), 'k');
title(['FFT com M = 2N pontos']);
ylabel('|Xjanela Retangular(f))|');
grid on; 
%
subplot(224)
plot(ff2, 20*log10(abs(TFDNXjanelaRetangular2)), 'k');
title(['FFT com M = 2N pontos']);
ylabel('Xjanela Retangular(f)| dB');
grid on; 

% Janela de Hamming
figure(2);
subplot(221);
plot(ff1, abs(TFDNXjanelaHamming1),'k');
title(['FFT com M = N pontos']);
ylabel('|Xjanela de Hamming(f)|');
grid on;
%
subplot(223);
plot(ff1,20*log10(abs(TFDNXjanelaHamming1)),'k'); 
title(['FFT com M = N pontos']);
ylabel('|Xjanela de Hamming(f)| dB');
grid on; 
%
subplot(222)
plot(ff2, abs(TFDNXjanelaHamming2), 'k');
title(['FFT com M = 2N pontos']);
ylabel('|Xjanela de Hamming(f))|');
grid on;
%
subplot(224)
plot(ff2, 20*log10(abs(TFDNXjanelaHamming2)), 'k');
title(['FFT com M = 2N pontos']);
ylabel('Xjanela de Hamming(f)| dB');
grid on;

% Janela de Blackman
figure(3); 
subplot(221);
plot(ff1, abs(TFDNXjanelaBlackman1),'k');
title(['FFT com M = N pontos']);
ylabel('|Xjanela de Blackman(f)|');
grid on;
%
subplot(223);
plot(ff1,20*log10(abs(TFDNXjanelaBlackman1)),'k'); 
title(['FFT com M = N pontos']);
ylabel('|Xjanela de Blackman(f)| dB');
grid on; 
%
subplot(222)
plot(ff2, abs(TFDNXjanelaBlackman2), 'k');
title(['FFT com M = 2N pontos']);
ylabel('|Xjanela de Blackman(f))|');
grid on;
%
subplot(224)
plot(ff2, 20*log10(abs(TFDNXjanelaBlackman2)), 'k');
title(['FFT com M = 2N pontos']);
ylabel('Xjanela de Blackman(f)| dB');
grid on;








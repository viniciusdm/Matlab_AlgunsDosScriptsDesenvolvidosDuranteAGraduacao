%////////////////////////////////////////////////////////////////////////
%///        Escola Polit�cnica da Universidade de S�o Paulo           ///
%/// EN3 - PTC3432 - 04/07/2021 - Interpreta��o Geom�trica de Sinais  ///
%///         Nome: Vinicius Bueno de Moraes, NUSP: 10256432           ///   
%///                                                                  ///
%///             Simula��o para obten��o da taxa de erro              ///
%////////////////////////////////////////////////////////////////////////

%Inicializa��es.
SNR=[3  : 14];                           %Valores de SNR a serem simulados (dB).
EbN0=[1 : 11];                           %Valores de Eb/N0 a serem simulados(dB).
a=[-3 -1 1 1i*sqrt(3)];                  %Define a constela��o no vetor a.
N=10^5;                                  %N�mero de s�mbolos a serem transmitidos e recebidos -
                                         %Considerada menor probabilidade
                                         %de erro para os Valores de SNR e
                                         %Eb/N0 * 1^2 = 10^6 (aprox.).
                 
dim=2;                                   %N�mero de dimens�es.                          
                                                                        
erros_SNR=zeros(1,length(SNR));          %Vetor que cont�m o n�mero de s�mbolos errados.
erros_EbN0=zeros(1,length(EbN0));        %Vetor que cont�m o n�mero de s�mbolos errados.

%Gera��o de �ndices 1, 2, 3 ou 4 aleat�rios e equiprov�veis para escolher
%os s�mbolos transmitidos em cada instante de tempo kT.
ponteiro=randi(4,1,N);
ponteiro2=randi(4, 1, N);

%Gera o vetor x j� considerando o modelo discreto.
x=a(ponteiro);  
w=a(ponteiro2);

%C�culo da pot�ncia m�dia dos s�mbolos (supondo Eq=1).
Psimbolos=(1/length(a))*(abs(a(1)))^2 + (1/length(a))*(abs(a(2)))^2 + (1/length(a))*(abs(a(3)))^2 + (1/length(a))*(abs(a(4)))^2; 

%C�lculo da energia m�dia de bits.
Eb=(1/(log2(length(a))))*Psimbolos; 

%C�lculo da taxa SNR. 
for m=1:length(SNR)
    N0=(2*Psimbolos)/(dim*(10^(SNR(m)/10)));                                           %C�lculo de N0.
    y=x+sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));                                         %C�lculo da sinal recebido ap�s adi��o de ruido. 
                                                                                     
    %Compara, para cada instante k, o simbolo recebido e o simbolo transmitido.
    for k=1:N
        rx = x(k) * ones(1, length(a))
        ry = y(k) * ones(1, length(a))
        d_sq1=abs(rx-a).^2;
        d_sq2=abs(ry-a).^2;
        [dmin1,ind1]=min(d_sq1);
        [dmin2, ind2]=min(d_sq2);
        if ind1~=ind2                                                                   %if SIMBOLO_ESTIMADO~=x(k)  %Se o seu SIMBOLO_ESTIMADO no instante k � diferente do transmitido x(k) ent�o h� um erro. � poss�vel verificar tamb�m pelo �ndice estimado em rela��o ao �ndice transmitido, similar no programa que tra�a as regi�es de decis�o. 
            erros_SNR(m)=erros_SNR(m)+1;                                                %Incrementa de 1 o erro para a m-�sima SNR (explicitado no estudo te�rico).
        end
    end
end

%Comparativo - SNR.
Taxa_erros_SNR=erros_SNR/N;                                                             %Obt�m a taxa de erro.
Pe_SNR=2*qfunc(sqrt((4*(10.^(SNR/10))/7)))+(1/2)*qfunc(sqrt((12*(10.^(SNR/10))/7)));    %Probabildade de erro calculada a partir do limitante da uni�o a relacionado com a SNR.
figure
semilogy(SNR,Pe_SNR,'bv-');                                                             %Plota o limitante da uni�o.
hold on;
semilogy(SNR,Taxa_erros_SNR,'ko-')                                                      %Plota a taxa obtida na simula��o (SNR).
title('Comparativo de Desempenho - Simulado vs Te�rico (SNR)');
xlabel('SNR (dB)');  
ylabel('Pe - Probabilidade de Erro');  
legend('Prob. de erro de simb.','Taxa de erro de simb.');
grid on;

%C�lculo da taxa Eb/N0. 
for m=1:length(EbN0)                                                                    
    N0=(2*Psimbolos)/(dim*2*(10^(EbN0(m)/10)));                                        %C�lculo de N0.
    y=x+sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));                                         %C�lculo da sinal recebido ap�s adi��o de ruido (explicitado no estudo te�rico).  
        
    %Compara, para cada instante k, o simbolo recebido e o simbolo transmitido.
    for k=1:N
        rx = x(k) * ones(1, length(a))
        ry = y(k) * ones(1, length(a))
        d_sq1=abs(rx-a).^2;
        d_sq2=abs(ry-a).^2;
        [dmin1,ind1]=min(d_sq1);
        [dmin2, ind2]=min(d_sq2);
        if ind1~=ind2                                                                  %Se o SIMBOLO_ESTIMADO no instante k � diferente do transmitido x(k) ent�o h� um erro. 
            erros_EbN0(m)=erros_EbN0(m)+1;                                             %Incrementa de 1 o erro para a m-�sima Eb/N0.
        end
    end
end   

%Comparativo - Eb/N0.
Taxa_erros_EbN0=erros_EbN0/N;                                                          %Obt�m a taxa de erro.
Pe_EbN0=2*qfunc(sqrt((8/7)*(10.^(EbN0/10))))+(1/2)*qfunc(sqrt((24/7)*(10.^(EbN0/10))));%Probabildade de erro calculado a partir do limitante da uni�o a relacionando com Eb/N0.
figure
semilogy(EbN0,Pe_EbN0,'bv-');                                                          %Plota o limitante da uni�o.
hold on;
semilogy(EbN0,Taxa_erros_EbN0,'ko-');                                                  %Plota a taxa obtida na simula��o (Eb/N0).
title('Comparativo de Desempenho - Simulado vs Te�rico (Eb/N0)');  
xlabel('Eb/N0 (dB)');     
ylabel('Pe - Probabilidade de Erro');  
legend('Prob. de erro de simb.','Taxa de erro de simb.');
grid on;

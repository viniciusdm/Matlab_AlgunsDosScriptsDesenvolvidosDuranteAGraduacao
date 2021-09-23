%////////////////////////////////////////////////////////////////////////
%///        Escola Politécnica da Universidade de São Paulo           ///
%/// EN3 - PTC3432 - 04/07/2021 - Interpretação Geométrica de Sinais  ///
%///         Nome: Vinicius Bueno de Moraes, NUSP: 10256432           ///   
%///                                                                  ///
%///             Simulação para obtenção da taxa de erro              ///
%////////////////////////////////////////////////////////////////////////

%Inicializações.
SNR=[3  : 14];                           %Valores de SNR a serem simulados (dB).
EbN0=[1 : 11];                           %Valores de Eb/N0 a serem simulados(dB).
a=[-3 -1 1 1i*sqrt(3)];                  %Define a constelação no vetor a.
N=10^5;                                  %Número de símbolos a serem transmitidos e recebidos -
                                         %Considerada menor probabilidade
                                         %de erro para os Valores de SNR e
                                         %Eb/N0 * 1^2 = 10^6 (aprox.).
                 
dim=2;                                   %Número de dimensões.                          
                                                                        
erros_SNR=zeros(1,length(SNR));          %Vetor que contém o número de símbolos errados.
erros_EbN0=zeros(1,length(EbN0));        %Vetor que contém o número de símbolos errados.

%Geração de índices 1, 2, 3 ou 4 aleatórios e equiprováveis para escolher
%os símbolos transmitidos em cada instante de tempo kT.
ponteiro=randi(4,1,N);
ponteiro2=randi(4, 1, N);

%Gera o vetor x já considerando o modelo discreto.
x=a(ponteiro);  
w=a(ponteiro2);

%Cáculo da potência média dos símbolos (supondo Eq=1).
Psimbolos=(1/length(a))*(abs(a(1)))^2 + (1/length(a))*(abs(a(2)))^2 + (1/length(a))*(abs(a(3)))^2 + (1/length(a))*(abs(a(4)))^2; 

%Cálculo da energia média de bits.
Eb=(1/(log2(length(a))))*Psimbolos; 

%Cálculo da taxa SNR. 
for m=1:length(SNR)
    N0=(2*Psimbolos)/(dim*(10^(SNR(m)/10)));                                           %Cálculo de N0.
    y=x+sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));                                         %Cálculo da sinal recebido após adição de ruido. 
                                                                                     
    %Compara, para cada instante k, o simbolo recebido e o simbolo transmitido.
    for k=1:N
        rx = x(k) * ones(1, length(a))
        ry = y(k) * ones(1, length(a))
        d_sq1=abs(rx-a).^2;
        d_sq2=abs(ry-a).^2;
        [dmin1,ind1]=min(d_sq1);
        [dmin2, ind2]=min(d_sq2);
        if ind1~=ind2                                                                   %if SIMBOLO_ESTIMADO~=x(k)  %Se o seu SIMBOLO_ESTIMADO no instante k é diferente do transmitido x(k) então há um erro. É possível verificar também pelo índice estimado em relação ao índice transmitido, similar no programa que traça as regiões de decisão. 
            erros_SNR(m)=erros_SNR(m)+1;                                                %Incrementa de 1 o erro para a m-ésima SNR (explicitado no estudo teórico).
        end
    end
end

%Comparativo - SNR.
Taxa_erros_SNR=erros_SNR/N;                                                             %Obtém a taxa de erro.
Pe_SNR=2*qfunc(sqrt((4*(10.^(SNR/10))/7)))+(1/2)*qfunc(sqrt((12*(10.^(SNR/10))/7)));    %Probabildade de erro calculada a partir do limitante da união a relacionado com a SNR.
figure
semilogy(SNR,Pe_SNR,'bv-');                                                             %Plota o limitante da união.
hold on;
semilogy(SNR,Taxa_erros_SNR,'ko-')                                                      %Plota a taxa obtida na simulação (SNR).
title('Comparativo de Desempenho - Simulado vs Teórico (SNR)');
xlabel('SNR (dB)');  
ylabel('Pe - Probabilidade de Erro');  
legend('Prob. de erro de simb.','Taxa de erro de simb.');
grid on;

%Cálculo da taxa Eb/N0. 
for m=1:length(EbN0)                                                                    
    N0=(2*Psimbolos)/(dim*2*(10^(EbN0(m)/10)));                                        %Cálculo de N0.
    y=x+sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));                                         %Cálculo da sinal recebido após adição de ruido (explicitado no estudo teórico).  
        
    %Compara, para cada instante k, o simbolo recebido e o simbolo transmitido.
    for k=1:N
        rx = x(k) * ones(1, length(a))
        ry = y(k) * ones(1, length(a))
        d_sq1=abs(rx-a).^2;
        d_sq2=abs(ry-a).^2;
        [dmin1,ind1]=min(d_sq1);
        [dmin2, ind2]=min(d_sq2);
        if ind1~=ind2                                                                  %Se o SIMBOLO_ESTIMADO no instante k é diferente do transmitido x(k) então há um erro. 
            erros_EbN0(m)=erros_EbN0(m)+1;                                             %Incrementa de 1 o erro para a m-ésima Eb/N0.
        end
    end
end   

%Comparativo - Eb/N0.
Taxa_erros_EbN0=erros_EbN0/N;                                                          %Obtém a taxa de erro.
Pe_EbN0=2*qfunc(sqrt((8/7)*(10.^(EbN0/10))))+(1/2)*qfunc(sqrt((24/7)*(10.^(EbN0/10))));%Probabildade de erro calculado a partir do limitante da união a relacionando com Eb/N0.
figure
semilogy(EbN0,Pe_EbN0,'bv-');                                                          %Plota o limitante da união.
hold on;
semilogy(EbN0,Taxa_erros_EbN0,'ko-');                                                  %Plota a taxa obtida na simulação (Eb/N0).
title('Comparativo de Desempenho - Simulado vs Teórico (Eb/N0)');  
xlabel('Eb/N0 (dB)');     
ylabel('Pe - Probabilidade de Erro');  
legend('Prob. de erro de simb.','Taxa de erro de simb.');
grid on;

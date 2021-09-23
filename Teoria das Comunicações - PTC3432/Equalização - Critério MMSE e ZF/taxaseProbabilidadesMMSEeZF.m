%////////////////////////////////////////////////////////////////////////
%///        Escola Politécnica da Universidade de São Paulo           ///
%///            EN4 - PTC3432 - 14/07/2021 - Equalização              ///
%///         Nome: Vinicius Bueno de Moraes, NUSP: 10256432           ///   
%///                                                                  ///
%///       Item B) Taxas e Probabilidades de Erro - MMSE e ZF         ///
%////////////////////////////////////////////////////////////////////////

Nsimbolos = 10^6;                                 %Número de símbolos a serem simulados - critério de escola explanado no documento teórico
a = sign(randn(Nsimbolos,1));                     %Geração dos símbolos 2-PAM
sigma_a2 = 1;                                     %Definicão de sigma_a2 - Explanado no documento teórico
SNR = 3:2:19;                                     %Faixa de SNR simulada
d = 7;                                            %Atraso do sinal desejado, veja que o valor máximo possível é o de N
h=[0.8 1].';                                      %Define Canal
M = length(h);                                    
N = 9;                                            %Número de coeficientes do FIR do equalizador

sigma_n2 = [];                                    %Converte as SNRs para linear e define sigma_n2
for i = SNR
    linear = (norm(h).^2) / (10.^(i/10));
    sigma_n2 = [sigma_n2 linear];
end

BER = zeros(1,length(SNR));                       %Vetor para guardar a taxa de erro
BER_zf = zeros(1,length(SNR));                    %Vetor para guardar a taxa de erro - caso ZF (para poder comparar em paralelo)

x=conv(h,a);                                      %Saída filtro casado (sem ruído)
Hc=convmtx(h.',N);                                %Define a matriz de convolução 

%Calculo da BER para cada valor de SNR

for k = 1:length(SNR)
    u = x+sqrt(sigma_n2(k))*randn(length(x),1);                                                                              %Saída do canal
    R = sigma_a2 * Hc * Hc' + sigma_n2(k) * eye(N);                                                                          %R
    R_zf = sigma_a2 * Hc * Hc';                                                                                              %R_zf (considerando o ruído nulo)
    
    pn = Hc * [zeros(d, 1); sigma_a2; zeros(N-d, 1)];                                                                        %Define p de acordo com o atraso desejado
    pnt = transpose(pn);
    
    w = inv(R) * pn;                                                                                                         %Obtém os coeficientes do filtro de acordo com o critério - MMSE
    y = conv(w',u);                                                                                                          %Obtém a saída do equalizador - MMSE
    BER(k) = sum(.5*abs(sign(y([(N+M):(Nsimbolos-(N+M-1))]+d))-a([(N+M):(Nsimbolos-(N+M-1))])))/(Nsimbolos-2*(M+N-1));       %Calcular a BER - MMSE
    
    w_zf = inv(R_zf) * pn;                                                                                                   %Obtém os coeficientes do filtro de acordo com o critério - ZF (MMSE com ruido zero)
    y_zf = conv(w_zf', u);                                                                                                   %Obtém a saída do equalizador - ZF (MMSE com ruido zero)
    BER_zf(k) = sum(.5*abs(sign(y_zf([(N+M):(Nsimbolos-(N+M-1))]+d))-a([(N+M):(Nsimbolos-(N+M-1))])))/(Nsimbolos-2*(M+N-1)); %Calcula a BER - ZF (MMSE com ruido zero
end

%Cálculo da Probabilidade de erro sem equalizador e a Probabilidade de erro por MMSE em função de certo atraso                                                     
Pe_apx = (1/2) * qfunc((0.2*(sqrt((10.^(SNR/10))/(norm(h)^2)))));                                                            %Probabilidade de erro sem equalizador - aproximada (limitante da união)
Pe_real = (1/2) * qfunc((1.8*(sqrt((10.^(SNR/10))/(norm(h)^2))))) + (1/2) * qfunc((0.2*(sqrt((10.^(SNR/10))/(norm(h)^2))))); %Plota a Probabilidade de erro sem equalizador - real (limitante da união)
h_eq = conv(h, w');                                                                                                          %Define h_qe - como visto na teoria e explanado no documento teórico
sigma_n = sqrt((norm(h_eq)^2 - norm(h_eq(d+1))^2) * sigma_a2 + norm(w)^2 * sigma_n2);                                        %Sigma_n - Probabilidade de Erro - Fórmula vista em aula e explanada no documento teórico
Pe_mmse = qfunc(abs(h_eq(d+1))./sigma_n);                                                                                    %Probabildade de erro - MMSE - Deduzida e explanada no documento teórico

%Plots Comparativos 
figure;
semilogy(SNR, Pe_apx, 'bo-');                                             %Plota a Probabilidade de erro sem equalizador - aproximada (limitante da união)
hold on;
semilogy(SNR, Pe_real, 'yo-');                                            %Plota a Probabilidade de erro sem equalizador - real (limitante da união)
hold on;
semilogy(SNR, Pe_mmse, 'ro-');                                            %Plota a Probabilidade - MMSE
hold on;
semilogy(SNR, BER,'gx-');                                                 %Plota a Taxa de erro obtida - MMSE
hold on;
semilogy(SNR, BER_zf,'kx-');                                              %Plota a Taxa de erro - ZF (MMSE com ruido zero)
hold on; 
title('Comparativo de desempenho para diferentes Métodos');
suptitle('Probabilidade de erro vs Taxa de erro');
xlabel('SNR (dB)');  
ylabel('Pb erro | Tx erro');  
legend('Prob. de erro de simb. - Sem Equa. - Apróx','Prob. de erro de simb. - Sem Equa. - Real','Prob. de erro de simb. - MMSE','Taxa de erro de simb. - MMSE','Taxa de erro de simb. - ZF');
fprintf('\n EN4 - Item B - PTC3432');
fprintf('\n\n Executado com Sucesso');
grid on;

figure;
zplane(h_eq');
title('Diagrama de Póloz e Zeros do Canal Equivalente');
xlabel('Re');    
ylabel('Im');
grid on;

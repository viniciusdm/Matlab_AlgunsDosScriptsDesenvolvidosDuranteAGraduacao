%////////////////////////////////////////////////////////////////////////
%///        Escola Polit�cnica da Universidade de S�o Paulo           ///
%///            EN4 - PTC3432 - 14/07/2021 - Equaliza��o              ///
%///         Nome: Vinicius Bueno de Moraes, NUSP: 10256432           ///   
%///                                                                  ///
%///       Item B) Taxas e Probabilidades de Erro - MMSE e ZF         ///
%////////////////////////////////////////////////////////////////////////

Nsimbolos = 10^6;                                 %N�mero de s�mbolos a serem simulados - crit�rio de escola explanado no documento te�rico
a = sign(randn(Nsimbolos,1));                     %Gera��o dos s�mbolos 2-PAM
sigma_a2 = 1;                                     %Definic�o de sigma_a2 - Explanado no documento te�rico
SNR = 3:2:19;                                     %Faixa de SNR simulada
d = 7;                                            %Atraso do sinal desejado, veja que o valor m�ximo poss�vel � o de N
h=[0.8 1].';                                      %Define Canal
M = length(h);                                    
N = 9;                                            %N�mero de coeficientes do FIR do equalizador

sigma_n2 = [];                                    %Converte as SNRs para linear e define sigma_n2
for i = SNR
    linear = (norm(h).^2) / (10.^(i/10));
    sigma_n2 = [sigma_n2 linear];
end

BER = zeros(1,length(SNR));                       %Vetor para guardar a taxa de erro
BER_zf = zeros(1,length(SNR));                    %Vetor para guardar a taxa de erro - caso ZF (para poder comparar em paralelo)

x=conv(h,a);                                      %Sa�da filtro casado (sem ru�do)
Hc=convmtx(h.',N);                                %Define a matriz de convolu��o 

%Calculo da BER para cada valor de SNR

for k = 1:length(SNR)
    u = x+sqrt(sigma_n2(k))*randn(length(x),1);                                                                              %Sa�da do canal
    R = sigma_a2 * Hc * Hc' + sigma_n2(k) * eye(N);                                                                          %R
    R_zf = sigma_a2 * Hc * Hc';                                                                                              %R_zf (considerando o ru�do nulo)
    
    pn = Hc * [zeros(d, 1); sigma_a2; zeros(N-d, 1)];                                                                        %Define p de acordo com o atraso desejado
    pnt = transpose(pn);
    
    w = inv(R) * pn;                                                                                                         %Obt�m os coeficientes do filtro de acordo com o crit�rio - MMSE
    y = conv(w',u);                                                                                                          %Obt�m a sa�da do equalizador - MMSE
    BER(k) = sum(.5*abs(sign(y([(N+M):(Nsimbolos-(N+M-1))]+d))-a([(N+M):(Nsimbolos-(N+M-1))])))/(Nsimbolos-2*(M+N-1));       %Calcular a BER - MMSE
    
    w_zf = inv(R_zf) * pn;                                                                                                   %Obt�m os coeficientes do filtro de acordo com o crit�rio - ZF (MMSE com ruido zero)
    y_zf = conv(w_zf', u);                                                                                                   %Obt�m a sa�da do equalizador - ZF (MMSE com ruido zero)
    BER_zf(k) = sum(.5*abs(sign(y_zf([(N+M):(Nsimbolos-(N+M-1))]+d))-a([(N+M):(Nsimbolos-(N+M-1))])))/(Nsimbolos-2*(M+N-1)); %Calcula a BER - ZF (MMSE com ruido zero
end

%C�lculo da Probabilidade de erro sem equalizador e a Probabilidade de erro por MMSE em fun��o de certo atraso                                                     
Pe_apx = (1/2) * qfunc((0.2*(sqrt((10.^(SNR/10))/(norm(h)^2)))));                                                            %Probabilidade de erro sem equalizador - aproximada (limitante da uni�o)
Pe_real = (1/2) * qfunc((1.8*(sqrt((10.^(SNR/10))/(norm(h)^2))))) + (1/2) * qfunc((0.2*(sqrt((10.^(SNR/10))/(norm(h)^2))))); %Plota a Probabilidade de erro sem equalizador - real (limitante da uni�o)
h_eq = conv(h, w');                                                                                                          %Define h_qe - como visto na teoria e explanado no documento te�rico
sigma_n = sqrt((norm(h_eq)^2 - norm(h_eq(d+1))^2) * sigma_a2 + norm(w)^2 * sigma_n2);                                        %Sigma_n - Probabilidade de Erro - F�rmula vista em aula e explanada no documento te�rico
Pe_mmse = qfunc(abs(h_eq(d+1))./sigma_n);                                                                                    %Probabildade de erro - MMSE - Deduzida e explanada no documento te�rico

%Plots Comparativos 
figure;
semilogy(SNR, Pe_apx, 'bo-');                                             %Plota a Probabilidade de erro sem equalizador - aproximada (limitante da uni�o)
hold on;
semilogy(SNR, Pe_real, 'yo-');                                            %Plota a Probabilidade de erro sem equalizador - real (limitante da uni�o)
hold on;
semilogy(SNR, Pe_mmse, 'ro-');                                            %Plota a Probabilidade - MMSE
hold on;
semilogy(SNR, BER,'gx-');                                                 %Plota a Taxa de erro obtida - MMSE
hold on;
semilogy(SNR, BER_zf,'kx-');                                              %Plota a Taxa de erro - ZF (MMSE com ruido zero)
hold on; 
title('Comparativo de desempenho para diferentes M�todos');
suptitle('Probabilidade de erro vs Taxa de erro');
xlabel('SNR (dB)');  
ylabel('Pb erro | Tx erro');  
legend('Prob. de erro de simb. - Sem Equa. - Apr�x','Prob. de erro de simb. - Sem Equa. - Real','Prob. de erro de simb. - MMSE','Taxa de erro de simb. - MMSE','Taxa de erro de simb. - ZF');
fprintf('\n EN4 - Item B - PTC3432');
fprintf('\n\n Executado com Sucesso');
grid on;

figure;
zplane(h_eq');
title('Diagrama de P�loz e Zeros do Canal Equivalente');
xlabel('Re');    
ylabel('Im');
grid on;

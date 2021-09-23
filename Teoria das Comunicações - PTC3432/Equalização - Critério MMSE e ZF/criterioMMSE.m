%////////////////////////////////////////////////////////////////////////
%///        Escola Politécnica da Universidade de São Paulo           ///
%///            EN4 - PTC3432 - 14/07/2021 - Equalização              ///
%///         Nome: Vinicius Bueno de Moraes, NUSP: 10256432           ///   
%///                                                                  ///
%/// Item A) Obtenção de Tabela MSE - Critério MMSE, p/ todos atrasos ///
%////////////////////////////////////////////////////////////////////////

h = [0.8 1].';                     %Coeficientes do canal FIR s ser utilizado (definido como vetor coluna por conveniência na área de comunicações)
N = 9;                             %Número de coeficientes de do equalizador linear FIR
sigma_a2 = 1;                      %Variância simbólica calculada no documento teórico complemetar a este
SNR = [4 15 20];                   %Valores de SNR a serem considerados

Hc=convmtx(h.',N);                 %define a matriz de convolução 

sigma_n2 = [];                     %Converte as SNRs para linear e define sigma_n2
for i = SNR
    linear = (norm(h).^2) / (10.^(i/10));
    sigma_n2 = [sigma_n2 linear];
end

mp = [];
mw = [];
mc = [];
mm = [];
f = 0; 
for i = sigma_n2                   %Laço responsável pelos cálculos
    R = sigma_a2 * Hc * Hc' + i * eye(N);
    f = f + 1;
    for j = 1 : N+1                %p0,p1, ..., pn
        p = Hc * [zeros(j-1, 1); sigma_a2; zeros(length(h)+N-(j+1), 1)];
        mp = [mp p];
        j = j + 1;
        if j == N+2
            for k = 1+((N+1)*(f-1)) : (N+1)*f %w0, w1, ..., wn
               w = inv(R) * mp(:, k);
               mw = [mw w];
               k = k + 1;
               if k == (N+2)+((f-1)*(N+1))    
                   for l = 1+((N+1)*(f-1)) : (N+1)*f   %Convolução com o filtro (extra)
                       c = conv(h, mw(:, l));
                       mc = [mc c];
                       l = l + 1;
                   end
                   for n = 1+((N+1)*(f-1)) : (N+1)*f   %Cálculo final (MSE - Critério MMSE)
                       m = sigma_a2 - mw(:, n)' * mp(:, n);
                       mm = [mm m];
                       n = n + 1;
                   end
               end
            end
       end
    end
end

%Exibição de resultados
mse = reshape(mm, N+1, [])';
[val, loc] = min(mse');
fprintf('\n EN4 - Item A - PTC3432')
fprintf('\n\n Valores (Tabela MSE - Critério MMSE) para SNRs iguais a: ')
disp(SNR) 
disp(mse)
menores_valores = (val)'
atrasos_otimos_d = (loc-1)'

%Plot do filtro (para avaliar localização de seus polos e zeros se for de
%interesse) 

figure;
subplot(1,2,1);
zplane(h');
title('Diagrama de Pólos e Zeros - Canal');
xlabel('Re(Z)');     
ylabel('Im(Z)');
grid on;

%Plot da variação do atraso ótimo com a SNR

subplot(1,2,2);
stem(SNR, atrasos_otimos_d, 'rs-');
title('Variação da Posição do Atraso Ótimo');
xlabel('SNR(dB)');    
ylabel('Atraso Ótimo');
grid on;


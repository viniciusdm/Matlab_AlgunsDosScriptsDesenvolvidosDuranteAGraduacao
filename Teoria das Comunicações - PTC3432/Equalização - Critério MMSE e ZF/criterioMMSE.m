%////////////////////////////////////////////////////////////////////////
%///        Escola Polit�cnica da Universidade de S�o Paulo           ///
%///            EN4 - PTC3432 - 14/07/2021 - Equaliza��o              ///
%///         Nome: Vinicius Bueno de Moraes, NUSP: 10256432           ///   
%///                                                                  ///
%/// Item A) Obten��o de Tabela MSE - Crit�rio MMSE, p/ todos atrasos ///
%////////////////////////////////////////////////////////////////////////

h = [0.8 1].';                     %Coeficientes do canal FIR s ser utilizado (definido como vetor coluna por conveni�ncia na �rea de comunica��es)
N = 9;                             %N�mero de coeficientes de do equalizador linear FIR
sigma_a2 = 1;                      %Vari�ncia simb�lica calculada no documento te�rico complemetar a este
SNR = [4 15 20];                   %Valores de SNR a serem considerados

Hc=convmtx(h.',N);                 %define a matriz de convolu��o 

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
for i = sigma_n2                   %La�o respons�vel pelos c�lculos
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
                   for l = 1+((N+1)*(f-1)) : (N+1)*f   %Convolu��o com o filtro (extra)
                       c = conv(h, mw(:, l));
                       mc = [mc c];
                       l = l + 1;
                   end
                   for n = 1+((N+1)*(f-1)) : (N+1)*f   %C�lculo final (MSE - Crit�rio MMSE)
                       m = sigma_a2 - mw(:, n)' * mp(:, n);
                       mm = [mm m];
                       n = n + 1;
                   end
               end
            end
       end
    end
end

%Exibi��o de resultados
mse = reshape(mm, N+1, [])';
[val, loc] = min(mse');
fprintf('\n EN4 - Item A - PTC3432')
fprintf('\n\n Valores (Tabela MSE - Crit�rio MMSE) para SNRs iguais a: ')
disp(SNR) 
disp(mse)
menores_valores = (val)'
atrasos_otimos_d = (loc-1)'

%Plot do filtro (para avaliar localiza��o de seus polos e zeros se for de
%interesse) 

figure;
subplot(1,2,1);
zplane(h');
title('Diagrama de P�los e Zeros - Canal');
xlabel('Re(Z)');     
ylabel('Im(Z)');
grid on;

%Plot da varia��o do atraso �timo com a SNR

subplot(1,2,2);
stem(SNR, atrasos_otimos_d, 'rs-');
title('Varia��o da Posi��o do Atraso �timo');
xlabel('SNR(dB)');    
ylabel('Atraso �timo');
grid on;


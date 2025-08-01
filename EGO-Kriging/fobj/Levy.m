function f = Levy(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exemplo de função n-dimensional para o otimizador EGO-Kriging.
%
%   Entrada:
%       d - vetor de variáveis de projeto.
%
%   Saída:
%       f - vetor de valores da função (cálculo vetorizado)
%
%   Autor: Fábio Felipe dos Santos Nascentes
%   Orientador: Rafael Holdorf Lopez
%   Universidade Federal de Santa Catarina
%   Data Última Atualização: 05/11/2018
%
%   Copyright (c) 2018 Fábio Nascentes
%   Todos os direitos autorais reservados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A dimensão da função Levy será definida de acordo com a dimensão do vetor x (ou
% número de colunas da matriz x)
d = size(x, 2);

% Cria o vetor que armazenará os valores da função
f = zeros(size(x, 1), 1);

% Cálculo estocástico da função
for i = 1:size(x, 1)
    w = zeros(1, d);
    
    for j = 1:d
        w(j) = 1 + (x(i, j) - 1)/4;
    end
    
    termo1 = (sin(pi*w(1)))^2;
    
    termo3 = (w(d)-1)^2 * (1+(sin(2*pi*w(d)))^2);
    
    soma = 0;
    for j = 1:(d-1)
        wi = w(j);
        novo = (wi-1)^2 * (1+10*(sin(pi*wi+1))^2);
        soma = soma + novo;
    end
    
    f(i) = (termo1 + soma + termo3);
end
end
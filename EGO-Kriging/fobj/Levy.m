function f = Levy(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exemplo de fun��o n-dimensional para o otimizador EGO-Kriging.
%
%   Entrada:
%       d - vetor de vari�veis de projeto.
%
%   Sa�da:
%       f - vetor de valores da fun��o (c�lculo vetorizado)
%
%   Autor: F�bio Felipe dos Santos Nascentes
%   Orientador: Rafael Holdorf Lopez
%   Universidade Federal de Santa Catarina
%   Data �ltima Atualiza��o: 05/11/2018
%
%   Copyright (c) 2018 F�bio Nascentes
%   Todos os direitos autorais reservados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A dimens�o da fun��o Levy ser� definida de acordo com a dimens�o do vetor x (ou
% n�mero de colunas da matriz x)
d = size(x, 2);

% Cria o vetor que armazenar� os valores da fun��o
f = zeros(size(x, 1), 1);

% C�lculo estoc�stico da fun��o
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
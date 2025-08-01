function f = Fabio1D(d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exemplo de função 1D para o otimizador EGO-Kriging.
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
a1 = (d.^3) ./ 6;
a2 = (17 / 16) .* (d.^2);
a3 = (7 / 12) .* d;
a4 = 37 / 16;
A = a1 - a2 + a3 + a4;

b1 = -1 / 200;
b2 = d .^ 4;
b3 = 8 .* (d.^3);
b4 = 14 .* (d.^2);
b5 = 8 .* d;
b6 = 15;
B = b1 .* (b2 - b3 + b4 + b5 - b6);

c1 = ((2/5) * (d.^2)) + ((5/2) * d);
C = sin(c1);

f = (-abs(A .* exp(B))) .* C;
end
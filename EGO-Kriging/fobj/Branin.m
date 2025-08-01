function f = Branin(d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exemplo de função 2D para o otimizador EGO-Kriging.
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
% Constantes que são utilizadas pela função Branin
a1 = 1;
a2 = 5.1 / (4 * (pi^2));
a3 = 5 / pi;
a4 = 6;
a5 = 10;
a6 = 1 / (8 * pi);

% Cálculo da função
f = a1*(d(:,2) - (a2*(d(:,1).^2)) + a3*d(:,1) - a4).^2 + a5*((1 - a6)*cos(d(:,1)))...
     + a5 + 5*d(:,1);
end
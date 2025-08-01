function f = Branin(d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exemplo de fun��o 2D para o otimizador EGO-Kriging.
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
% Constantes que s�o utilizadas pela fun��o Branin
a1 = 1;
a2 = 5.1 / (4 * (pi^2));
a3 = 5 / pi;
a4 = 6;
a5 = 10;
a6 = 1 / (8 * pi);

% C�lculo da fun��o
f = a1*(d(:,2) - (a2*(d(:,1).^2)) + a3*d(:,1) - a4).^2 + a5*((1 - a6)*cos(d(:,1)))...
     + a5 + 5*d(:,1);
end
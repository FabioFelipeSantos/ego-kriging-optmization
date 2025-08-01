function f = Hartman(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exemplo de função 3D ou 6D para o otimizador EGO-Kriging.
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
% Escolhe com base no número da dimensão k, se a função avaliada é a Hartman 3D ou
% Hartman 6D.
if size(x, 2) == 3
    f = Hartman3D(x);
else
    f = Hartman6D(x);
end
end

function f = Hartman3D(x)
% Cálculo da função Hartman 3D

% Parâmetros da função
a = [1.0 1.2 3.0 3.2];

B = [3.0 10.0 30.0;
     0.0 10.0 35.0;
     3.0 10.0 30.0;
     0.1 10.0 35.0];
 
D = [0.3689  0.1170 0.2673;
     0.4699  0.4387 0.7470;
     0.1091  0.8732 0.5547;
     0.03815 0.5743 0.8828];

% Definição das somas referente ao valor da função e o cálculo efetivo da função
soma1 = 0;
soma2 = 0;
f = zeros(size(x,1), 1);

% Determina a média da distribuição normal
mu = 1;

% Determina o desvio padrão da distribução normal
dp = 0;

% Cálculo da variável estocástica que causa o erro ou ruído na saída da função
X = mu + dp * randn([1, size(x, 2)]);

% Cálculo estocástico da função
for n = 1:size(x,1)
    for i = 1:4
        for j = 1:3
            soma1 = soma1 + B(i,j) * ((x(n, j)*X(j) - D(i,j))^2);
        end
        soma2 = soma2 + a(i) * exp(-soma1);
        soma1 = 0;
    end
    f(n, 1) = -soma2;
    soma2 = 0;
end
end

function f = Hartman6D(x)
% Cálculo da função Hartman 3D

% Parâmetros da função
a = [1.0 1.2 3.0 3.2];

B = [10.0 3.0  17.0 3.5  1.7  8.0
     0.05 10.0 17.0 0.1  8.0  14.0
     3.0  3.5  1.7  10.0 17.0 8.0
     17.0 8.0  0.05 10.0 0.1  14.0];
D = [0.1312 0.1696 0.5569 0.0124 0.8283 0.5886
     0.2329 0.4135 0.8307 0.3736 0.1004 0.9991
     0.2348 0.1451 0.3522 0.2883 0.3047 0.6650
     0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];

% Definição das somas referente ao valor da função e o cálculo efetivo da função
soma1 = 0;
soma2 = 0;
f = zeros(size(x,1), 1);

% Cálculo estocástico da função
for n = 1:size(x,1)
    for i = 1:4
        for j = 1:6
            soma1 = soma1 + B(i,j) * ((x(n, j) - D(i,j))^2);
        end
        soma2 = soma2 + a(i) * exp(-soma1);
        soma1 = 0;
    end
    f(n, 1) = -soma2;
    soma2 = 0;
end
end
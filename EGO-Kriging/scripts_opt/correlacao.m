function r = correlacao(theta, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Correlação Gaussiana para o ajuste do Kriging
%
%   Entrada:
%       theta - parâmetros da função de correlação;
%       d - matriz (n x k) com as diferenças entre os pontos amostrais.
%
%   Saída:
%       r - correlação desses pontos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula o número de diferenças e a dimensão do problema
n = size(d, 1);

% Ajusta um vetor com o parâmetro theta para cada distância
auxTheta = repmat(-theta(:)', n, 1);

% Calcula o produto da distância pelo parâmetro theta
auxFinal = (d.^2) .* auxTheta;

% Calcula a correlação
r = exp(sum(auxFinal, 2));
return
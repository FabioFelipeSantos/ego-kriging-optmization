function [yHat, predVar] = predKrg(x, Modelo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Rotina de predição utilizando o metamodelo em Kriging.
%
%   Entrada:
%       x - vetor de k variáveis que será predita pelo metamodelo;
%       Modelo - informações sobre o modelo.
%
%   Saída:
%       yHat - valor escalar da predição realizada pelo Kriging;
%       predVar - valor do erro médio cometido na predição.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula a quantidade pontos amostrais e o número de dimensões do problema
[n, k] = size(Modelo.S);

% Número de pontos que precisam ser preditos
numPredicoes = size(x, 1);

if numPredicoes == 1        % Somente uma predição é necessária
    % Calcula a distância entre o ponto e o espaço do modelo
    dx = repmat(x, n, 1) - Modelo.S;
    
    % Calcula a correlação do ponto a ser predito com os pontos amostrais
    r = feval(@correlacao, (10.^Modelo.theta(1, :)), dx);
    
    % Calcula a predição
    yHat = Modelo.mu + (Modelo.gamma * r);
    
    % Calcula o erro médio na predição
    rt = Modelo.U' \ r;
    u = Modelo.Ft'*rt - 1;
    v = Modelo.G \ u;
    
    predVar = repmat(Modelo.sigma2, 1, 1) .* repmat((1 + sum(v.^2) - sum(rt.^2))', 1, 1);
else
    % Inicializa o vetor de distâncias entre cada um dos pontos a ser predito e o espaço
    % amostral S
    dx = zeros(numPredicoes * n, k);
    
    % Inicializa um vetor auxiliar para as posições das distâncias
    kk = 1:n;
    
    % Calcula as distâncias
    for k = 1:numPredicoes
        dx(kk, :) = repmat(x(k, :), n, 1) - Modelo.S;
        kk = kk + n;
    end
    
    % Inicializa o vetor de uns
    one = ones(numPredicoes, 1);
    
    % Calcula a correlação do ponto a ser predito com os pontos amostrais
    r = feval(@correlacao, (10.^Modelo.theta(1, :)), dx);
    r = reshape(r, n, numPredicoes);
    
    % Calcula o vetor de predições
    yHat = repmat(Modelo.mu, numPredicoes, 1) + (Modelo.gamma * r)';
    
    % Calcula o vetor de erros médios nas precições
    rt = Modelo.U' \ r;
    u = Modelo.G \ (Modelo.Ft' * rt - one');
    
    predVar = repmat(Modelo.sigma2, numPredicoes, 1) .* ...
              repmat((1 + u.^2 - sum(rt.^2))', 1, 1);
end
return
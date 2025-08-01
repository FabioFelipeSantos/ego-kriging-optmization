function [yHat, predVar] = predKrg(x, Modelo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Rotina de predi��o utilizando o metamodelo em Kriging.
%
%   Entrada:
%       x - vetor de k vari�veis que ser� predita pelo metamodelo;
%       Modelo - informa��es sobre o modelo.
%
%   Sa�da:
%       yHat - valor escalar da predi��o realizada pelo Kriging;
%       predVar - valor do erro m�dio cometido na predi��o.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula a quantidade pontos amostrais e o n�mero de dimens�es do problema
[n, k] = size(Modelo.S);

% N�mero de pontos que precisam ser preditos
numPredicoes = size(x, 1);

if numPredicoes == 1        % Somente uma predi��o � necess�ria
    % Calcula a dist�ncia entre o ponto e o espa�o do modelo
    dx = repmat(x, n, 1) - Modelo.S;
    
    % Calcula a correla��o do ponto a ser predito com os pontos amostrais
    r = feval(@correlacao, (10.^Modelo.theta(1, :)), dx);
    
    % Calcula a predi��o
    yHat = Modelo.mu + (Modelo.gamma * r);
    
    % Calcula o erro m�dio na predi��o
    rt = Modelo.U' \ r;
    u = Modelo.Ft'*rt - 1;
    v = Modelo.G \ u;
    
    predVar = repmat(Modelo.sigma2, 1, 1) .* repmat((1 + sum(v.^2) - sum(rt.^2))', 1, 1);
else
    % Inicializa o vetor de dist�ncias entre cada um dos pontos a ser predito e o espa�o
    % amostral S
    dx = zeros(numPredicoes * n, k);
    
    % Inicializa um vetor auxiliar para as posi��es das dist�ncias
    kk = 1:n;
    
    % Calcula as dist�ncias
    for k = 1:numPredicoes
        dx(kk, :) = repmat(x(k, :), n, 1) - Modelo.S;
        kk = kk + n;
    end
    
    % Inicializa o vetor de uns
    one = ones(numPredicoes, 1);
    
    % Calcula a correla��o do ponto a ser predito com os pontos amostrais
    r = feval(@correlacao, (10.^Modelo.theta(1, :)), dx);
    r = reshape(r, n, numPredicoes);
    
    % Calcula o vetor de predi��es
    yHat = repmat(Modelo.mu, numPredicoes, 1) + (Modelo.gamma * r)';
    
    % Calcula o vetor de erros m�dios nas preci��es
    rt = Modelo.U' \ r;
    u = Modelo.G \ (Modelo.Ft' * rt - one');
    
    predVar = repmat(Modelo.sigma2, numPredicoes, 1) .* ...
              repmat((1 + u.^2 - sum(rt.^2))', 1, 1);
end
return
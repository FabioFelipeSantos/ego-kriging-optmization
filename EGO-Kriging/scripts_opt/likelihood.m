function [fLike, krgFit] = likelihood(theta, Modelo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcula o negativo da função de Verossimilhança (ln-concentrada)
%
%   Entrada:       
%       theta - Vetor com os parâmetros do Kriging;
%       Modelo - Informações acerca do modelo a ser substituído.
%
%   Saída:
%       fLike - Verossimilhança Ln-Concentrada * (- 1) para a minimização;
%       krgFit - parâmetros da modelagem em Kriging da função:
%                S - espaço amostral de substituição;
%                y - valores da função objetivo em S;
%			 theta - parâmetro da função de correlação entre os pontos;
%			    mu - média otimizada pela verossimilhança;
%			 gamma - produto (R^(-1))*(y - one*mu);
%			sigma2 - ariância otimizada pela verossimilhança;
%			   Psi - Matriz de correlação dos termos (forma esparsa);
%			     U - Fatoração da matriz de correlação (forma esparsa);
%			    Ft - Solução de (R^(-1))*one;
%			     G - Matriz R da fatoração QR de Ft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parte Inicial da Função
% Recolhe os pontos iniciais e os valores da função objetivo
X = Modelo.S;
y = Modelo.y;

% Calcula o número de pontos iniciais amostrais
[n, k] = size(X);

% Determina o número de valores de Theta a serem analisados
numPontos = size(theta, 1);

% Armazena os valores de theta sem a normalização logarítmica
theta = 10.^theta;

%% Calcula a distância entre os pontos da amostra
% Determina o número de distâncias existentes no espaço amostral
numDist = n*(n-1) / 2;

% Inicializa as matrizes de índices e de distâncias
ij = zeros(numDist, 2);
D = zeros(numDist, k);

% Cálculo das distâncias
auxPosicaoi = 0;        % Auxiliar no índice do primeiro elemento
for ii = 1:(n - 1)       % Percorre todos os linhas da matriz da amostra
    % Monta o vetor com o primeiro ponto da distância
    auxPosicaoi = auxPosicaoi(end) + (1:(n - ii));
    
    % Combina o primeiro ponto com os demais e organiza na matriz de índices
    ij(auxPosicaoi, :) = [repmat(ii, n - ii, 1) ((ii+1):n)'];
    
    % Calcula a diferença entre as coordenadas dos pontos
    D(auxPosicaoi, :) = repmat(X(ii,:), n-ii, 1) - X((ii+1):n, :);
end

%% Calcula o valor negativo da logaritmo da verossimilhança
if ~isfield(Modelo, 'fitKrg')
    krgFit = [];
end

% Inicialização do vetor de uns e do vetor de valores da função
one = ones(n,1);
fLike = zeros(numPontos, 1);

for i = 1:numPontos         % Percorre todos os thetas que precisam ser avaliados
    % Calcula a correlação entre as distâncias
    r = feval(@correlacao, theta(i, :), D);
    
    % Variáveis de auxílio para montagem da matriz esparsa
    idx = find(r > 0);
    o = (1:n)';
    
    % Monta a matriz esparsa de correlação dos pontos
    Psi = sparse([ij(idx,1); o], [ij(idx,2); o], [r(idx) ;one + (10 + n)*eps]);
    
    % Calcula a fatoração de Cholesky da matriz
    [U, pp] = chol(Psi);
    
    % Verifica se a matriz não é positiva definida
    if pp, fLike(i)=1e5; krgFit = []; continue, end
    
    % Calcula a média otimizada e a variância otimizada
    Ft = U' \ one;
    
    % Realiza a fatoração QR
    [Q, G] = qr(Ft, 0);
    
    % Verifica o mal condicionamento da matriz G da fatoração QR
    if rcond(G) < 1e-10, fLike(i)=1e5; krgFit = []; continue, end
    
    % Calcula a média otimizada
    Yt = U' \ y;
    mu = G \ (Q' * Yt);
    
    % Calcula a variância otimizada
    rho = Yt - Ft * mu;
    Sigma2 = sum(rho.^2) / n;
    
    if isfield(Modelo, 'fitKrg')
        % Realiza a saída dos dados o modelo Kriging para a fu
        gamma = (U \ rho)';
        krgFit = struct('S', X,...
                        'y', y,...
                        'theta', log10(theta(i,:)),...
                        'mu', mu,...
                        'gamma', gamma,...
                        'sigma2', Sigma2,...
                        'Psi', Psi,...
                        'U', U,...
                        'Ft', Ft,...
                        'G', G);
        return
    end
    
    % Soma dos lns da diagonal principal para se encontrar ln(abs(det(Psi)))
    LnDetPsi = 2 * sum(log(abs(diag(U))));
    
    % Cálculo do valor negativo da função Ln-Concentrada do Likelihood
    fLike(i) = -1 * (-(n / 2) * log(Sigma2) - 0.5 * LnDetPsi);
end
return
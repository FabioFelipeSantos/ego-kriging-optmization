function [fLike, krgFit] = likelihood(theta, Modelo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcula o negativo da fun��o de Verossimilhan�a (ln-concentrada)
%
%   Entrada:       
%       theta - Vetor com os par�metros do Kriging;
%       Modelo - Informa��es acerca do modelo a ser substitu�do.
%
%   Sa�da:
%       fLike - Verossimilhan�a Ln-Concentrada * (- 1) para a minimiza��o;
%       krgFit - par�metros da modelagem em Kriging da fun��o:
%                S - espa�o amostral de substitui��o;
%                y - valores da fun��o objetivo em S;
%			 theta - par�metro da fun��o de correla��o entre os pontos;
%			    mu - m�dia otimizada pela verossimilhan�a;
%			 gamma - produto (R^(-1))*(y - one*mu);
%			sigma2 - ari�ncia otimizada pela verossimilhan�a;
%			   Psi - Matriz de correla��o dos termos (forma esparsa);
%			     U - Fatora��o da matriz de correla��o (forma esparsa);
%			    Ft - Solu��o de (R^(-1))*one;
%			     G - Matriz R da fatora��o QR de Ft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parte Inicial da Fun��o
% Recolhe os pontos iniciais e os valores da fun��o objetivo
X = Modelo.S;
y = Modelo.y;

% Calcula o n�mero de pontos iniciais amostrais
[n, k] = size(X);

% Determina o n�mero de valores de Theta a serem analisados
numPontos = size(theta, 1);

% Armazena os valores de theta sem a normaliza��o logar�tmica
theta = 10.^theta;

%% Calcula a dist�ncia entre os pontos da amostra
% Determina o n�mero de dist�ncias existentes no espa�o amostral
numDist = n*(n-1) / 2;

% Inicializa as matrizes de �ndices e de dist�ncias
ij = zeros(numDist, 2);
D = zeros(numDist, k);

% C�lculo das dist�ncias
auxPosicaoi = 0;        % Auxiliar no �ndice do primeiro elemento
for ii = 1:(n - 1)       % Percorre todos os linhas da matriz da amostra
    % Monta o vetor com o primeiro ponto da dist�ncia
    auxPosicaoi = auxPosicaoi(end) + (1:(n - ii));
    
    % Combina o primeiro ponto com os demais e organiza na matriz de �ndices
    ij(auxPosicaoi, :) = [repmat(ii, n - ii, 1) ((ii+1):n)'];
    
    % Calcula a diferen�a entre as coordenadas dos pontos
    D(auxPosicaoi, :) = repmat(X(ii,:), n-ii, 1) - X((ii+1):n, :);
end

%% Calcula o valor negativo da logaritmo da verossimilhan�a
if ~isfield(Modelo, 'fitKrg')
    krgFit = [];
end

% Inicializa��o do vetor de uns e do vetor de valores da fun��o
one = ones(n,1);
fLike = zeros(numPontos, 1);

for i = 1:numPontos         % Percorre todos os thetas que precisam ser avaliados
    % Calcula a correla��o entre as dist�ncias
    r = feval(@correlacao, theta(i, :), D);
    
    % Vari�veis de aux�lio para montagem da matriz esparsa
    idx = find(r > 0);
    o = (1:n)';
    
    % Monta a matriz esparsa de correla��o dos pontos
    Psi = sparse([ij(idx,1); o], [ij(idx,2); o], [r(idx) ;one + (10 + n)*eps]);
    
    % Calcula a fatora��o de Cholesky da matriz
    [U, pp] = chol(Psi);
    
    % Verifica se a matriz n�o � positiva definida
    if pp, fLike(i)=1e5; krgFit = []; continue, end
    
    % Calcula a m�dia otimizada e a vari�ncia otimizada
    Ft = U' \ one;
    
    % Realiza a fatora��o QR
    [Q, G] = qr(Ft, 0);
    
    % Verifica o mal condicionamento da matriz G da fatora��o QR
    if rcond(G) < 1e-10, fLike(i)=1e5; krgFit = []; continue, end
    
    % Calcula a m�dia otimizada
    Yt = U' \ y;
    mu = G \ (Q' * Yt);
    
    % Calcula a vari�ncia otimizada
    rho = Yt - Ft * mu;
    Sigma2 = sum(rho.^2) / n;
    
    if isfield(Modelo, 'fitKrg')
        % Realiza a sa�da dos dados o modelo Kriging para a fu
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
    
    % C�lculo do valor negativo da fun��o Ln-Concentrada do Likelihood
    fLike(i) = -1 * (-(n / 2) * log(Sigma2) - 0.5 * LnDetPsi);
end
return
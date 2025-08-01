function [Resultado, fitKrg] = otimEGO(Problema)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Programa principal do EGO utilizando o Kriging. Essa rotina � a respons�vel por
%   organizar e executar todas as etapas necess�rias para a otimiza��o de problemas com 
%   restri��o utilizando o Kriging como modelo substituto para as restri��es.
%
%   Entrada:
%       Problema - struct com todas as inform��es sobre o problema;
%       dadosPSO - struct com os par�metros do PSO (argumento opcional).
%
%   Sa�da:
%       Modelo - todas as informa��es sobre o problema otimizado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recolhe as vari�veis de entrada do problema
% PAR�METROS OBRIGAT�RIOS
% Fun��o avaliada
if isfield(Problema, 'func')
    func = str2func(Problema.func);
else
    error('� preciso definir o nome da fun��o objetivo a ser otimizada!');
end
% N�mero m�ximo de avalia��es permitido
if isfield(Problema, 'maxFE')
    maxFE = Problema.maxFE;
else
    error('� preciso definir o n�mero m�ximo de avalia��es da fun��o objetivo!');
end
% N�mero de dimens�es
if isfield(Problema, 'k')
    k = Problema.k;
else
    error('� preciso definir o n�mero de dimens�es do problema!');
end
% Lower bound
if isfield(Problema, 'LB')
    LB = Problema.LB;
else
    error('� preciso definir o Lower Bound do problema (limite inferior)!');
end
% Upper bound
if isfield(Problema, 'UB')
    UB = Problema.UB;
else
    error('� preciso definir o Upper Bound do problema (limite superior)!');
end
% Verifica o par�metro c�lculo vetorizado
if Problema.calculoVetorizado ~= 0 && Problema.calculoVetorizado ~= 1
    error('O valor do par�metro calculoVetorizado deve ser 0 ou 1!');
end
if isfield(Problema, 'tolMinEI')
    tolMinEI = Problema.tolMinEI;
else
    tolMinEI = 'NaoTemToleranciaMinima';
end
% N�mero de pontos amostrais iniciais
if isfield(Problema, 'n')
    n = Problema.n;
else
    if k <= 5
        n = ceil(6.5 * k);
    else
        n = ceil (4.5 * k);
    end
end
% Define os limites inferiores e superiores de pesquisa do par�metro theta
lbT = -4*ones(1, k);
ubT =  4*ones(1, k);

if isfield(Problema, 'graf')
    if Problema.graf ~= 0 && Problema.graf ~= 1
        error('O par�metro graf s� pode assumir os valores 0 ou 1.');
    end
else
    Problema.graf = 0;
end

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%                                       PARTE 1
%              C�LCULO DA AMOSTRA INICIAL E DOS METAMODELOS DAS RESTRI��ES
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%
% Montagem do Plano da Amostra inicial
fprintf('Efficient Global Optimization...\n')
if isfield(Problema, 'S')
    S = Problema.S;
else
    S = lhsdesign(n, k,'iterations',20);
end

% Avalia��o das fun��es no plano de amostra inicial
% Calcula os valores da fun��o a ser substitu�da sobre todos os pontos de interpola��o
% para o Kriging Determin�stico
% Inicializa o vetor de valores da fun��o objetivo e das restri��es
f = zeros(n, 1);

% Cria��o de uma barra de progresso para o n�mero de avalia��es das restri��es
barra_progresso = waitbar(0, sprintf('Avalia��es da Fun��o = %d', 0),...
                                      'Name', 'C�lculo de avalia��es da fun��o objetivo');

% Inicializa a contagem de avalia��es da fun��o objetivo
numFE = 0;           % N�mero de Function Evaluations (FE)

if Problema.calculoVetorizado
    % Converte o espa�o de busca
    d = LB + (UB - LB).*S;
    
    % Calcula a fun��o objetivo
    f = feval(func, d);
else
    for i = 1:n         % Percorre todos os pontos amostrais iniciais
        % Converte o espa�o de busca
        d = LB + (UB - LB).*S(i, :);
        
        % Calcula a fun��o objetivo no ponto dado
        f(i) = func(d);
    end
end

% Cria o vetor que armazenar� as melhores solu��es e os pontos onde foram obtidos
% Encontra o melhor m�nimo at� o momento
[yMin, ind] = min(f);
Ymin = [ind, yMin];
Dmin = LB + (UB - LB).*S(ind, :);
fprintf('\tValor Otimo(%d) = %.10f\n', Ymin(end,1), Ymin(end,2));

% Atualiza o n�mero de fun��es calculadas
numFE = numFE + n;

waitbar(numFE / maxFE, barra_progresso, sprintf('Avalia��es da Fun��o = %d', numFE));
%% C�lculo da superf�cie de resposta utilizando o Kriging
% Busca dos par�metros do metamodelo para as restri��es
if ~isfield(Problema, 'nPop')
    optGA = optimoptions('ga','display','off','UseVectorized',true,'PopulationSize',20*k);
else
    optGA = optimoptions('ga','display','off','UseVectorized',true,...
            'PopulationSize',Problema.nPop);
end

Modelo.S = S;
Modelo.y = f;
Theta = ga(@(x)likelihood(x,Modelo),k,[],[],[],[],lbT,ubT,[],optGA);
Modelo.fitKrg = [];    
[~, fitKrg] = likelihood(Theta, Modelo);
clear Modelo;
 
%% Plot da parte inicial dos resultados
if k == 1 && Problema.graf == 1
    nnn = 200;
    x = (LB:(UB - LB)/(nnn-1):UB)';
    pontos = (x - LB) ./ (UB - LB);
    fobj = feval(func, x);
    [fres, ~] = predKrg(pontos, fitKrg);
    
    figure('units','normalized','outerposition',[0 0 1 1],'Name','Fun��o');
    hold on
    box on
    grid on
    ax = gca();
    ax.XLim = [LB UB];
    ax.XTick = LB:((UB - LB) / 20):UB;
    xlabel('x','FontSize',14);
    ylabel('f(x)','FontSize',14);
    plot(x, fobj,'--b','LineWidth', 1.8);
    plot(x, fres,'-r','LineWidth', 2);
    dd = LB + (UB-LB).*S;
    plot(dd,f,'o','MarkerSize',9.5,'MarkerFaceColor',[4/5, 0, 4/5],'MarkerEdgeColor',[77/255, 0, 77/255]);
    drawnow;
end
if k == 2 && Problema.graf == 1
    nnn = 80;
    x = LB(1):(UB(1) - LB(1))/(nnn-1):UB(1);
    y = LB(2):(UB(2) - LB(2))/(nnn-1):UB(2);
    xx = 0:1/(nnn-1):1;
    yy = xx;
    [Y, X] = meshgrid(yy, xx);
    Y = reshape(Y, [length(xx)^2, 1]);
    X = reshape(X, [length(xx)^2, 1]);
    pontos = [X Y];
    [fres, ~] = predKrg(pontos, fitKrg);
    fres = reshape(fres, nnn, nnn)';
    
    figure('units','normalized','outerposition',[0 0.1 0.5 0.85],'Name','Curvas de N�vel');
    hax1 = axes;
    hold on
    box on
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    title('Curvas de N�vel da Fun��o Objetivo','FontSize',20);
    contour(x, y, fres, 40,'LineWidth', 1.8)
    dd = LB + (UB-LB).*S;
    plot(dd(:,1),dd(:,2),'o','MarkerSize',9.5,'MarkerFaceColor',[4/5, 0, 4/5],'MarkerEdgeColor',[77/255, 0, 77/255]);
    drawnow;
    
    figure('units','normalized','outerposition',[0.5 0.1 0.5 0.85],'Name','Superf�cie');
    hax2 = axes;
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('f(x,y)','FontSize',14)
    title('Superf�cie da fun��o objetivo')
    surf(x,y,fres,'FaceAlpha',0.99)
end
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%                                       PARTE 2
%            ADI��O DE INFILL POINTS (IP) E EFFICIENT GLOBAL OPTMIZATION (EGO)
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%% Rotina de Otimiza��o
numPontosIniciais = n;

% Inicializa os controladores do processo de adi��o de IP
it = 0;                 % N�mero de Itera��es Realizadas

% Inicializa os m�ltiplicadores de Lagrange e o coeficiente de penaliza��o
criterioParada = [];

while numFE < maxFE          % Realiza as itera��es de otimiza��o    
    % Atualiza o n�mero de itera��es
    it = it + 1;
    %fprintf('\tItera��o %d:', it)
    
    % Procura pelo m�ximo da m�trica de adi��o do IP
    %fprintf('\n\t\tCalculando o m�ximo do Expected Improvement: ')    
    [dNovo, EI] = ga(@(x)expImp(x,fitKrg),k,[],[],[],[],zeros(1, k),ones(1, k),[],optGA);
    %fprintf('Resolvido!\n\t\tNovo Ponto Calculado.')
    
    d = LB + (UB - LB).*dNovo;
    
    % Calcula a fun��o objetivo no novo ponto
    fnew = func(d); 
    
    if k == 1 && Problema.graf == 1
        p1 = plot(d,fnew,'s','MarkerSize',12,'MarkerFaceColor',[0.3, 0.65, 1],'MarkerEdgeColor',[2/5, 0, 2/5]);
        drawnow;
    end
    
    if k == 2 && Problema.graf == 1
        p1 = plot(hax1,d(1),d(2),'s','MarkerSize',12,'MarkerFaceColor',[0.3, 0.65, 1],'MarkerEdgeColor',[2/5, 0, 2/5]);
        drawnow;
    end
        
    % Acumula o n�mero de avalia��es da fun��o
    numFE = numFE + 1;
    
    waitbar(numFE / maxFE, barra_progresso, sprintf('Avalia��es da Fun��o = %d', numFE));

    % Adiciona esses valores na amostra inicial
    S = [S; dNovo]; %#ok<AGROW>
    f = [f; fnew];  %#ok<AGROW>
    
    % Verifica se o valor obtido nessa itera��o � melhor que o anterior
    if fnew < Ymin(end, 2)
        % Atualiza a posi��o do m�nimo
        pos = numPontosIniciais + it;
        % Adiciona o m�nimo obtido
        Ymin = [Ymin; [pos, fnew]];
        % Adiciona o novo ponto minimizador
        Dmin = [Dmin; d];
        fprintf('\tValor Otimo(%d) = %.10f\n', Ymin(end,1), Ymin(end,2));
        
        if k == 1 && Problema.graf == 1
            delete(p1)
            p1 = plot(d,fnew,'p','MarkerSize',18,'MarkerFaceColor',[1,4/5,0],'MarkerEdgeColor',[2/5, 0, 2/5]);
            drawnow;
        end
        
        if k == 2 && Problema.graf == 1
            delete(p1)
            p1 = plot(hax1,d(1),d(2),'p','MarkerSize',18,'MarkerFaceColor',[1,4/5,0],'MarkerEdgeColor',[2/5, 0, 2/5]);
            drawnow;
        end
    end
    
    %fprintf('\n\t\tReavaliando Modelo: ')
    % Atualiza n e as matrizes de correla��o
    n = n + 1;
    
    % Cria o Modelo
    Modelo.S = S;
    Modelo.y = f;
    Theta = ga(@(x)likelihood(x,Modelo),k,[],[],[],[],lbT,ubT,[],optGA);
    Modelo.fitKrg = [];
    [~, fitKrg] = likelihood(Theta, Modelo);
    clear Modelo;
    
    if k == 1 && Problema.graf == 1
        fres = predKrg(pontos, fitKrg);
        h = findobj('LineWidth', 2);
        delete(h);
        plot(x, fres,'-r','LineWidth', 2);
        drawnow;
    end
    
    if k == 2 && Problema.graf == 1
        [fres, ~] = predKrg(pontos, fitKrg);
        fres = reshape(fres, nnn, nnn)';
        h = findobj('LineWidth', 1.8);
        delete(h);
        contour(hax1,x, y, fres, 40,'LineWidth', 1.8)
        drawnow;
        
        h2 = findobj('FaceAlpha',0.99);
        delete(h2);
        surf(hax2,x,y,fres,'FaceAlpha',0.99);
    end
    
    if ~strcmp(tolMinEI, 'NaoTemToleranciaMinima')
        if abs(EI) < tolMinEI
            criterioParada = 'Toler�ncia do EI alcan�ada';
            break;
        end
    end
end
delete(barra_progresso)
if isempty(criterioParada)
    criterioParada = 'N�mero M�ximo de Avalia��es das Restri��es';
end

fprintf('T�rmino das itera��es!\n');

%% Aramazena a struct ModelInfo na vari�vel Modelo
Resultado.MelhorMinimizador = Dmin(end, :);
Resultado.MelhorValorOtimo = Ymin(end);
Resultado.PosMelhorMin = ind;
Resultado.numFE = numFE;
Resultado.numIPAdd = n - numPontosIniciais;
Resultado.numPontosIniciais = numPontosIniciais;
Resultado.fobj = f;
Resultado.Histfmin = Ymin;
Resultado.Histdmin = Dmin;
Resultado.CriterioParada = criterioParada;
return
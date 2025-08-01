function [Resultado, fitKrg] = otimEGO(Problema)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Programa principal do EGO utilizando o Kriging. Essa rotina é a responsável por
%   organizar e executar todas as etapas necessárias para a otimização de problemas com 
%   restrição utilizando o Kriging como modelo substituto para as restrições.
%
%   Entrada:
%       Problema - struct com todas as informções sobre o problema;
%       dadosPSO - struct com os parâmetros do PSO (argumento opcional).
%
%   Saída:
%       Modelo - todas as informações sobre o problema otimizado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recolhe as variáveis de entrada do problema
% PARÂMETROS OBRIGATÓRIOS
% Função avaliada
if isfield(Problema, 'func')
    func = str2func(Problema.func);
else
    error('É preciso definir o nome da função objetivo a ser otimizada!');
end
% Número máximo de avaliações permitido
if isfield(Problema, 'maxFE')
    maxFE = Problema.maxFE;
else
    error('É preciso definir o número máximo de avaliações da função objetivo!');
end
% Número de dimensões
if isfield(Problema, 'k')
    k = Problema.k;
else
    error('É preciso definir o número de dimensões do problema!');
end
% Lower bound
if isfield(Problema, 'LB')
    LB = Problema.LB;
else
    error('É preciso definir o Lower Bound do problema (limite inferior)!');
end
% Upper bound
if isfield(Problema, 'UB')
    UB = Problema.UB;
else
    error('É preciso definir o Upper Bound do problema (limite superior)!');
end
% Verifica o parâmetro cálculo vetorizado
if Problema.calculoVetorizado ~= 0 && Problema.calculoVetorizado ~= 1
    error('O valor do parâmetro calculoVetorizado deve ser 0 ou 1!');
end
if isfield(Problema, 'tolMinEI')
    tolMinEI = Problema.tolMinEI;
else
    tolMinEI = 'NaoTemToleranciaMinima';
end
% Número de pontos amostrais iniciais
if isfield(Problema, 'n')
    n = Problema.n;
else
    if k <= 5
        n = ceil(6.5 * k);
    else
        n = ceil (4.5 * k);
    end
end
% Define os limites inferiores e superiores de pesquisa do parâmetro theta
lbT = -4*ones(1, k);
ubT =  4*ones(1, k);

if isfield(Problema, 'graf')
    if Problema.graf ~= 0 && Problema.graf ~= 1
        error('O parâmetro graf só pode assumir os valores 0 ou 1.');
    end
else
    Problema.graf = 0;
end

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%                                       PARTE 1
%              CÁLCULO DA AMOSTRA INICIAL E DOS METAMODELOS DAS RESTRIÇÕES
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%
% Montagem do Plano da Amostra inicial
fprintf('Efficient Global Optimization...\n')
if isfield(Problema, 'S')
    S = Problema.S;
else
    S = lhsdesign(n, k,'iterations',20);
end

% Avaliação das funções no plano de amostra inicial
% Calcula os valores da função a ser substituída sobre todos os pontos de interpolação
% para o Kriging Determinístico
% Inicializa o vetor de valores da função objetivo e das restrições
f = zeros(n, 1);

% Criação de uma barra de progresso para o número de avaliações das restrições
barra_progresso = waitbar(0, sprintf('Avaliações da Função = %d', 0),...
                                      'Name', 'Cálculo de avaliações da função objetivo');

% Inicializa a contagem de avaliações da função objetivo
numFE = 0;           % Número de Function Evaluations (FE)

if Problema.calculoVetorizado
    % Converte o espaço de busca
    d = LB + (UB - LB).*S;
    
    % Calcula a função objetivo
    f = feval(func, d);
else
    for i = 1:n         % Percorre todos os pontos amostrais iniciais
        % Converte o espaço de busca
        d = LB + (UB - LB).*S(i, :);
        
        % Calcula a função objetivo no ponto dado
        f(i) = func(d);
    end
end

% Cria o vetor que armazenará as melhores soluções e os pontos onde foram obtidos
% Encontra o melhor mínimo até o momento
[yMin, ind] = min(f);
Ymin = [ind, yMin];
Dmin = LB + (UB - LB).*S(ind, :);
fprintf('\tValor Otimo(%d) = %.10f\n', Ymin(end,1), Ymin(end,2));

% Atualiza o número de funções calculadas
numFE = numFE + n;

waitbar(numFE / maxFE, barra_progresso, sprintf('Avaliações da Função = %d', numFE));
%% Cálculo da superfície de resposta utilizando o Kriging
% Busca dos parâmetros do metamodelo para as restrições
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
    
    figure('units','normalized','outerposition',[0 0 1 1],'Name','Função');
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
    
    figure('units','normalized','outerposition',[0 0.1 0.5 0.85],'Name','Curvas de Nível');
    hax1 = axes;
    hold on
    box on
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    title('Curvas de Nível da Função Objetivo','FontSize',20);
    contour(x, y, fres, 40,'LineWidth', 1.8)
    dd = LB + (UB-LB).*S;
    plot(dd(:,1),dd(:,2),'o','MarkerSize',9.5,'MarkerFaceColor',[4/5, 0, 4/5],'MarkerEdgeColor',[77/255, 0, 77/255]);
    drawnow;
    
    figure('units','normalized','outerposition',[0.5 0.1 0.5 0.85],'Name','Superfície');
    hax2 = axes;
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('f(x,y)','FontSize',14)
    title('Superfície da função objetivo')
    surf(x,y,fres,'FaceAlpha',0.99)
end
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%                                       PARTE 2
%            ADIÇÃO DE INFILL POINTS (IP) E EFFICIENT GLOBAL OPTMIZATION (EGO)
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%% Rotina de Otimização
numPontosIniciais = n;

% Inicializa os controladores do processo de adição de IP
it = 0;                 % Número de Iterações Realizadas

% Inicializa os múltiplicadores de Lagrange e o coeficiente de penalização
criterioParada = [];

while numFE < maxFE          % Realiza as iterações de otimização    
    % Atualiza o número de iterações
    it = it + 1;
    %fprintf('\tIteração %d:', it)
    
    % Procura pelo máximo da métrica de adição do IP
    %fprintf('\n\t\tCalculando o máximo do Expected Improvement: ')    
    [dNovo, EI] = ga(@(x)expImp(x,fitKrg),k,[],[],[],[],zeros(1, k),ones(1, k),[],optGA);
    %fprintf('Resolvido!\n\t\tNovo Ponto Calculado.')
    
    d = LB + (UB - LB).*dNovo;
    
    % Calcula a função objetivo no novo ponto
    fnew = func(d); 
    
    if k == 1 && Problema.graf == 1
        p1 = plot(d,fnew,'s','MarkerSize',12,'MarkerFaceColor',[0.3, 0.65, 1],'MarkerEdgeColor',[2/5, 0, 2/5]);
        drawnow;
    end
    
    if k == 2 && Problema.graf == 1
        p1 = plot(hax1,d(1),d(2),'s','MarkerSize',12,'MarkerFaceColor',[0.3, 0.65, 1],'MarkerEdgeColor',[2/5, 0, 2/5]);
        drawnow;
    end
        
    % Acumula o número de avaliações da função
    numFE = numFE + 1;
    
    waitbar(numFE / maxFE, barra_progresso, sprintf('Avaliações da Função = %d', numFE));

    % Adiciona esses valores na amostra inicial
    S = [S; dNovo]; %#ok<AGROW>
    f = [f; fnew];  %#ok<AGROW>
    
    % Verifica se o valor obtido nessa iteração é melhor que o anterior
    if fnew < Ymin(end, 2)
        % Atualiza a posição do mínimo
        pos = numPontosIniciais + it;
        % Adiciona o mínimo obtido
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
    % Atualiza n e as matrizes de correlação
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
            criterioParada = 'Tolerância do EI alcançada';
            break;
        end
    end
end
delete(barra_progresso)
if isempty(criterioParada)
    criterioParada = 'Número Máximo de Avaliações das Restrições';
end

fprintf('Término das iterações!\n');

%% Aramazena a struct ModelInfo na variável Modelo
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
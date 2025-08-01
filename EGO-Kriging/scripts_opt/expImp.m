function EI = expImp(d, fitKrg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Rotina para o cálculo do valor do Expected Improvement utilizando a Integração de 
%   Monte Carlo.
%
%   Entrada:
%       d - indíviduo sob avaliação;
%       Modelo - informações sobre o modelo a ser explorado.
%
%   Saída:
%       valor - o valor do esperança do Lagrangeano Aumentado.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula o valor mínimo obtido até o momento
ymin = min(fitKrg.y);

% Calcula o valor da predição bem como o valor de Ssqr (a variância do erro)
[f, s] = predKrg(d, fitKrg);

% Verifica se há algum erro menor que zero (devido a erro numérico)
aux = find(s(:) < 0, 1);
if ~isempty(aux)
    % Se tiver, transforma o erro em 0 (valor negativo não é permitido)
    s(aux) = 0;
end

% Calcula o RMSE (Root Mean Square Error)
s = sqrt(s);

% Cálculo do Expected Improvement
termo1 = (ymin - f) .* (0.5 + 0.5*erf((ymin - f) ./ (s * sqrt(2))));
termo2 = (1 / sqrt(2*pi)) * s .* exp(-(((ymin - f) .^ 2) ./ (2 * (s .^ 2))));
EI = -(termo1 + termo2);
return
function EI = expImp(d, fitKrg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Rotina para o c�lculo do valor do Expected Improvement utilizando a Integra��o de 
%   Monte Carlo.
%
%   Entrada:
%       d - ind�viduo sob avalia��o;
%       Modelo - informa��es sobre o modelo a ser explorado.
%
%   Sa�da:
%       valor - o valor do esperan�a do Lagrangeano Aumentado.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula o valor m�nimo obtido at� o momento
ymin = min(fitKrg.y);

% Calcula o valor da predi��o bem como o valor de Ssqr (a vari�ncia do erro)
[f, s] = predKrg(d, fitKrg);

% Verifica se h� algum erro menor que zero (devido a erro num�rico)
aux = find(s(:) < 0, 1);
if ~isempty(aux)
    % Se tiver, transforma o erro em 0 (valor negativo n�o � permitido)
    s(aux) = 0;
end

% Calcula o RMSE (Root Mean Square Error)
s = sqrt(s);

% C�lculo do Expected Improvement
termo1 = (ymin - f) .* (0.5 + 0.5*erf((ymin - f) ./ (s * sqrt(2))));
termo2 = (1 / sqrt(2*pi)) * s .* exp(-(((ymin - f) .^ 2) ./ (2 * (s .^ 2))));
EI = -(termo1 + termo2);
return
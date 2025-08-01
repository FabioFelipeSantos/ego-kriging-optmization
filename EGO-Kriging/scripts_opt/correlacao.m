function r = correlacao(theta, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Correla��o Gaussiana para o ajuste do Kriging
%
%   Entrada:
%       theta - par�metros da fun��o de correla��o;
%       d - matriz (n x k) com as diferen�as entre os pontos amostrais.
%
%   Sa�da:
%       r - correla��o desses pontos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula o n�mero de diferen�as e a dimens�o do problema
n = size(d, 1);

% Ajusta um vetor com o par�metro theta para cada dist�ncia
auxTheta = repmat(-theta(:)', n, 1);

% Calcula o produto da dist�ncia pelo par�metro theta
auxFinal = (d.^2) .* auxTheta;

% Calcula a correla��o
r = exp(sum(auxFinal, 2));
return
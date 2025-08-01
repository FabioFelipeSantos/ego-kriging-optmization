%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Otimiza��o de Fun��es Determin�sticas Utilizando o Kriging com a abordagem Efficient 
%   Global Optimization (EGO).
%   Autor: F�bio Felipe dos Santos Nascentes
%   Orientador: Rafael Holdorf Lopez
%   Universidade Federal de Santa Catarina
%   Data �ltima Atualiza��o: 31/10/2018
%
%   Copyright (c) 2018 F�bio Nascentes
%   Todos os direitos autorais reservados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	Parte Inicial (N�O NECESSITA DE ALTERA��O)
% Adiciona as pastas com arquivos .m
addpath(strcat(pwd, '\fobj'))
addpath(strcat(pwd, '\scripts_opt'))

clear;              % Limpa a �rea de trabalho das vari�veis
clc;                % Limpa a janela de resultados
close all;          % Fecha todas as janelas abertas

%%  DEFINI��O DOS PAR�METROS DE CONTROLE DO PROGRAMA E DO MODELO
% Define qual ser� a fun��o objetivo
Problema.func = 'costmtmd';

% Define o n�mero m�ximo de vezes em que as restri��es podem ser avaliadas
Problema.maxFE = 60;

% N�mero de Vari�veis
Problema.k = 4;

% Define o Lower Bound das vari�veis
Problema.LB = [10 10 10 10];

% Define o Upper Bound das vari�veis
Problema.UB =  [133000 133000 10000/2 10000/2]/2;

% Define se a fun��o objetivo poder� ser calculada de forma vetorizada
Problema.calculoVetorizado = 0;

% Par�metro opcional graf. Controla se haver� sa�da gr�fica
Problema.graf = 0;

%% Chamada da rotina de otimiza��o
[Resultado, fitKrg] = otimEGO(Problema);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Otimização de Funções Determinísticas Utilizando o Kriging com a abordagem Efficient 
%   Global Optimization (EGO).
%   Autor: Fábio Felipe dos Santos Nascentes
%   Orientador: Rafael Holdorf Lopez
%   Universidade Federal de Santa Catarina
%   Data Última Atualização: 31/10/2018
%
%   Copyright (c) 2018 Fábio Nascentes
%   Todos os direitos autorais reservados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	Parte Inicial (NÃO NECESSITA DE ALTERAÇÃO)
% Adiciona as pastas com arquivos .m
addpath(strcat(pwd, '\fobj'))
addpath(strcat(pwd, '\scripts_opt'))

clear;              % Limpa a área de trabalho das variáveis
clc;                % Limpa a janela de resultados
close all;          % Fecha todas as janelas abertas

%%  DEFINIÇÃO DOS PARÂMETROS DE CONTROLE DO PROGRAMA E DO MODELO
% Define qual será a função objetivo
Problema.func = 'costmtmd';

% Define o número máximo de vezes em que as restrições podem ser avaliadas
Problema.maxFE = 60;

% Número de Variáveis
Problema.k = 4;

% Define o Lower Bound das variáveis
Problema.LB = [10 10 10 10];

% Define o Upper Bound das variáveis
Problema.UB =  [133000 133000 10000/2 10000/2]/2;

% Define se a função objetivo poderá ser calculada de forma vetorizada
Problema.calculoVetorizado = 0;

% Parâmetro opcional graf. Controla se haverá saída gráfica
Problema.graf = 0;

%% Chamada da rotina de otimização
[Resultado, fitKrg] = otimEGO(Problema);
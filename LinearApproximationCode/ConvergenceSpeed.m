function [ meanSpeed varSpeed ] = ConvergenceSpeed( Para, arg, argvalue )
%CONVERGENCESPEED Summary of this function goes here
%   Detailed explanation goes here
    eval(['Para.' arg '= argvalue;']);
    if strcmp(arg,'alpha_1')
        Para.alpha_2 = 1-Para.alpha_1;
    end
    [ A, XSS, B, BS ] = LinearApproximation( Para);
    Bbar = Para.P(1,1)*B{1}+Para.P(1,2)*B{2};
    meanSpeed = 1- max(eig(Bbar));
    varSpeed = 1- max(eig(BS));
end


function [ R ] = fR( p,pstring,Para )
%FR Summary of this function goes here
%   Detailed explanation goes here
    eval(['Para.' pstring ' = p;']);
    [~,R] = findSteadyState(0,3,Para);
    %R = R^(1/Para.sigma);
end


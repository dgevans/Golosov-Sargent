function xstar=getCrossingPoints(coeff,V,Para)
get_x_crossing_nag= @(x) resCrossingPoints(x,coeff,V,Para)  ;
 [xstar,~,exitflag]=fsolve(@(x) get_x_crossing_nag(x) ,Para.xMin/2);
end

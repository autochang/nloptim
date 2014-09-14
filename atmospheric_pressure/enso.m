function [f, g, B] = enso(b)
	[ensoFit,ensoRes,ensoData,ensoJac]=ensoFitResJac(b);
	f = dot(ensoRes, ensoRes)/2;
	g = ensoJac' * ensoRes;
	B = ensoJac' * ensoJac;


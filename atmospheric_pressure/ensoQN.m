function [f, g, Resi, Jaco] = ensoQN(b)
	[ensoFit,ensoRes,ensoData,ensoJac]=ensoFitResJac(b);
	f = dot(ensoRes, ensoRes)/2;
	g = ensoJac' * ensoRes;
	Resi = ensoRes;
	Jaco = ensoJac;
	


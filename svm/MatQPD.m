function H = MatQPD(Data)
	YY = Data(:,1)*Data(:,1)';
	XX = Data(:,2:end)*Data(:,2:end)';
	H = YY .* XX;
	

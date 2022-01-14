%path_string=''; %here add your relevant path with throt and freeprecess functions
addpath(genpath(path_string)) 

%define saturation train
FA=155;
dT=0.005; %in ms
dur=0.5;
sub=ones(1,length(0:dT:dur)-1);
b1=[-sub,sub,sub,-sub,-sub,sub,sub,-sub,-sub,sub,sub,-sub]; 
rf=b1*FA*pi*dT/(180);

%define parameters
T1 = 2000;       % ms.
freq = [-1000:20:1000];	% Hz. frequency offset space
T2_s=logspace(-4,6,100); % in ms. T2 space


sig=zeros(length(T2_s),length(freq));

for T2_l=1:length(T2_s)

	T2=T2_s(T2_l);
	for f=1:length(freq)
	    M = [0;0;1];
	    M = throt(abs(rf(1)),angle(rf(1))) * M;	% RF Rotation.
	    for k = 2:length(rf)
			[A,B] = freeprecess(dT,T1,T2,freq(f));
			M = A*M+B;				% Propagate to next pulse.
		    	M = throt(abs(rf(k)),angle(rf(k))) * M;	% RF Rotation.
	    end;
	    sig(T2_l,f) = M(3);  
	end;
end;


rowNames = cellstr(string(T2_s));
colNames =cellstr(string(freq));
sTable = array2table(sig,'RowNames',rowNames,'VariableNames',colNames);

surf(sig)
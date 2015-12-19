%% Read RNASeq data to matrix X
[num,fileContents,raw]=xlsread('ProteomicsFiles.xlsx');

X=[];
tcgaID=[];
for i = 2:size(fileContents,1)
    if (size(fileContents{i,2},2)>0)
        i
        tcgaID=[tcgaID; fileContents{i,1}(1:12)];
        proteomicsFile = tdfread(strcat('RPPA/',fileContents{i,2}));
        fNames = fieldnames(proteomicsFile);
    	x_new = proteomicsFile.(fNames{2});
    	x_new = x_new(2:size(x_new,1),:);
    	x_new = str2num(x_new);
    	X=[X x_new];
    end
end
XT=transpose(X);

csvwrite('proteomics.csv',XT);
dlmwrite('proteomicsTcgaID.csv',tcgaID,'delimiter','');

%% Read RNASeq data to matrix X
[num,fileContents,raw]=xlsread('RNAseqV2Files.xlsx');

X=[];
tcgaID=[];
for i = 2:size(fileContents,1)
    if (size(fileContents{i,2},2)>0)
        i
        tcgaID=[tcgaID; fileContents{i,1}(1:12)];
        RNASeqFile = tdfread(strcat('RNASeqV2/',fileContents{i,2}));
        new_item = RNASeqFile.normalized_count;
        X=[X new_item];
    end
end
XT=transpose(X);

csvwrite('RNAseq.csv',XT);
dlmwrite('RNAseqTcgaID.csv',tcgaID,'delimiter','');

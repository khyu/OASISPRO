%% Read Clinical data to matrix Y
[num,txt,raw]=xlsread('ClinicalSorted.xlsx');
Y=raw;
n=size(Y,1)-1;

tcgaID=cell(n,1);

age=zeros(n,1);
gender=zeros(n,1);
race=zeros(n,1);

pathT=zeros(n,1);
pathN=zeros(n,1);
pathM=zeros(n,1);
stages=zeros(n,1);

survival=zeros(n,1);
daysLastFollowUp=zeros(n,1);
daysToDeath=zeros(n,1);

tumorNormal=zeros(n,1);

Youtput=Y;
%% Assign numeric values for each classification
for i = 2:(n+1)
    % TCGA id
    tcgaID{i-1}=Y{i,1};
    % Gender
    switch Y{i,7}
        case 'MALE'
            gender(i-1) = 1;
            Youtput{i,7} = 1;
        case 'FEMALE'
            gender(i-1) = 0;
            Youtput{i,7} = 0;
        otherwise
            gender(i-1) = -1;
            Youtput{i,7} = -1;
    end
    % Race
    switch Y{i,8}
        case 'AMERICAN INDIAN OR ALASKA NATIVE'
            race(i-1) = 1;
            Youtput{i,8} = 1;
        case 'ASIAN'
            race(i-1) = 2;
            Youtput{i,8} = 2;
        case 'BLACK OR AFRICAN AMERICAN'
            race(i-1) = 3;
            Youtput{i,8} = 3;
        case 'WHITE'
            race(i-1) = 4;
            Youtput{i,8} = 4;
        otherwise
            race(i-1) = -1;
            Youtput{i,8} = -1;
    end
    % Pathologic T
    switch Y{i,37}
        case 'T1'
            pathT(i-1) = 10;
            Youtput{i,37} = 10;
        case 'T1a'
            pathT(i-1) = 11;
            Youtput{i,37} = 11;
        case 'T1b'
            pathT(i-1) = 12;
            Youtput{i,37} = 12;
        case 'T2'
            pathT(i-1) = 20;
            Youtput{i,37} = 20;
        case 'T3'
            pathT(i-1) = 30;
            Youtput{i,37} = 30;
        case 'T4'
            pathT(i-1) = 40;
            Youtput{i,37} = 40;
        case 'T4a'
            pathT(i-1) = 41;
            Youtput{i,37} = 41;
        otherwise
            pathT(i-1) = -1;
            Youtput{i,37} = -1;
    end
    % Pathologic N
    switch Y{i,38}
        case 'N0'
            pathN(i-1) = 0;
            Youtput{i,38} = 0;
        case 'N1'
            pathN(i-1) = 10;
            Youtput{i,38} = 10;
        case 'N1a'
            pathN(i-1) = 11;
            Youtput{i,38} = 11;
        case 'N1b'
            pathN(i-1) = 12;
            Youtput{i,38} = 12;
        case 'N2'
            pathN(i-1) = 20;
            Youtput{i,38} = 20;
        case 'N3'
            pathN(i-1) = 30;
            Youtput{i,38} = 30;
        otherwise
            pathN(i-1) = -1;
            Youtput{i,38} = -1;
    end
    
    % Pathologic M
    switch Y{i,39}
        case 'M0'
            pathM(i-1) = 0;
            Youtput{i,39} = 0;
        case 'M1'
            pathM(i-1) = 10;
            Youtput{i,39} = 10;
        otherwise
            pathM(i-1) = -1;
            Youtput{i,39} = -1;
    end
    % Pathologic stage
    switch Y{i,40}
        case 'Stage I'
            stages(i-1) = 10;
            Youtput{i,40} = 10;
        case 'Stage II'
            stages(i-1) = 20;
            Youtput{i,40} = 20;
        case 'Stage III'
            stages(i-1) = 30;
            Youtput{i,40} = 30;
        case 'Stage IVA'
            stages(i-1) = 41;
            Youtput{i,40} = 41;
        case 'Stage IVC'
            stages(i-1) = 43;
            Youtput{i,40} = 43;
        otherwise
            stages(i-1) = -1;
            Youtput{i,40} = -1;
    end
    % Days to Last Follow-up
    switch Y{i,14}
        case '#N/A'
            daysLastFollowUp(i-1) = -1;
            Youtput{i,14} = -1;
        case '[Not Available]'
            daysLastFollowUp(i-1) = -1;
            Youtput{i,14} = -1;
        otherwise
            daysLastFollowUp(i-1) = Y{i,14};
            Youtput{i,14} = Y{i,14};
    end
    % Days to Death
    switch Y{i,15}
        case '#N/A'
            daysToDeath(i-1) = -1;
            Youtput{i,15} = -1;
        case '[Not Available]'
            daysToDeath(i-1) = -1;
            Youtput{i,15} = -1;
        case '[Not Applicable]'
            daysToDeath(i-1) = -1;
            Youtput{i,15} = -1;
        otherwise
            daysToDeath(i-1) = Y{i,15};
            Youtput{i,15} = Y{i,15};
    end
    % Vital status
    switch Y{i,13}
        case 'Alive'
            survival(i-1) = 1;
            Youtput{i,13} = 1;
        case 'Dead'
            survival(i-1) = 0;
            Youtput{i,13} = 0;
        otherwise % '#N/A'
            survival(i-1) = -1;
            Youtput{i,13} = -1;
    end
    
    % Days to birth
    if (~isnumeric(Y{i,6}))
        Youtput{i,6} = 0;
    end
    
    % Histological type
    switch Y{i,20}
        case 'Thyroid Papillary Carcinoma - Classical/usual'
            Youtput{i,20} = 1;
        case 'Thyroid Papillary Carcinoma - Follicular (>= 99% follicular patterned)'
            Youtput{i,20} = 2;
        case 'Thyroid Papillary Carcinoma - Tall Cell (>= 50% tall cell features)'
            Youtput{i,20} = 3;
        otherwise % 'Other  specify'
            Youtput{i,20} = 4;
    end
    
    % Residual tumor
    switch Y{i,35}
        case 'R0'
            Youtput{i,35} = 0;
        case 'R1'
            Youtput{i,35} = 1;
        case 'R2'
            Youtput{i,35} = 2;
        otherwise % 'RX', '[Not Available]', '[Unknown]'
            Youtput{i,15} = -1;
    end
    
    
end
%% output files
%csvwrite('pathT.csv',pathT);
%csvwrite('pathN.csv',pathN);
%csvwrite('pathM.csv',pathM);

%% output Youtput
fid = fopen('ClinicalAnnot.tsv','wt'); 
[nrows,ncols] = size(Youtput);
for row = 1:nrows
    for col = 1:ncols
        if (isnumeric(Youtput{row,col}))
            fprintf(fid,'%f',Youtput{row,col});
        else
            fprintf(fid,'%s',Youtput{row,col});
        end
        if (col<ncols)
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid)

%% output survival information (newer) use XTumor.csv
survivedDays=max(daysLastFollowUp,daysToDeath);
event=-ones(size(survivedDays,1),1);     % censored = -1 has no survival information
for i = 1:size(survivedDays,1)
    if (survival(i)==1)         % alive
        event(i)=0;             % censored / non-event
        if (survivedDays(i)<0)  % alive, but days not available
            survivedDays(i)=0;
        end
    elseif (survival(i)==0)     % died
        event(i)=1;             % not censored / event
        if (survivedDays(i)<0)  % died, but days not available
            survivedDays(i)=0;
        end
    end
end

csvwrite('event.csv',event);
csvwrite('survivedDays.csv',survivedDays);

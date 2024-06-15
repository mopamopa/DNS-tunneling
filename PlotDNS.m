% 1ma riga di 'DtAQSK_stat.txt': mDt	mA	mQ	vDt	vA	vQ	sDt	sA	sQ	kDt
% kA	kQ	g, quindi:
% cd 'D:\sw analisi tunnel new'

clc;

if(1)      
    I=load('DtAQSK_stat_p2p_no1mariga.txt');    
    if(1) 
        figure; gscatter3(I(:,2),I(:,1),I(:,3),I(:,13)); title('m');
    end;
    if(1)
        figure; gscatter(I(:,1),I(:,2),I(:,13),'br','.+'); xlabel('mDt'); ylabel('mA');
        figure; gscatter(I(:,1),I(:,3),I(:,13),'br','.+'); xlabel('mDt'); ylabel('mQ');    
    end;
    g=I(:,13);
    if(0)
        figure; gscatter3(I(:,2+3),I(:,1+3),I(:,3+3),g); title('v');
        figure; gscatter3(I(:,2+6),I(:,1+6),I(:,3+6),g); title('s');
        figure; gscatter3(I(:,2+9),I(:,1+9),I(:,3+9),g); title('k');
    end;
    if(0)        
        figure; gscatter3(I(:,2+9),I(:,1+9),I(:,3+9),g); title('k');
    end;    
end;
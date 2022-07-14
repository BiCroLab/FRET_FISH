clear all
close all
clc

%% Input information
gene = {'Minar2 cLAD'; 'Grxcr2 cLAD'; '4930426D05Rik cLAD'; ...
    'Hspa9 iLAD'; 'Nars iLAD'; 'Atp5a1 iLAD'; ...
    'Pkm chr9'; 'P4hb chr11'; 'Hsp90ab1 chr17'};


FRET_ALL = [];
for replicate = 1:3
EffFRET_allele = nan(1500, numel(gene));
i = 1;
%% Extract data from file
for g = 1:numel(gene)
    if strcmp(gene{g}, 'Minar2 cLAD')
        filenum = {'iAM746'; 'iAM765'; 'iAM776'};
    elseif strcmp(gene{g}, 'Grxcr2 cLAD')
        filenum = {'iAM747'; 'iAM759'; 'iAM779'};
    elseif strcmp(gene{g}, 'Hspa9 iLAD')
        filenum = {'iAM748'; 'iAM761'; 'iAM778'};
    elseif strcmp(gene{g}, 'Nars iLAD')
        filenum = {'iAM743'; 'iAM758'; 'iAM777'};
    elseif strcmp(gene{g}, '4930426D05Rik cLAD')
        filenum = {'iAM744'; 'iAM757'; 'iAM786'};
    elseif strcmp(gene{g}, 'Atp5a1 iLAD')
        filenum = {'iAM745'; 'iAM762'; 'iAM780'};
    elseif strcmp(gene{g}, 'Pkm chr9')
        filenum = {'iAM773'; 'iAM785'; ''};
    elseif strcmp(gene{g}, 'P4hb chr11')
        filenum = {'iAM774'; 'iAM782'; ''};
    elseif strcmp(gene{g}, 'Hsp90ab1 chr17')
        filenum = {'iAM775'; 'iAM781'; ''};
    end

    
    try
        ID = readtable(['wCentr.out.dilate5.fret_' filenum{replicate} '_autoPick.csv']);
    catch
        disp('file does not exist')
        continue
    end
    
    r = 0;
    D = [];
    A = [];
    FRET = [];
    for Fields = 1:ID.File(end)

        NucleiAuto = unique(ID.Nuclei(ID.File == Fields));
    
        for n = 1:length(NucleiAuto)
            Indx = (ID.File == Fields & ID.Nuclei == NucleiAuto(n));
    
            r = r + 1;
            D(r, :) = unique(ID.Donor(Indx));
            A(r, :) = unique(ID.Acceptor(Indx));
            FRET(r, :) = unique(ID.Fret(Indx));
        end
    end
    EfficiencyFRET = FRET./(FRET+D)*100;
    EffFRET_allele(1:length(EfficiencyFRET), g) = abs(EfficiencyFRET(:, 1)-EfficiencyFRET(:, 2)); 
end
FRET_ALL = [FRET_ALL; EffFRET_allele];
end

figure('units','normalized','outerposition',[0 0 0.3 0.5])
boxplot(FRET_ALL, 'labels',gene)
ylabel('Allelic FRET efficiency Difference (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)


%cd '/Users/anamota/Desktop'
%print -painters -dsvg allele_diff_per_gene.svg


% Variable = table(EfficiencyFRET);
% writetable(Variable, 'ATAC-seq.txt', 'WriteVariableNames',0)
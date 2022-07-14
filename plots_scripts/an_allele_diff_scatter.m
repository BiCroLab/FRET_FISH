clear all
close all
clc

%% Input information
gene = {'Minar2 cLAD'; 'Grxcr2 cLAD'; '4930426D05Rik cLAD'; ...
    'Hspa9 iLAD'; 'Nars iLAD'; 'Atp5a1 iLAD'; ...
    'Pkm chr9'; 'P4hb chr11'; 'Hsp90ab1 chr17'};

figure('units','normalized','outerposition',[0 0 0.6 0.7])
%% Extract data from file
for g = 1:numel(gene)
    FRET_ALL = [];
    Lam_ALL = [];

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

    for replicate = 1:3
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
        Lam = [];
        for Fields = 1:ID.File(end)

            NucleiAuto = unique(ID.Nuclei(ID.File == Fields));

            for n = 1:length(NucleiAuto)
                Indx = (ID.File == Fields & ID.Nuclei == NucleiAuto(n));

                r = r + 1;
                D(r, :) = unique(ID.Donor(Indx));
                A(r, :) = unique(ID.Acceptor(Indx));
                FRET(r, :) = unique(ID.Fret(Indx));
                Lam(r, :) = unique(ID.lamin_dist_norm(Indx));
            end
        end
        EfficiencyFRET = FRET./(FRET+D)*100;
        Lam = abs(Lam(:, 1) - Lam(:, 2));
        FRET_ALL = [FRET_ALL; EfficiencyFRET];
        Lam_ALL = [Lam_ALL; Lam];
    end
subplot(3, 3, g)
scatter(FRET_ALL(:,1), FRET_ALL(:,2), 10, Lam_ALL, 'filled')
h = colorbar;
ylabel(h, 'Lam distance difference norm', 'FontSize', 8)

p = polyfit(FRET_ALL(:,1),FRET_ALL(:,2),1);
lm = fitlm(FRET_ALL(:,1),FRET_ALL(:,2));
text(70, 20, ['R^2 = ' num2str(lm.Rsquared.Ordinary)],'HorizontalAlignment', 'center','FontSize', 8)
text(60, 70, ['P = ' num2str(lm.Coefficients.pValue(2))],'HorizontalAlignment', 'center','FontSize', 8)

A = FRET_ALL(:,1);
B = FRET_ALL(:,2);
[RHOscc,~] = corr(A(~isnan(A)), B(~isnan(B)),'Type','Spearman');
text(60, 60, ['SCC = ' num2str(RHOscc)],'HorizontalAlignment', 'center','FontSize', 8)

hold on
x1 = linspace(min(FRET_ALL(:,1)), max(FRET_ALL(:,1)));
y1 = polyval(p,x1);
plot(x1,y1)

title(gene{g})
xlim([0 80])
ylim([0 80])
text(70, 10, ['n = ' num2str(sum(~isnan(FRET_ALL(:,1))))],'HorizontalAlignment', 'center','FontSize', 8)
ylabel('Allele FRET efficiency (%)', 'Color', 'k');
xlabel('Allele FRET efficiency (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
end




%cd '/Users/anamota/Desktop'
%print -painters -dsvg scatter_allele_LADs.svg


% Variable = table(EfficiencyFRET);
% writetable(Variable, 'ATAC-seq.txt', 'WriteVariableNames',0)
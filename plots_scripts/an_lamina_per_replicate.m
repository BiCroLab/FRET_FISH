clear all
close all
clc


gene = {'Minar2 cLAD'; 'Grxcr2 cLAD'; '4930426D05Rik cLAD'; ...
    'Hspa9 iLAD'; 'Nars iLAD'; 'Atp5a1 iLAD'; ...
    'Pkm chr9'; 'P4hb chr11'; 'Hsp90ab1 chr17'};
%Minar2	chr18 cLAD
%Grxcr2 chr18 cLAD
%4930426D05Rik chr18 cLAD
%Hspa9	chr18 ciLADs
%Atp5a1	chr18 ciLADs
%Nars chr18 ciLADs


i = 0;
for g = 1:numel(gene)
    figure('units','normalized','outerposition',[0 0 0.3 0.6])
    AlleleG = [];
    geneID = [];
    %% Input information
    filenum = {};
 
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
    %% Extract data from file
    for file = 1:numel(filenum)

        if isempty(filenum) || isempty(filenum{file})
            disp(filenum{file})
            disp('file does not exist for the gene')
            continue
        end

        try
            ID = readtable(['wCentr.out.dilate5.fret_' filenum{file} '_autoPick.csv']);
        catch
            disp(filenum{file})
            disp('file does not exist in the folder')
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
        EfficiencyFRET = EfficiencyFRET(:);
        Lam = Lam(:);

        Int = [0 0.25 0.50 0.75 1.05];
        G = [];
        grp = [];
        for q = 1:4
            G = [G; EfficiencyFRET(Lam >= Int(q) & Lam < Int(q+1))];
            grp = [grp; q*ones(length(EfficiencyFRET(Lam >= Int(q) & Lam < Int(q+1))), 1)];
        end
        AlleleG = [AlleleG; G];
        geneID = [geneID; grp];


%         figure
%         counts = [];
%         for up = 1:4
%             subplot(4, 1, up)
%             plotDist(AlleleG(geneID==up));
%             counts(up) = sum(geneID==up);
%             text(50, 0.055, ['n = ' num2str(sum(geneID==up))],'HorizontalAlignment', 'center','FontSize', 8)
%             xlim([0 61])
%             ylim([0 0.06])
%             title([ gene{g} ' | Lam dist: ' num2str(round(Int(up), 2)) '-' num2str(round(Int(up+1), 2)) ])
%         end

        subplot(3, 1, file)
        plotDist(AlleleG)
        text(70, 0.025, ['n = ' num2str(length(AlleleG))],'HorizontalAlignment', 'center','FontSize', 8)
        title([gene{g} ' replicate ' num2str(file)])
    end

    if strcmp(gene{g}, 'Minar2 cLAD')
        print -painters -dsvg Minar2.svg
    elseif strcmp(gene{g}, 'Grxcr2 cLAD')
        print -painters -dsvg Grxcr2.svg
    elseif strcmp(gene{g}, 'Hspa9 iLAD')
        print -painters -dsvg Hspa9.svg
    elseif strcmp(gene{g}, 'Nars iLAD')
        print -painters -dsvg Nars.svg
    elseif strcmp(gene{g}, '4930426D05Rik cLAD')
        print -painters -dsvg 4930426D05Rik.svg
    elseif strcmp(gene{g}, 'Atp5a1 iLAD')
        print -painters -dsvg Atp5a1.svg
    elseif strcmp(gene{g}, 'Pkm chr9')
        print -painters -dsvg Pkm.svg
    elseif strcmp(gene{g}, 'P4hb chr11')
        print -painters -dsvg P4hb.svg
    elseif strcmp(gene{g}, 'Hsp90ab1 chr17')
        print -painters -dsvg Hsp90ab1.svg
    end
    
end




%cd '/Users/anamota/Desktop'
%print -painters -dsvg distribution_per_layer.svg
%print -painters -dsvg counts_per_layer.svg
%print -painters -dsvg boxplot_per_layer.svg
%print -painters -dsvg distribution_per_gene.svg


function plotDist(EfficiencyFRET)
EfficiencyFRET = EfficiencyFRET(~isnan(EfficiencyFRET));

pts = linspace(0,100,100);
[f,xi] = ksdensity(EfficiencyFRET, pts);

plot(xi,f, '-b', 'LineWidth', 3);


ylabel('Probability Density', 'Color', 'k');
xlabel('FRET efficiency (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
end
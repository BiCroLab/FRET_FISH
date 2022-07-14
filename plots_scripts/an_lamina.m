clear all
close all
clc

figure('units','normalized','outerposition',[0 0 1 0.9])
figure('units','normalized','outerposition',[0 0 0.6 0.9])
figure('units','normalized','outerposition',[0 0 0.6 0.9])
figure('units','normalized','outerposition',[0 0 0.6 0.9])
gene = {'Minar2 cLAD'; 'Grxcr2 cLAD'; '4930426D05Rik cLAD'; ...
    'Hspa9 iLAD'; 'Nars iLAD'; 'Atp5a1 iLAD'; ...
    'Pkm chr9'; 'P4hb chr11'; 'Hsp90ab1 chr17'};
%Minar2	chr18 cLAD
%Grxcr2 chr18 cLAD
%4930426D05Rik chr18 cLAD
%Hspa9	chr18 ciLADs
%Atp5a1	chr18 ciLADs
%Nars chr18 ciLADs

ALL = nan(1000, numel(gene));
LAM = nan(1000, numel(gene));
i = 0;
for g = 1:numel(gene)
    AlleleG = [];
    geneID = [];
    L = [];
    F = [];
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
        L = [L; Lam];
        F = [F; EfficiencyFRET];
    end




    figure(1)
    counts = [];
    for up = 1:4
        subplot(4, numel(gene), g+((up-1)*numel(gene)))
        plotDist(AlleleG(geneID==up));
        counts(up) = sum(geneID==up);
        text(50, 0.055, ['n = ' num2str(sum(geneID==up))],'HorizontalAlignment', 'center','FontSize', 8)
        xlim([0 61])
        ylim([0 0.06])
        title([ gene{g} ' | Lam dist: ' num2str(round(Int(up), 2)) '-' num2str(round(Int(up+1), 2)) ])
    end

    figure(2)
    subplot(3, 3, g)

    plot([1 2 3 4], counts, '-b', 'LineWidth', 3)

    title(gene{g})
    ylabel('Counts in each layer', 'Color', 'k');
    xlabel('Concentric layer towards center', 'Color', 'k');


    figure(3)
    subplot(3, 3, g)
    boxplot(AlleleG, geneID, 'labels', {'0 - 25%', '25 - 50%', '50 - 75%','75 - 100%'},'Symbol', '.k', 'Widths',0.7)
    hold on
    for n = 1:3
        plot(n:n+1, [50-n, 50-n], '-k', 'Linewidth', 1)
        A = AlleleG(geneID == n);
        B = AlleleG(geneID == n+1);
        text((n+n+1)/2, 55-n, ['P = ' num2str(round(ranksum(A(~isnan(A)),...
            B(~isnan(B))), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8)
        text(n, 10, ['n = ' num2str(sum(~isnan(A)))],'HorizontalAlignment', 'center','FontSize', 8)
    end
    ylim([5 60])
    text(n+1, 10, ['n = ' num2str(sum(~isnan(B)))],'HorizontalAlignment', 'center','FontSize', 8)
    ylabel('FRET efficiency (%)')
    xlabel('Lamina distance quantiles (pixel)')
    title(gene{g})


figure(4)
subplot(3, 3, g)
plotDist(AlleleG)
text(70, 0.025, ['n = ' num2str(length(AlleleG))],'HorizontalAlignment', 'center','FontSize', 8)
title(gene{g})

ALL(1:length(AlleleG), g) = AlleleG;
LAM(1:length(L), g) = L;
end

figure('units','normalized','outerposition',[0 0 0.5 0.4])
p = violin(ALL, 'xlabel', gene);
ylabel('FRET efficiency (%)')
hold on
for n=1:8
    plot([1 n+1], [70+n*6 70+n*6], '-k', 'Linewidth', 1)
    text((1+n+1)/2, 72+n*6, ...
        ['P = ' num2str(ranksum(ALL(~isnan(ALL(:, 1)), 1),ALL(~isnan(ALL(:, n+1)), n+1)))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
end
ylim([0 100])



figure('units','normalized','outerposition',[0 0 0.5 0.3])
subplot(1,2, 1)
gene_representation(ALL(:,1:3), 'cLAD');
subplot(1,2, 2)
gene_representation(ALL(:,4:7), 'iLAD');


T_1 = array2table(ALL,'VariableNames',gene);
writetable(T_1, 'FRET_eff.txt', 'WriteVariableNames',1)

T_2 = array2table(LAM,'VariableNames',gene);
writetable(T_2, 'Lam.txt', 'WriteVariableNames',1)




figure('units','normalized','outerposition',[0 0 0.5 0.4])
p = violin(LAM, 'xlabel', gene);
ylabel('Lamina distance norm.')
hold on
for n=1:8
    plot([1 n+1], [1+n*0.1 1+n*0.1], '-k', 'Linewidth', 1)
    text((1+n+1)/2, 1.02+n*0.1, ...
        ['P = ' num2str(ranksum(LAM(~isnan(LAM(:, 1)), 1),LAM(~isnan(LAM(:, n+1)), n+1)))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
end
ylim([0 1.5])
%cd '/Users/anamota/Desktop'
%print -painters -dsvg distribution_per_layer.svg
%print -painters -dsvg counts_per_layer.svg
%print -painters -dsvg boxplot_per_layer.svg
%print -painters -dsvg distribution_per_gene.svg
%print -painters -dsvg violin_distribution_per_gene.svg
%print -painters -dsvg violin_distribution_laminas.svg


ALL = reshape(ALL(:,1:6),[], 2);
figure('units','normalized','outerposition',[0 0 0.3 0.4])
p = violin(ALL, 'xlabel', {'cLADs', 'iLADs'});
ylabel('FRET efficiency (%)')
hold on
for n=1:1
    plot([1 n+1], [70+n*6 70+n*6], '-k', 'Linewidth', 1)
    text((1+n+1)/2, 72+n*6, ...
        ['P = ' num2str(ranksum(ALL(~isnan(ALL(:, 1)), 1),ALL(~isnan(ALL(:, n+1)), n+1)))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
end
ylim([0 100])


LAM = reshape(LAM(:,1:6),[], 2);
figure('units','normalized','outerposition',[0 0 0.3 0.4])
p = violin(LAM, 'xlabel', {'cLADs', 'iLADs'});
ylabel('Lamina distance norm.')
hold on
for n=1:1
    plot([1 n+1], [1+n*0.1 1+n*0.1], '-k', 'Linewidth', 1)
    text((1+n+1)/2, 1.02+n*0.1, ...
        ['P = ' num2str(ranksum(LAM(~isnan(LAM(:, 1)), 1),LAM(~isnan(LAM(:, n+1)), n+1)))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
end
ylim([0 1.5])

T_1 = array2table(ALL,'VariableNames',{'cLADs', 'iLADs'});
writetable(T_1, 'FRET_eff_LAD.txt', 'WriteVariableNames',1)

T_2 = array2table(LAM,'VariableNames',{'cLADs', 'iLADs'});
writetable(T_2, 'Lam_LAD.txt', 'WriteVariableNames',1)



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


function gene_representation(EfficiencyFRET, gene)
if sum(~isnan(EfficiencyFRET)) == 0
    return
end

EfficiencyFRET = EfficiencyFRET(~isnan(EfficiencyFRET));

yyaxis left
histogram(EfficiencyFRET, 30,'FaceColor','none', 'LineWidth', 1)
ylabel('# FRET pairs', 'Color', 'k')
set(gca,'YColor', 'k')
text(58, 20, ['n = ' num2str(length(EfficiencyFRET))],'FontSize', 8, 'Color', 'k')

hold on

yyaxis right
[f,xi] = ksdensity(EfficiencyFRET, 0:0.5:100);
plot(xi,f, '-r', 'LineWidth', 1);

xlim([0 85])
set(gca,'linew',1)
xlabel('FRET efficiency (%)', 'Color', 'k')
ylabel('Probability Density', 'Color', 'k')
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
title(gene, 'Color', 'k')
set(gca,'FontSize', 8)
end
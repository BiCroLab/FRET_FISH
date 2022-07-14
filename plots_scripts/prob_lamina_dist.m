clear all
close all
clc

figure('units','normalized','outerposition',[0 0 1 0.9])
gene = {'Minar2 cLAD'; 'Grxcr2 cLAD'; '4930426D05Rik cLAD'; ...
    'Hspa9 iLAD'; 'Nars iLAD'; 'Atp5a1 iLAD'};
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

end



%cd '/Users/anamota/Desktop'
%print -painters -dsvg distribution_per_layer.svg









function plotDist(EfficiencyFRET)
EfficiencyFRET = EfficiencyFRET(~isnan(EfficiencyFRET));

pts = linspace(0,100,100);
[f,xi] = ksdensity(EfficiencyFRET, pts);

th = df_dapiThDialog(xi, f);

LowCondensation = sum(EfficiencyFRET < th)/length(EfficiencyFRET) *100;
HighCondensation = sum(EfficiencyFRET > th)/length(EfficiencyFRET) *100;

area(xi(xi<th), f(xi<th), 'FaceColor', [0.3010 0.7450 0.9330])
hold on
area(xi(xi>th), f(xi>th), 'FaceColor', [0.8500 0.3250 0.0980])
text(15, 0.015, [num2str(round(LowCondensation)) '%'],'FontSize', 8, 'Color', 'k')
text(40, 0.015, [num2str(round(HighCondensation)) '%'],'FontSize', 8, 'Color', 'k')


xlim([0 100])
ylim([0 0.065])

ylabel('Probability Density', 'Color', 'k');
xlabel('FRET efficiency (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
text(80, 0.035, ['n = ' num2str(sum(sum(~isnan(EfficiencyFRET))))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
end



function th = df_dapiThDialog(D, f)
% function th = df_dapiThDialog(D)
% Pick an upper threshold for DAPI
%


if ~exist('th', 'var')
    th = median(D);
end

gui.isok = 0; % set to 1 if 'Ok' was pressed at exit
gui.f = figure('Name', 'Set upper threshold for FRET efficiency', 'NumberTitle', 'off');
gui.a = axes('Position', [0.1,0.25,.8,.7]);
% gui.h = histogram('Parent', gui.a, D, 50);
gui.h = plot(D, f);
hold on
ax = axis();
gui.thLine = plot([th, th], [ax(3), ax(4)], 'LineWidth', 2);

set(gui.f, 'WindowButtonDownFcn', @interStart);

gui.thValue = uicontrol('Style', 'text', ...
    'String', '', ...
    'Units', 'Normalized', ...
    'Position', [0.1,0,.8,.2], ...
    'Callback', @ok, ...
    'Parent', gui.f, ...
    'HorizontalAlignment','left', ...
    'FontName', get(0,'FixedWidthFontName'));

gui.ok = uicontrol('Style', 'pushbutton', ...
    'String', 'Ok', ...
    'Units', 'Normalized', ...
    'Position', [0.85,0.05,.1,.1], ...
    'Callback', @ok, ...
    'Parent', gui.f);

setTh(th);

uiwait(gui.f);

if gui.isok == 0
    th = [];
end

if isvalid(gui.f)
    close(gui.f);
end

function ok(varargin)
    gui.isok = 1;
    uiresume();
end

    function interStart(varargin)
        gco
        if gco == gui.h | gco == gui.a
            x = get(gui.a, 'CurrentPoint'); x = x(1);        
            setTh(x);          
        end
        if gco == gui.thLine
            set(gui.f, 'WindowButtonMotionFcn', @lineDrag);  
            set(gui.f, 'WindowButtonUpFcn', @stopDrag);
        end
    end

    function stopDrag(varargin)
            set(gui.f, 'WindowButtonMotionFcn', []);  
    end

    function lineDrag(varargin)
           x = get(gui.a, 'CurrentPoint'); x = x(1);
           setTh(x);
    end

    function setTh(x)
        gui.thLine.XData = ones(1,2)*x;
        th = x;
        gui.thValue.String = sprintf('Nuclei: %d\nTh: %.2e\nAbove: %d\nBelow: %d', numel(D), th, sum(D>th), sum(D<th));     
    end


end
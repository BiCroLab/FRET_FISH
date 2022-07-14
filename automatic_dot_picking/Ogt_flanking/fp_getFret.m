function T = getFret(s, MM, NN, nn, prefix)
% function O = getFret(s, MM, NN, nn)
% Get the FRET values for the nuclei number nn
%
% Return values:
% O is a struct where Acceptor, Donor, ...

% Extract the nuclei and the associated metadata
N = NN{nn};
M = MM{N.metaNo};

O.('nuclei') = nn; % all nuclei regardless the FoV
O.('File') = [sprintf('%03d', N.metaNo) '.NM']; % FoV
O.('Nuclei') = N.nucleiNr; % Original nuclei

s.donorchannel;
donoridx = find(cellfun(@(v) sum(strfind(v, s.donorchannel)), M.channels)>0);
acceptoridx = find(cellfun(@(v) sum(strfind(v, s.acceptorchannel)), M.channels)>0);

adots = N.userDots{acceptoridx};
ddots = N.userDots{donoridx};

% x,y,z,#, x,y,z,#
if isempty(adots)
    P = [];
else
    P = fp_getPairs(s, adots, ddots);
end

if numel(P) == 0
    T = [];
    return
end

A = []; % Pixel values for the acceptor dots
D = []; % Pixel values for the donor dots
assert(isequal(M.dotsMeta{5}, 'pixel'));
for kk = 1:size(P,1)    
    A(kk, 1) = adots(P(kk, 4), 5); % Col 4 is for DoG and Col 5 is intensity pixel
    D(kk, 1) = ddots(P(kk, 8), 5);
end

O.('n_pairs') = size(P,1);

% save getPairs adots ddots
% Returns rows where each row contains
%    acceptor    donor
% [x1, y1, z1, x1, x2, y2]

% Use mean position of each pair
%
% Extract the ending '007.tif' or similar
ending = regexp(MM{N.metaNo}.dapifile, '[0-9][0-9][0-9]\.tiff?', 'match');
ending = ending{1};

folder = s.imagefolder; prestring = prefix;
AcceptorFile = sprintf('%s%s%s_%s', folder, prestring, s.acceptorchannel, ending);
assert(isfile(AcceptorFile));
DonorFile = sprintf('%s%s%s_%s', folder, prestring, s.donorchannel, ending);
assert(isfile(DonorFile));

O.('tmrFile') = DonorFile;
O.('Cy5File') = AcceptorFile;

s.donorfile = DonorFile;
s.acceptorfile = AcceptorFile;



T = table();
for kk = 1:numel(A)
    row = O;
    row.('X_Cy5') = P(kk, 1);
    row.('Y_Cy5') = P(kk, 2);
    row.('Z_Cy5') = P(kk, 3);
    row.('X_tmr') = P(kk, 5); % X_D
    row.('Y_tmr') = P(kk, 6); % Y_D
    row.('Z_tmr') = P(kk, 7); % Z_D
    row.('Pair') = kk;
    T = [T; struct2table(row)];
end

end
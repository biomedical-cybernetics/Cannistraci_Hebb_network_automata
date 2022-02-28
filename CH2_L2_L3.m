function scores = CH2_L2_L3(x, L, w, par)

% code to compute the Cannistraci-Hebb (CH) network automata scores
% for network links considering paths of length two (L2) or three (L3)
%
% Authors:
% Alessandro Muscoloni, 2020-01-16
%
% Reference:
% "Local-community network automata modelling based on length-three-paths
% for prediction of complex network structures in protein interactomes, food webs and more"
% A. Muscoloni, I. Abdelhamid, C. V. Cannistraci, bioRxiv, 2018
% https://doi.org/10.1101/346916
%
% Released under MIT License
% Copyright (c) 2020 A. Muscoloni, C. V. Cannistraci

%%% INPUT %%%
% x - adjacency matrix of the network;
%     the network is considered unweighted, undirected and zero-diagonal
%
% L - [optional] integer to indicate the path length to compute:
%     0 -> compute both CH2-L2 and CH2-L3
%     2 -> compute only CH2-L2
%     3 -> compute only CH2-L3
%     if not given or empty, the option 0 is considered
%
% w - [optional] 2-columns matrix (id1,id2) indicating the links for which the score should be calculated;
%     if not given or empty, the scores for all the missing links are computed
%
% par - [optional] 1 or 0 to indicate whether the function should use parallel computation or not;
%     if not given or empty, parallel computation is used
%
%%% OUTPUT %%%
% scores - 3-columns or 4-columns matrix depending on the input parameter L:
%     L=0 -> 4-columns matrix containing the values (id1,id2,score_CH2_L2,score_CH2_L3)
%     L=2 -> 3-columns matrix containing the values (id1,id2,score_CH2_L2)
%     L=3 -> 3-columns matrix containing the values (id1,id2,score_CH2_L3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
narginchk(1,4)
validateattributes(x, {'logical','numeric'}, {'square','binary'});
x = double(sparse(x));
x = max(x,x');
x(speye(size(x))==1) = 0;
n = size(x,1);
if ~exist('L', 'var') || isempty(L)
    L = 0;
elseif ~any(L==[0,2,3])
    error('Possible values for L are: [0,2,3]')
end
if ~exist('w', 'var') || isempty(w)
    [i,j] = find(triu(x==0,1));
    w = [i,j]; clear i j;
else
    validateattributes(w, {'numeric'}, {'2d','ncols',2,'integer','>=',1,'<=',n});
end
if ~exist('par', 'var') || isempty(par)
    par = Inf;
else
    validateattributes(par, {'numeric'}, {'scalar','binary'});
    if par == 1
        par = Inf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
m = size(w,1);
d = full(sum(x,1));
if any(L==[0,2])
    scores_CH2_L2 = zeros(m,1);
end
if any(L==[0,3])
    scores_CH2_L3 = zeros(m,1);
end

% adjacency list
A = cell(n,1);
parfor (i = 1:n, par)
    A{i} = find(x(i,:));
end

% main code
parfor (i = 1:m, par)
    u = w(i,1); v = w(i,2);
    if x(u,v)==0
        Au = A{u}; Av = A{v};
    else
        Au = setdiff(A{u},v); Av = setdiff(A{v},u);
    end
    
    % L2
    cn = intersect(Au,Av);
    if any(L==[0,2]) && ~isempty(cn)
        di = full(sum(x(cn,cn),1));
        de = d(cn) - di - 2;
        scores_CH2_L2(i) = sum((1+di)./(1+de));
    end
    
    % L3
    [e1,e2] = find(x(Au,Av));
    paths = [Au(e1); Av(e2)];
    if any(L==[0,3]) && ~isempty(paths)
        paths_size = size(paths);
        paths = reshape(paths, 1, numel(paths));
        [cn,~,idx] = unique(paths);
        di = full(sum(x(cn,cn),1));
        di = di(idx);
        de = d(paths) - di - full(x(u,paths)) - full(x(v,paths));
        di = reshape(di, paths_size);
        de = reshape(de, paths_size);
        scores_CH2_L3(i) = sum(geomean(1+di,1)./geomean(1+de,1));
    end
end

% output scores
if L==0
    scores = [w scores_CH2_L2 scores_CH2_L3];
elseif L==2
    scores = [w scores_CH2_L2];
elseif L==3
    scores = [w scores_CH2_L3];
end

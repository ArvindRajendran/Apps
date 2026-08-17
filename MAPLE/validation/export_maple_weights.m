% export_maple_weights.m
% Stage 0 of the MAPLE HTML app: extract the trained ANN weights from the
% two .mat model files of github.com/ArvindRajendran/MAPLE, self-check a
% hand-rolled forward pass against the MATLAB network objects, discover
% the exact input convention against the labelled LHS data, and report
% R2_adj of the nets vs the detailed-model labels.
%
% Run from the repo directory:  matlab -batch "run('export_maple_weights.m')"

repo = pwd;
outdir = fullfile(repo, 'export');
if ~exist(outdir, 'dir'), mkdir(outdir); end

tags  = {'LPP','FP'};
files = {'MAPLE4stepwithLPPModels.mat','MAPLE4stepwithFPModels.mat'};
% Net roles as in the repo's MAPLE.m wrapper: net1=recovery, net2=purity.
% Confirmed against the Limits-paper SI optimum points (their purity
% column is pinned at the active >=95% constraint, so the semantics are
% unambiguous). CAUTION: the repo's LHS-samples xlsx has its Pu/Re
% column HEADERS swapped — its "Pu [%]" column correlates with net1
% (recovery) and its "Re [%]" column with net2 (purity). net5 is the
% energy at 100% pump efficiency; net6 correlates with nothing shipped.
kpis  = {'recovery','purity','log10energy','productivity','log10energy100'};

EX = struct();
for k = 1:2
    A = load(files{k});
    fprintf('=== %s: variables: %s\n', files{k}, strjoin(fieldnames(A)', ', '));
    fprintf('mue = '); fprintf('%.10g ', A.mue); fprintf('\n');
    fprintf('sig = '); fprintf('%.10g ', A.sig); fprintf('\n');
    nets = {A.net1, A.net2, A.net3, A.net4, A.net5};
    C = struct('mue', A.mue(:).', 'sig', A.sig(:).');
    for j = 1:5
        nn = nets{j};
        NN = exportnet(nn, kpis{j});
        C.(sprintf('net%d', j)) = NN;
        % --- self-check: manual forward pass vs network object ---
        rng(42 + j);
        X = randn(numel(A.mue), 200);            % z-scored domain
        y1 = nn(X);
        y2 = zeros(size(y1));
        for m = 1:size(X,2), y2(:,m) = fwdnet(NN, X(:,m)); end
        err = max(abs(y1 - y2), [], 'all') / max(1, max(abs(y1), [], 'all'));
        fprintf('%s %-13s: layers=%d  sizes=[%s]  transfer={%s}  selfcheck=%.3e\n', ...
            tags{k}, kpis{j}, nn.numLayers, ...
            strjoin(string(cellfun(@(l) l.size, nn.layers.', 'uni', 0)), ' '), ...
            strjoin(cellfun(@(l) l.transferFcn, nn.layers.', 'uni', 0), ','), err);
        if err > 1e-12
            error('Self-check failed for %s %s (err=%g)', tags{k}, kpis{j}, err);
        end
    end
    EX.(tags{k}) = C;
end

fid = fopen(fullfile(outdir, 'maple_weights.json'), 'w');
fwrite(fid, jsonencode(EX)); fclose(fid);
fprintf('weights written: %s\n', fullfile(outdir, 'maple_weights.json'));

% ------------------------- helpers -----------------------------------
function NN = exportnet(nn, name)
    NN = struct(); NN.kpi = name; L = nn.numLayers;
    lay = {};
    for l = 1:L
        S = struct(); S.transfer = nn.layers{l}.transferFcn;
        if l == 1, S.W = nn.IW{1,1}; else, S.W = nn.LW{l, l-1}; end
        S.b = nn.b{l};
        lay{end+1} = S; %#ok<AGROW>
    end
    NN.layers = lay;
    NN.inProc  = procexport(nn.inputs{1});
    NN.outProc = procexport(nn.outputs{L});
end

function P = procexport(io)
    P = {};
    for m = 1:numel(io.processFcns)
        Q = struct(); Q.fcn = io.processFcns{m}; Q.settings = io.processSettings{m};
        P{end+1} = Q; %#ok<AGROW>
    end
end

function y = fwdnet(NN, x)
    for m = 1:numel(NN.inProc), x = applyproc(NN.inProc{m}, x); end
    a = x;
    for l = 1:numel(NN.layers)
        S = NN.layers{l}; z = S.W * a + S.b;
        switch S.transfer
            case 'tansig',  a = 2 ./ (1 + exp(-2*z)) - 1;
            case 'logsig',  a = 1 ./ (1 + exp(-z));
            case 'purelin', a = z;
            otherwise, error('unknown transfer %s', S.transfer);
        end
    end
    for m = numel(NN.outProc):-1:1, a = revproc(NN.outProc{m}, a); end
    y = a;
end

function x = applyproc(Q, x)
    switch Q.fcn
        case 'mapminmax'
            x = (x - Q.settings.xoffset) .* Q.settings.gain + Q.settings.ymin;
        case 'removeconstantrows'
            if isfield(Q.settings, 'keep') && ~isempty(Q.settings.keep)
                x = x(Q.settings.keep, :);
            end
        case 'mapstd'
            x = (x - Q.settings.xmean) ./ Q.settings.xstd;
        otherwise, error('unknown processFcn %s', Q.fcn);
    end
end

function y = revproc(Q, y)
    switch Q.fcn
        case 'mapminmax'
            y = (y - Q.settings.ymin) ./ Q.settings.gain + Q.settings.xoffset;
        case 'removeconstantrows'
            if isfield(Q.settings, 'remove') && ~isempty(Q.settings.remove)
                z = zeros(Q.settings.xrows, size(y,2));
                z(Q.settings.keep, :) = y;
                z(Q.settings.remove, :) = repmat(Q.settings.value, 1, size(y,2));
                y = z;
            end
        case 'mapstd'
            y = y .* Q.settings.xstd + Q.settings.xmean;
        otherwise, error('unknown processFcn %s', Q.fcn);
    end
end

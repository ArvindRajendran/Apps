% net_role_grid.m — identify which net predicts which label, and pin the
% exact input transform, by brute R2 grid.
xf = 'MAPLE_LHCSamples_with_Labels_LPP_FP_Cycles.xlsx';
files = {'MAPLE4stepwithLPPModels.mat','MAPLE4stepwithFPModels.mat'};
sheets = {'SamplesLPP','SamplesFP'};

for k = 1:2
    A = load(files{k});
    T = readmatrix(xf, 'Sheet', sheets{k});
    fprintf('=== %s size %dx%d; row1: ', sheets{k}, size(T,1), size(T,2));
    fprintf('%.4g ', T(1,1:min(18,size(T,2)))); fprintf('\n');
    T = T(all(isfinite(T(:,1:16)), 2), :);
    fprintf('finite rows: %d\n', size(T,1));

    for Tref = [298 298.15]
        N = T(:,1:12);
        N(:,1) = N(:,1) .* N(:,6);                       % qsat*rho [mol/m3]
        N(:,2) = log10(T(:,2) .* exp(-T(:,4)*1000/(8.314*Tref)));
        N(:,3) = log10(T(:,3) .* exp(-T(:,5)*1000/(8.314*Tref)));
        N(:,10) = log10(N(:,10)./N(:,9));
        N(:,11) = log10(N(:,11)./N(:,9));
        Z = ((N - A.mue)./A.sig)';
        % candidate labels
        lab = struct('Pu',T(:,13), 'Re',T(:,14), 'En',T(:,15), ...
                     'log10En',log10(max(T(:,15),1e-6)), 'Prod',T(:,16));
        if size(T,2) >= 17 && any(isfinite(T(:,17)))
            lab.En100 = T(:,17); lab.log10En100 = log10(max(T(:,17),1e-6));
        end
        nets = {A.net1 A.net2 A.net3 A.net4 A.net5 A.net6};
        fn = fieldnames(lab);
        fprintf('--- Tref=%.2f: R2 grid (rows=net, cols=%s)\n', Tref, strjoin(fn', ' '));
        for j = 1:6
            p = nets{j}(Z)';
            fprintf('net%d: ', j);
            for m = 1:numel(fn)
                L = lab.(fn{m});
                ok = isfinite(L);
                r2 = 1 - sum((L(ok)-p(ok)).^2)/sum((L(ok)-mean(L(ok))).^2);
                fprintf('%8.4f ', r2);
            end
            fprintf('\n');
        end
    end
end

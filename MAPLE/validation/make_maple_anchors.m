% make_maple_anchors.m — generate validation artifacts for the JS engine:
%   maple_anchors.json  : 200 LHS rows/cycle + 4 presets x 3 operating
%                         conditions, with all 6 net outputs (MATLAB truth)
%   maple_labels_LPP.json / maple_labels_FP.json : all labelled rows
%                         (12 raw inputs + Pu,Re,En,Prod,En100 labels)
%   R2 reference values printed for the harness README.
% Net roles (per the repo wrapper, confirmed vs the Limits-paper SI
% optima where purity is pinned by the >=95% constraint):
%   net1=recovery  net2=purity  net3=log10(En)  net4=productivity
%   net5=log10(En at 100% pump efficiency)  net6=probed below
% The LHS xlsx label columns have SWAPPED Pu/Re headers: its col 13
% ("Pu") is recovery, col 14 ("Re") is purity. Everything below uses the
% corrected semantics.
xf = 'MAPLE_LHCSamples_with_Labels_LPP_FP_Cycles.xlsx';
files = {'MAPLE4stepwithLPPModels.mat','MAPLE4stepwithFPModels.mat'};
sheets = {'SamplesLPP','SamplesFP'};
tags = {'LPP','FP'};
outdir = fullfile(pwd, 'export');

% presets: [qsat(mol/kg) b0CO2 b0N2 dUCO2 dUN2 rho]  (dU negative, kJ/mol)
pnames = {'Zeolite13X','UTSA16','IISERPMOF2','MgMOF74'};
pdata = [4.390  2.50e-6 2.70e-6 -31.19 -16.38 1130;
         4.478  4.70e-7 1.40e-6 -30.57  -9.91 1000;
         5.000  2.02e-7 2.64e-7 -31.13 -11.89 1000;
         5331.29/1130 6.38e-7 2.06e-6 -33.73 -18.32 1130];
% operating conditions: [y tads PH PI PL vF]
ocnames = {'paper1','mid','highy'};
ocdata = [0.15 92.4 1.0 0.08 0.03 0.64;
          0.15 60.0 2.0 0.30 0.05 0.80;
          0.30 40.0 1.5 0.20 0.04 0.50];

OUT = struct();
for k = 1:2
    A = load(files{k});
    nets = {A.net1 A.net2 A.net3 A.net4 A.net5 A.net6};
    T = readmatrix(xf, 'Sheet', sheets{k});
    T = T(all(isfinite(T(:,1:16)), 2), :);
    n = size(T,1);

    % --- probe net6 against remaining columns (times?) ---
    Z = zin(T(:,1:12), A);
    p6 = nets{6}(Z)';
    fprintf('%s net6 vs extra cols: ', tags{k});
    for c = 17:size(T,2)
        L = T(:,c); ok = isfinite(L);
        if nnz(ok) > 100
            r2 = 1 - sum((L(ok)-p6(ok)).^2)/sum((L(ok)-mean(L(ok))).^2);
            fprintf('col%d:%.4f ', c, r2);
            Ll = log10(max(L(ok),1e-9));
            r2l = 1 - sum((Ll-p6(ok)).^2)/sum((Ll-mean(Ll)).^2);
            fprintf('log:%.4f  ', r2l);
        end
    end
    fprintf('\n');

    % --- R2 reference values over all labelled rows ---
    % (label col 13 carries the swapped "Pu" header but IS recovery;
    %  col 14 carries "Re" but IS purity — matched accordingly)
    prd = [nets{2}(Z)' nets{1}(Z)' 10.^(nets{3}(Z))' nets{4}(Z)' 10.^(nets{5}(Z))'];
    lbl = [T(:,14) T(:,13) T(:,15) T(:,16) T(:,17)];
    kn = {'purity','recovery','energy','productivity','energy100'};
    R2 = zeros(1,5);
    for j = 1:5
        ok = isfinite(lbl(:,j));
        R2(j) = 1 - sum((lbl(ok,j)-prd(ok,j)).^2)/sum((lbl(ok,j)-mean(lbl(ok,j))).^2);
        fprintf('%s R2(%s) = %.6f over %d rows\n', tags{k}, kn{j}, R2(j), nnz(ok));
    end

    % --- 200 anchor rows, deterministic spread ---
    idx = round(linspace(1, n, 200));
    S = struct('inputs', T(idx,1:12), 'outputs', struct( ...
        'purity', prd(idx,1), 'recovery', prd(idx,2), 'energy', prd(idx,3), ...
        'productivity', prd(idx,4), 'energy100', prd(idx,5), 'net6', p6(idx)));

    % --- presets x operating conditions ---
    P = {};
    for a = 1:4
        for b = 1:3
            row = [pdata(a,1) pdata(a,2) pdata(a,3) pdata(a,4) pdata(a,5) pdata(a,6) ...
                   ocdata(b,1) ocdata(b,2) ocdata(b,3) ocdata(b,4) ocdata(b,5) ocdata(b,6)];
            Zi = zin(row, A);
            E = struct('name', sprintf('%s_%s', pnames{a}, ocnames{b}), ...
                'inputs', row, ...
                'purity', nets{2}(Zi), 'recovery', nets{1}(Zi), ...
                'energy', 10.^(nets{3}(Zi)), 'productivity', nets{4}(Zi), ...
                'energy100', 10.^(nets{5}(Zi)));
            P{end+1} = E; %#ok<AGROW>
        end
    end
    OUT.(tags{k}) = struct('samples', S, 'presets', {P}, 'R2', R2);

    % --- full label file for the JS-side R2 reproduction ---
    L = struct('cols', {{'qsat_molkg','b0CO2','b0N2','dUCO2','dUN2','rho', ...
        'yF','tads','PH','PI','PL','vF','Pu','Re','En','Prod','En100'}}, ...
        'rows', T(:,1:17));
    fid = fopen(fullfile(outdir, sprintf('maple_labels_%s.json', tags{k})), 'w');
    fwrite(fid, jsonencode(L)); fclose(fid);
end
fid = fopen(fullfile(outdir, 'maple_anchors.json'), 'w');
fwrite(fid, jsonencode(OUT)); fclose(fid);
fprintf('anchors + labels written to %s\n', outdir);

function Z = zin(N, A)
    % physical inputs -> z-scored net inputs
    % N cols: qsat[mol/kg] b0CO2 b0N2 dUCO2 dUN2 rho yF tads PH PI PL vF
    M = N;
    M(:,1) = N(:,1) .* N(:,6);
    M(:,2) = log10(N(:,2) .* exp(-N(:,4)*1000/(8.314*298)));
    M(:,3) = log10(N(:,3) .* exp(-N(:,5)*1000/(8.314*298)));
    M(:,10) = log10(N(:,10) ./ N(:,9));
    M(:,11) = log10(N(:,11) ./ N(:,9));
    Z = ((M - A.mue) ./ A.sig)';
end

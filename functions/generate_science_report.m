function report_path = generate_science_report(results, config)
% GENERATE_SCIENCE_REPORT  Plain-text discovery report for NEUROSCOPE.
%
% Replaces generate_clinical_report() with a science-first narrative that
% frames findings as communication biology, not therapeutic targets.
%
% Sections:
%   1. Dataset overview
%   2. Channel capacity landscape (novel vs. existing tool comparison)
%   3. Top differential channels (strength + capacity tests)
%   4. Capacity-only enriched channels (the novel finding)
%   5. Causal validation summary (IV/2SLS)
%   6. Communication programs (NMF)
%   7. Spectral communication communities
%   8. Communication states vs. transcriptomic clusters
%   9. Downstream response profiling
%  10. Methods note (for paper writing)

out_dir = config.output_dir;
if ~exist(out_dir, 'dir'), out_dir = '.'; end
report_path = fullfile(out_dir, [config.dataset_name '_science_report.txt']);

fid = fopen(report_path, 'w');
if fid < 0
    warning('generate_science_report: could not open %s for writing', report_path);
    return;
end

cs   = results.condition_specific;
net  = results.condition_network;
cr   = results.channel_rankings;
fdr  = config.fdr_threshold;

line  = @(s) fprintf(fid, '%s\n', s);
blank = @()  fprintf(fid, '\n');
hdr   = @(s) fprintf(fid, '\n%s\n%s\n', s, repmat('=', 1, length(s)));
sub   = @(s) fprintf(fid, '\n%s\n%s\n', s, repmat('-', 1, length(s)));

%% Header
fprintf(fid, 'NEUROSCOPE v5.0 — Scientific Discovery Report\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Dataset:   %s\n', config.dataset_name);
fprintf(fid, 'Condition: %s\n', config.condition_name);
fprintf(fid, 'Reference: %s  (empty = one-sample mode)\n', config.reference_name);
fprintf(fid, 'FDR cutoff: %.2f\n', fdr);
line(repmat('=',1,60));


%% Spearman fallback warning (shown prominently when KSG was not used)
using_spearman = false;
if ~isempty(net) && isfield(net, 'mi_mode') && ~isempty(net.mi_mode)
    using_spearman = any(strcmp(net.mi_mode, 'spearman_copula'));
end

if using_spearman
    blank();
    fprintf(fid, '*** DATA QUALITY WARNING: SPEARMAN FALLBACK ACTIVE ***\n');
    fprintf(fid, '%s\n', repmat('*', 1, 60));
    fprintf(fid, '\n');
    fprintf(fid, 'This run used the SPEARMAN COPULA MI approximation, not\n');
    fprintf(fid, 'the KSG k-NN estimator. Fewer than 10 donors were present.\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Consequences:\n');
    fprintf(fid, '\n');
    fprintf(fid, '  1. CondCapacity values are INFLATED and uninterpretable.\n');
    fprintf(fid, '     Spearman MI saturates to ~1.0 for broadly expressed\n');
    fprintf(fid, '     genes regardless of whether real information is\n');
    fprintf(fid, '     transmitted (e.g. CADM1 capacity=1.0 is an artifact).\n');
    fprintf(fid, '\n');
    fprintf(fid, '  2. The headline claim that "X%% of co-expressed pairs\n');
    fprintf(fid, '     transmit I~=0" CANNOT be made from this run.\n');
    fprintf(fid, '     With Spearman fallback, capacity tracks expression\n');
    fprintf(fid, '     strength almost perfectly — capacity-only enriched\n');
    fprintf(fid, '     channels will always be near zero as a result.\n');
    fprintf(fid, '\n');
    fprintf(fid, '  3. FDRCapacity and FDR tests are not independent here.\n');
    fprintf(fid, '     Do not report FDRCapacity results as a novel finding.\n');
    fprintf(fid, '\n');
    fprintf(fid, 'What IS valid from this run:\n');
    fprintf(fid, '  - Network topology (which L-R pairs are expressed)\n');
    fprintf(fid, '  - Communication hub identification\n');
    fprintf(fid, '  - NMF communication programs\n');
    fprintf(fid, '  - Spectral community structure\n');
    fprintf(fid, '  - IV/2SLS causal validation results\n');
    fprintf(fid, '  - KSG AUROC on synthetic validation (unaffected by\n');
    fprintf(fid, '    donor count; this IS the publishable number)\n');
    fprintf(fid, '\n');
    fprintf(fid, 'To obtain publication-quality KSG MI estimates, re-run\n');
    fprintf(fid, 'on a dataset with >=20 donors (e.g. OneK1K: 982 donors).\n');
    fprintf(fid, '%s\n', repmat('*', 1, 60));
    blank();
end


%% Section 1 — Network overview
hdr('1. COMMUNICATION NETWORK OVERVIEW');
if ~isempty(net) && isfield(net, 'weights')
    fprintf(fid, 'Total L-R interactions (channel capacity > 0):  %d\n', length(net.weights));
    fprintf(fid, 'Unique ligands:                                  %d\n', length(unique(net.ligands)));
    fprintf(fid, 'Unique receptors:                                %d\n', length(unique(net.receptors)));
    if isfield(net,'channel_capacity') && ~isempty(net.channel_capacity)
        cap = net.channel_capacity;
        fprintf(fid, 'Mean channel capacity:                           %.4f bits\n', mean(cap));
        fprintf(fid, 'Median channel capacity:                         %.4f bits\n', median(cap));
        fprintf(fid, 'High-capacity channels (>0.1 bits):              %d (%.1f%%)\n', ...
            sum(cap > 0.1), 100*mean(cap > 0.1));
        fprintf(fid, 'High-capacity channels (>0.5 bits):              %d (%.1f%%)\n', ...
            sum(cap > 0.5), 100*mean(cap > 0.5));
    end
else
    line('Network data unavailable.');
end

%% Section 2 — Differential channels
hdr('2. DIFFERENTIAL CHANNEL ANALYSIS');
if ~isempty(cs) && height(cs) > 0
    n_sig_str = sum(cs.FDR < fdr);
    n_sig_cap = sum(cs.FDRCapacity < fdr);
    n_novel   = sum(cs.FDRCapacity < fdr & cs.FDR >= fdr);
    fprintf(fid, 'Total L-R pairs tested:                          %d\n', height(cs));
    fprintf(fid, 'Significant by strength (Wilcoxon, FDR<%.2f):   %d\n', fdr, n_sig_str);
    fprintf(fid, 'Significant by capacity (KSG MI, FDR<%.2f):     %d  [NOVEL TEST]\n', fdr, n_sig_cap);
    fprintf(fid, 'Capacity-only enriched (missed by strength):     %d  [KEY FINDING]\n\n', n_novel);

    if using_spearman
        fprintf(fid, 'NOTE: Capacity-only = 0 is an ARTIFACT of Spearman\n');
        fprintf(fid, 'fallback, not a scientific finding. Spearman MI tracks\n');
        fprintf(fid, 'expression strength almost perfectly so capacity-only\n');
        fprintf(fid, 'enriched channels are always near zero. Re-run with\n');
        fprintf(fid, '>=20 donors for genuine KSG-based results.\n');
    else
    fprintf(fid, 'Interpretation: Capacity-only enriched channels show elevated\n');
    fprintf(fid, 'information transmission without proportional expression increase.\n');
    fprintf(fid, 'These would be INVISIBLE to CellChat/NicheNet/CellPhoneDB.\n');
    end

    has_ci = ismember('BootCI_lo', cr.Properties.VariableNames) && ...
             ismember('BootCI_hi', cr.Properties.VariableNames);

    sub('Top 10 Channels by Channel Capacity');
    if has_ci
        fprintf(fid, '%-20s %-20s %8s %14s %8s %8s %20s\n', ...
            'Ligand','Receptor','Capacity','95%CI','FDR_cap','FDR_str','Pathway');
        fprintf(fid, '%s\n', repmat('-',1,115));
    else
        fprintf(fid, '%-20s %-20s %8s %8s %8s %20s %20s\n', ...
            'Ligand','Receptor','Capacity','FDR_cap','FDR_str','ChannelType','Pathway');
        fprintf(fid, '%s\n', repmat('-',1,110));
    end
    n_print = min(10, height(cr));
    for i = 1:n_print
        if has_ci && ~isnan(cr.BootCI_lo(i)) && cr.BootCI_lo(i) ~= cr.BootCI_hi(i)
            ci_str = sprintf('[%.3f,%.3f]', cr.BootCI_lo(i), cr.BootCI_hi(i));
            fprintf(fid, '%-20s %-20s %8.4f %14s %8.4f %8.4f %20s\n', ...
                cr.Ligand{i}, cr.Receptor{i}, cr.CondCapacity(i), ci_str, ...
                cr.FDRCapacity(i), cr.FDR(i), cr.Pathway{i});
        else
            fprintf(fid, '%-20s %-20s %8.4f %8.4f %8.4f %20s %20s\n', ...
                cr.Ligand{i}, cr.Receptor{i}, cr.CondCapacity(i), ...
                cr.FDRCapacity(i), cr.FDR(i), ...
                cr.ChannelType{i}, cr.Pathway{i});
        end
    end
    if has_ci
        n_wide = sum(cr.BootCI_hi(1:n_print) - cr.BootCI_lo(1:n_print) > 0.15, 'omitnan');
        if n_wide > 0
            fprintf(fid, 'NOTE: %d/%d top channels have 95%% CI width >0.15 bits — expected at n=%d\n', ...
                n_wide, n_print, net.n_samples_used);
        end
    end

    % Channel type distribution
    sub('Channel Type Distribution');
    types = {'HighCapacity_Broadcast','HighCapacity_Selective',...
             'LowCapacity_Broadcast','LowCapacity_Selective'};
    for t = 1:length(types)
        n_t = sum(strcmp(cs.ChannelType, types{t}));
        fprintf(fid, '  %-30s  %d (%.1f%%)\n', types{t}, n_t, 100*n_t/height(cs));
    end

    % Pathway distribution
    sub('Pathway Annotation (top 5)');
    [pw, ~, pidx] = unique(cs.Pathway);
    pw_counts = accumarray(pidx, 1);
    [pw_sorted, pw_ord] = sort(pw_counts, 'descend');
    for p = 1:min(5, length(pw))
        fprintf(fid, '  %-25s  %d channels\n', pw{pw_ord(p)}, pw_sorted(p));
    end
else
    line('No differential channel data available.');
end

%% Section 3 — Causal validation
hdr('3. CAUSAL INTERACTION VALIDATION (IV/2SLS)');
causal = results.causal_results;
if isfield(causal,'n_causal') && isfield(causal,'n_tested') && causal.n_tested > 0
    fprintf(fid, 'Instrument type:        %s\n', causal.instrument_type);
    fprintf(fid, 'Interactions tested:    %d\n', causal.n_tested);
    fprintf(fid, 'Causally validated:     %d (%.1f%%)\n', ...
        causal.n_causal, 100*causal.causal_fraction);
    blank();
    fprintf(fid, 'Interpretation: Causal validation uses donor_id as a natural\n');
    fprintf(fid, 'instrumental variable (analogous to Mendelian Randomization).\n');
    fprintf(fid, 'A causally validated interaction means: natural inter-donor variation\n');
    fprintf(fid, 'in ligand level predicts receptor-response variation in receiver cells,\n');
    fprintf(fid, 'independently of shared confounders (batch, cell state, etc.).\n');
else
    line('Causal validation data unavailable or not run.');
end

%% Section 4 — NMF programs
hdr('4. LATENT COMMUNICATION PROGRAMS (NMF)');
progs = results.comm_programs;
if isfield(progs,'rank') && progs.rank > 0
    fprintf(fid, 'NMF rank selected:  %d programs\n', progs.rank);
    if isfield(progs,'program_table') && istable(progs.program_table) && height(progs.program_table) > 0
        fprintf(fid, '\nTop L-R pairs per program:\n');
        for p = 1:progs.rank
            prog_rows = progs.program_table(progs.program_table.ProgramID == p, :);
            if height(prog_rows) > 0
                if ismember('TopLRPairs', prog_rows.Properties.VariableNames) && height(prog_rows) > 0
                    raw_pairs = strsplit(prog_rows.TopLRPairs{1}, ';');
                    raw_pairs = strtrim(raw_pairs);
                    pairs = strjoin(raw_pairs(1:min(3,length(raw_pairs))), ', ');
                else
                    pairs = '(no pairs)';
                end
                fprintf(fid, '  Program %d: %s\n', p, pairs);
            end
        end
    end
    blank();
    fprintf(fid, 'Interpretation: NMF decomposes the [cluster-pair × L-R] communication\n');
    fprintf(fid, 'tensor into co-regulated interaction modules. Each program represents\n');
    fprintf(fid, 'a coordinated signaling circuit operating across multiple cell types.\n');
else
    line('Communication program data unavailable or not run.');
end

%% Section 5 — Spectral communities
hdr('5. SPECTRAL COMMUNICATION COMMUNITIES');
entropy = results.entropy_results;
if isfield(entropy,'n_communities') && entropy.n_communities > 0
    fprintf(fid, 'Communities detected:  %d\n', entropy.n_communities);
    fprintf(fid, 'Spectral gap:          %.4f  (larger = more modular)\n', entropy.spectral_gap);
    if isfield(entropy,'cluster_roles') && istable(entropy.cluster_roles) && height(entropy.cluster_roles) > 0
        blank();
        fprintf(fid, '%-12s %-25s %8s %8s\n', 'Cluster','Role','H_send','H_recv');
        fprintf(fid, '%s\n', repmat('-',1,60));
        for r = 1:min(10, height(entropy.cluster_roles))
            fprintf(fid, '%-12s %-25s %8.4f %8.4f\n', ...
                num2str(entropy.cluster_roles.Cluster(r)), ...
                entropy.cluster_roles.Role{r}, ...
                entropy.cluster_roles.H_sending(r), ...
                entropy.cluster_roles.H_receiving(r));
        end
    end
    blank();
    fprintf(fid, 'Interpretation: Spectral decomposition of the communication Laplacian\n');
    fprintf(fid, 'reveals functional modules — groups of cell types that communicate\n');
    fprintf(fid, 'more within the group than across groups.\n');
else
    line('Entropy/community data unavailable or not run.');
end

%% Section 6 — Communication identity
hdr('6. COMMUNICATION IDENTITY vs. TRANSCRIPTOMIC IDENTITY');
prof = results.cell_comm_profiles;
id   = [];
if isfield(results, 'identity_results'), id = results.identity_results; end

if isfield(prof,'n_comm_clusters') && prof.n_comm_clusters > 0
    n_trans = -1;
    try
        n_trans = length(unique(results.condition_network.source_clusters));
    catch
    end
    fprintf(fid, 'Communication states (k-means on L-R sending/receiving PCs):  %d\n', prof.n_comm_clusters);
    if n_trans > 0
        fprintf(fid, 'Transcriptomic clusters:                                      %d\n', n_trans);
    end
else
    line('Communication profile data unavailable.');
end

if ~isempty(id) && isfield(id,'ari') && ~isnan(id.ari)
    blank();
    fprintf(fid, 'Adjusted Rand Index (comm_state vs trans_cluster):  %.4f\n', id.ari);
    fprintf(fid, 'Normalised Mutual Information:                       %.4f\n', id.nmi);
    fprintf(fid, 'Interpretation:                                      %s\n',  id.interpretation);
    fprintf(fid, 'Divergent cells (CDS > 0.5):                        %d / %d  (%.1f%%)\n', ...
        id.n_divergent, length(id.divergence), 100*id.divergent_fraction);
    blank();
    if strcmp(id.interpretation, 'Decoupled')
        fprintf(fid, 'KEY FINDING: Communication identity is ORTHOGONAL to transcriptomic\n');
        fprintf(fid, 'identity (ARI = %.4f). This means a cell''s communication behaviour\n', id.ari);
        fprintf(fid, 'CANNOT be predicted from its transcriptomic cluster membership alone.\n');
        fprintf(fid, 'This population of %.1f%% divergent cells is invisible to all\n', ...
            100*id.divergent_fraction);
        fprintf(fid, 'conventional sc-RNA-seq clustering approaches.\n');
    elseif strcmp(id.interpretation, 'Partial')
        fprintf(fid, 'Partial decoupling: some cell types communicate differently from\n');
        fprintf(fid, 'what their transcriptomic cluster predicts. Divergent cells (CDS>0.5)\n');
        fprintf(fid, 'are candidates for targeted follow-up.\n');
    else
        fprintf(fid, 'Communication identity largely mirrors transcriptomic identity in\n');
        fprintf(fid, 'this dataset. Consider whether cluster resolution is sufficient.\n');
    end
    blank();
    fprintf(fid, 'Top 10 most divergent cells:\n');
    fprintf(fid, '%-20s %15s %15s %10s\n', 'CellID','TransCluster','CommState','CDS');
    fprintf(fid, '%s\n', repmat('-',1,65));
    n_print = min(10, height(id.divergent_cells_table));
    for r = 1:n_print
        row = id.divergent_cells_table(r,:);
        cid_str = num2str(row.CellID);
        if iscell(row.CellID), cid_str = row.CellID{1}; end
        fprintf(fid, '%-20s %15s %15d %10.4f\n', cid_str, ...
            num2str(row.TranscriptomicCluster), row.CommunicationState, row.DivergenceScore);
    end
else
    blank();
    line('Communication identity analysis not available (requires cell_comm_profiles).');
end

%% Section 7 — Downstream response
hdr('7. DOWNSTREAM TRANSCRIPTIONAL RESPONSE');
dsr = results.downstream_response;
if isstruct(dsr) && ~isempty(fieldnames(dsr))
    fnames = fieldnames(dsr);
    n_downstream = sum(cellfun(@(f) isfield(dsr.(f),'source') && strcmp(dsr.(f).source,'downstream'), fnames));
    n_fallback   = sum(cellfun(@(f) isfield(dsr.(f),'source') && strcmp(dsr.(f).source,'receptor_fallback'), fnames));
    fprintf(fid, 'Receptors profiled:             %d\n', length(fnames));
    fprintf(fid, 'With downstream response genes: %d\n', n_downstream);
    fprintf(fid, 'Using receptor fallback:        %d\n', n_fallback);
    blank();
    fprintf(fid, 'Interpretation: For receptors with downstream response genes, channel\n');
    fprintf(fid, 'capacity is computed on the PC1 of the response gene matrix rather\n');
    fprintf(fid, 'than receptor expression. This measures actual transcriptional change\n');
    fprintf(fid, 'in receiver cells, not just receptor abundance.\n');
else
    line('Downstream response data unavailable or not run.');
end

%% Section 8 — Methods note
hdr('8. METHODS NOTE (for paper writing)');
if using_spearman
    fprintf(fid, 'Channel capacity: SPEARMAN COPULA FALLBACK (not KSG).\n');
    fprintf(fid, '  MI = -0.5*log(1-rho^2)  (Spearman rank correlation).\n');
    fprintf(fid, '  Used because n_donors < 10. Capacity values inflated.\n');
else
fprintf(fid, 'Channel capacity estimation: KSG k-NN mutual information estimator\n');
end
fprintf(fid, '  (Kraskov, Stoegbauer & Grassberger 2004, Phys Rev E 69:066138)\n');
fprintf(fid, '  I(X;Y) = psi(k) + psi(n) - <psi(nx+1) + psi(ny+1)>\n');
fprintf(fid, '  Capacity = I(ligand ; response) / H(ligand)  [0..1]\n');
fprintf(fid, '  Permutation test: 199 shuffles of Y, p = fraction(null >= observed)\n\n');
fprintf(fid, 'Bootstrap confidence intervals:\n');
fprintf(fid, '  500 resamples of the n paired donor observations with replacement.\n');
fprintf(fid, '  95%% CI = [2.5th, 97.5th] percentile of bootstrap mi_norm distribution.\n');
fprintf(fid, '  Adaptive k used consistently across all bootstrap resamples.\n');
fprintf(fid, '  Wide CIs at n=16 (~0.10-0.20 bits) reflect irreducible small-sample\n');
fprintf(fid, '  variance; relative ranking between channels remains valid.\n\n');
fprintf(fid, 'Differential testing: two-sample Wilcoxon rank-sum on per-cluster-pair\n');
fprintf(fid, '  strengths (strength test) and capacities (capacity test). BH-FDR.\n\n');
fprintf(fid, 'Causal validation: instrumental variable / 2SLS regression.\n');
fprintf(fid, '  Instrument: donor_id (natural experiment, analogous to MR in GWAS).\n');
fprintf(fid, '  Weak instrument filter: first-stage F > 10.\n\n');
fprintf(fid, 'NMF: non-negative matrix factorization of the [cluster-pair × L-R]\n');
fprintf(fid, '  communication tensor. Rank selected by elbow (<10%% improvement).\n\n');
fprintf(fid, 'Normalization: library-size median-ratio (scran-style), log1p.\n');
fprintf(fid, 'HVG selection: degree-2 polynomial regression log(var)~log(mean);\n');
fprintf(fid, '  top 3000 genes by residual variance.\n');

fclose(fid);
fprintf('Science report written: %s\n\n', report_path);

end

function generate_paper_figures(condition_network, condition_specific, ...
    identity_results, causal_results, comm_programs, cell_comm_profiles, ...
    landscape, config)
% GENERATE_PAPER_FIGURES  Publication-quality figures for the NEUROSCOPE paper.
%
% Each figure panel corresponds to a specific claim in the manuscript.
% Figures are saved as PDF (vector) and PNG (300 dpi) to results/paper_figures/.
%
% FIGURE PANELS GENERATED:
%
%   Fig 1A — Capacity distribution (histogram)
%            CLAIM: most co-expressed L-R pairs have capacity ≈ 0
%
%   Fig 1B — Capacity vs co-expression scatter
%            CLAIM: low Spearman r shows capacity is a different ranking
%
%   Fig 1C — 2×2 enrichment bar chart
%            CLAIM: capacity-only subset is missed by all existing tools
%
%   Fig 2A — Channel type composition (pie/bar)
%            CLAIM: HighCapacity_Broadcast channels are the novel finding
%
%   Fig 2B — Pathway distribution per channel type
%            CLAIM: capacity-only channels are pathway-distinct
%
%   Fig 3A — Communication states vs transcriptomic clusters (heatmap)
%            CLAIM: ARI < threshold means communication ≠ transcriptome
%
%   Fig 3B — Communication Divergence Score (CDS) distribution
%            CLAIM: a population of high-CDS cells exists
%
%   Fig 4A — Causal validation: IV F-statistic distribution
%            CLAIM: donor instrument is strong (F > 10 for validated pairs)
%
%   Fig 4B — Capacity in causally validated vs non-validated interactions
%            CLAIM: causally validated pairs have higher capacity
%
%   Fig 5A — NMF program loadings heatmap
%            CLAIM: co-regulated communication programs exist
%
%   Fig 5B — Program cluster-pair participation bubble chart
%            CLAIM: programs span multiple cell types (not cluster-specific)

fprintf('Generating publication figures...\n\n');

out_dir = fullfile(config.output_dir, 'paper_figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

cs  = condition_specific;
net = condition_network;

%% =========================================================================
%% Figure 1 — The I≈0 result  (Core Claim Figure)
%% =========================================================================

fig1 = figure('Name','Fig1_Information_Landscape', ...
    'Position',[100 100 1200 400], 'Color','w');

%--- 1A: Capacity histogram ---
ax1 = subplot(1,3,1, 'Parent', fig1);
if ~isempty(cs) && ismember('CondCapacity', cs.Properties.VariableNames)
    cap = cs.CondCapacity;
    fdr_cut = 0.05;
    sig_cap = cs.FDRCapacity < fdr_cut;

    % Histogram coloured by significance
    edges = linspace(0, max(max(cap)+0.01, 1), 30);
    n_nonsig = histcounts(cap(~sig_cap), edges);
    n_sig    = histcounts(cap(sig_cap),  edges);

    b = bar(ax1, edges(1:end-1), [n_nonsig(:), n_sig(:)], 'stacked', 'BarWidth', 1);
    b(1).FaceColor = [0.75 0.75 0.75];  % grey: I≈0
    b(2).FaceColor = [0.18 0.45 0.70];  % blue: significant
    b(1).EdgeColor = 'none';
    b(2).EdgeColor = 'none';

    xlabel(ax1, 'Channel capacity (bits)', 'FontSize', 11);
    ylabel(ax1, 'Number of L-R pairs',      'FontSize', 11);
    title(ax1,  'A   Capacity distribution', 'FontSize', 12, 'FontWeight', 'bold');

    zero_pct = 100 * mean(~sig_cap);
    text(ax1, max(cap)*0.5, max([n_nonsig n_sig])*0.85, ...
        sprintf('%.0f%% of pairs: I \\approx 0', zero_pct), ...
        'FontSize', 10, 'Color', [0.5 0.5 0.5]);

    legend(ax1, {'I≈0 (FDR≥0.05)', 'I>0 (FDR<0.05)'}, ...
        'Location','northeast', 'FontSize', 9);
    set(ax1, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
end

%--- 1B: Capacity vs strength scatter ---
ax2 = subplot(1,3,2, 'Parent', fig1);
if ~isempty(cs) && ismember('CondCapacity', cs.Properties.VariableNames) && ...
   ismember('CondStrength', cs.Properties.VariableNames)

    cap = cs.CondCapacity;
    str = cs.CondStrength;
    sig = cs.FDRCapacity < fdr_cut;

    scatter(ax2, str(~sig), cap(~sig), 15, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5);
    hold(ax2, 'on');
    scatter(ax2, str(sig),  cap(sig),  20, [0.18 0.45 0.70], 'filled', 'MarkerFaceAlpha', 0.8);

    % Rank correlation
    r = corr(tiedrank(cap), tiedrank(str), 'type', 'Pearson');
    text(ax2, min(str)+0.01*(max(str)-min(str)), max(cap)*0.9, ...
        sprintf('Spearman r = %.2f', r), 'FontSize', 10);

    xlabel(ax2, 'Co-expression strength',     'FontSize', 11);
    ylabel(ax2, 'Channel capacity (bits)',     'FontSize', 11);
    title(ax2,  'B   Capacity vs co-expression', 'FontSize', 12, 'FontWeight', 'bold');
    legend(ax2, {'I≈0', 'I>0'}, 'Location', 'northwest', 'FontSize', 9);
    set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
    hold(ax2, 'off');
end

%--- 1C: 2×2 bar chart ---
ax3 = subplot(1,3,3, 'Parent', fig1);
if ~isempty(cs) && ismember('FDRCapacity', cs.Properties.VariableNames)
    cap_sig = cs.FDRCapacity < fdr_cut;
    str_sig = cs.FDR         < fdr_cut;
    n_total = height(cs);

    groups = {'Both enriched', 'Capacity only [novel]', 'Strength only', 'Neither'};
    counts = [sum(cap_sig & str_sig), sum(cap_sig & ~str_sig), ...
              sum(~cap_sig & str_sig), sum(~cap_sig & ~str_sig)];
    pcts   = 100 * counts / n_total;

    colors = [0.30 0.60 0.30;   % green: both
              0.18 0.45 0.70;   % blue:  capacity only
              0.80 0.50 0.10;   % orange: strength only
              0.75 0.75 0.75];  % grey: neither

    bh = bar(ax3, 1:4, pcts, 'FaceColor', 'flat');
    for i = 1:4
        bh.CData(i,:) = colors(i,:);
    end
    bh.EdgeColor = 'none';

    set(ax3, 'XTickLabel', groups, 'XTick', 1:4, 'FontSize', 9, ...
        'Box', 'off', 'TickDir', 'out');
    xtickangle(ax3, 35);
    ylabel(ax3, '% of L-R pairs', 'FontSize', 11);
    title(ax3, 'C   Enrichment classification', 'FontSize', 12, 'FontWeight', 'bold');

    for i = 1:4
        if pcts(i) > 0.5
            text(ax3, i, pcts(i)+0.3, sprintf('%.1f%%', pcts(i)), ...
                'HorizontalAlignment','center', 'FontSize', 9);
        end
    end
end

savefig_both(fig1, out_dir, 'Fig1_information_landscape');

%% =========================================================================
%% Figure 2 — Channel type and pathway analysis
%% =========================================================================

fig2 = figure('Name','Fig2_Channel_Types', ...
    'Position',[100 100 1000 400], 'Color','w');

%--- 2A: Channel type bar ---
ax4 = subplot(1,2,1, 'Parent', fig2);
if ~isempty(cs) && ismember('ChannelType', cs.Properties.VariableNames)
    types  = {'HighCapacity_Broadcast','HighCapacity_Selective', ...
              'LowCapacity_Broadcast', 'LowCapacity_Selective'};
    labels = {'HC Broadcast','HC Selective','LC Broadcast','LC Selective'};
    counts = cellfun(@(t) sum(strcmp(cs.ChannelType, t)), types);
    colors4 = [0.18 0.45 0.70; 0.30 0.65 0.90; 0.80 0.50 0.10; 0.90 0.80 0.60];

    bh4 = bar(ax4, 1:4, counts, 'FaceColor', 'flat');
    for i = 1:4, bh4.CData(i,:) = colors4(i,:); end
    bh4.EdgeColor = 'none';
    set(ax4, 'XTickLabel', labels, 'XTick', 1:4, 'FontSize', 9, ...
        'Box','off','TickDir','out');
    xtickangle(ax4, 30);
    ylabel(ax4, 'Number of channels', 'FontSize', 11);
    title(ax4, 'A   Channel type composition', 'FontSize', 12, 'FontWeight', 'bold');
end

%--- 2B: Top pathways (capacity-only vs both) ---
ax5 = subplot(1,2,2, 'Parent', fig2);
if ~isempty(cs) && ismember('Pathway', cs.Properties.VariableNames) && ...
   ismember('FDRCapacity', cs.Properties.VariableNames)

    cap_only_mask = cs.FDRCapacity < fdr_cut & cs.FDR >= fdr_cut;
    both_mask_p   = cs.FDRCapacity < fdr_cut & cs.FDR < fdr_cut;

    pways = unique(cs.Pathway);
    n_cap = cellfun(@(p) sum(strcmp(cs.Pathway(cap_only_mask), p)), pways);
    n_bth = cellfun(@(p) sum(strcmp(cs.Pathway(both_mask_p),  p)), pways);

    % Keep pathways with at least 1 hit
    has_hit = n_cap > 0 | n_bth > 0;
    pways   = pways(has_hit);
    n_cap   = n_cap(has_hit);
    n_bth   = n_bth(has_hit);

    % Sort by capacity-only count
    [~, ord] = sort(n_cap, 'descend');
    top_n    = min(8, length(pways));
    pways    = pways(ord(1:top_n));
    n_cap    = n_cap(ord(1:top_n));
    n_bth    = n_bth(ord(1:top_n));

    bh5 = barh(ax5, 1:top_n, [n_cap(:), n_bth(:)]);
    bh5(1).FaceColor = [0.18 0.45 0.70];  % blue: capacity-only
    bh5(2).FaceColor = [0.30 0.60 0.30];  % green: both
    bh5(1).EdgeColor = 'none';
    bh5(2).EdgeColor = 'none';
    set(ax5, 'YTickLabel', pways, 'YTick', 1:top_n, 'FontSize', 9, ...
        'Box','off','TickDir','out');
    xlabel(ax5, 'Number of channels', 'FontSize', 11);
    title(ax5, 'B   Pathway distribution by type', 'FontSize', 12, 'FontWeight', 'bold');
    legend(ax5, {'Capacity-only [novel]','Both enriched'}, ...
        'Location','southeast', 'FontSize', 9);
end

savefig_both(fig2, out_dir, 'Fig2_channel_types');

%% =========================================================================
%% Figure 3 — Communication identity
%% =========================================================================

if isstruct(identity_results) && isfield(identity_results, 'ari') && ...
   ~isnan(identity_results.ari)

    fig3 = figure('Name','Fig3_Communication_Identity', ...
        'Position',[100 100 1000 400], 'Color','w');

    %--- 3A: Contingency heatmap ---
    ax6 = subplot(1,2,1, 'Parent', fig3);
    C = identity_results.contingency;
    if ~isempty(C) && numel(C) > 1
        % Normalise by column (transcriptomic cluster size)
        col_sums = sum(C, 1);
        C_norm   = C ./ max(col_sums, 1);

        imagesc(ax6, C_norm);
        colormap(ax6, [1 1 1; parula(64)]);
        colorbar(ax6);
        xlabel(ax6, 'Transcriptomic cluster', 'FontSize', 11);
        ylabel(ax6, 'Communication state',    'FontSize', 11);
        title(ax6, sprintf('A   ARI = %.3f (%s)', ...
            identity_results.ari, identity_results.interpretation), ...
            'FontSize', 12, 'FontWeight', 'bold');
        set(ax6, 'FontSize', 9);

        % Mark diagonal to show ideal alignment
        hold(ax6,'on');
        dim = min(size(C_norm));
        for d = 1:dim
            rectangle(ax6, 'Position',[d-0.5 d-0.5 1 1], ...
                'EdgeColor','r','LineWidth',0.5,'LineStyle',':');
        end
        hold(ax6,'off');
    end

    %--- 3B: CDS distribution ---
    ax7 = subplot(1,2,2, 'Parent', fig3);
    cds = identity_results.divergence;
    if ~isempty(cds)
        histogram(ax7, cds, 30, 'FaceColor',[0.18 0.45 0.70], 'EdgeColor','none');
        hold(ax7,'on');
        xline(ax7, 0.5, 'r--', 'LineWidth', 1.5);
        pct_high = 100 * identity_results.divergent_fraction;
        text(ax7, 0.52, max(histcounts(cds,30))*0.85, ...
            sprintf('%.1f%% with CDS > 0.5', pct_high), ...
            'FontSize', 10, 'Color', 'r');
        hold(ax7,'off');
        xlabel(ax7, 'Communication Divergence Score (CDS)', 'FontSize', 11);
        ylabel(ax7, 'Number of cells',                      'FontSize', 11);
        title(ax7, 'B   CDS distribution', 'FontSize', 12, 'FontWeight', 'bold');
        set(ax7, 'Box','off','TickDir','out','FontSize',10);
    end

    savefig_both(fig3, out_dir, 'Fig3_communication_identity');
end

%% =========================================================================
%% Figure 4 — Causal validation
%% =========================================================================

if isstruct(causal_results) && isfield(causal_results, 'table') && ...
   istable(causal_results.table) && height(causal_results.table) > 0

    fig4 = figure('Name','Fig4_Causal_Validation', ...
        'Position',[100 100 1000 400], 'Color','w');
    ct = causal_results.table;

    %--- 4A: F-statistic distribution ---
    ax8 = subplot(1,2,1, 'Parent', fig4);
    if ismember('FStatistic', ct.Properties.VariableNames)
        fstat = ct.FStatistic;
        fstat_valid = fstat(isfinite(fstat));
        histogram(ax8, fstat_valid, 20, 'FaceColor',[0.18 0.45 0.70], 'EdgeColor','none');
        hold(ax8,'on');
        xline(ax8, 10, 'r--', 'LineWidth', 1.5);
        text(ax8, 10.5, max(histcounts(fstat_valid,20))*0.9, ...
            'F=10 threshold', 'FontSize', 9, 'Color', 'r');
        hold(ax8,'off');
        xlabel(ax8, 'IV first-stage F-statistic', 'FontSize', 11);
        ylabel(ax8, 'Number of L-R pairs',        'FontSize', 11);
        title(ax8, 'A   Instrument strength', 'FontSize', 12, 'FontWeight', 'bold');
        set(ax8, 'Box','off','TickDir','out','FontSize',10);
    end

    %--- 4B: Capacity in causal vs non-causal ---
    ax9 = subplot(1,2,2, 'Parent', fig4);
    if ismember('Causal', ct.Properties.VariableNames) && ...
       ismember('CondCapacity', ct.Properties.VariableNames)
        causal_cap     = ct.CondCapacity(logical(ct.Causal));
        noncausal_cap  = ct.CondCapacity(~logical(ct.Causal));

        data  = {causal_cap; noncausal_cap};
        clrs  = [0.18 0.45 0.70; 0.75 0.75 0.75];
        lbls  = {'Causally validated','Not validated'};

        hold(ax9,'on');
        for g = 1:2
            d = data{g};
            d = d(isfinite(d));
            boxchart(ax9, g*ones(size(d)), d, 'BoxFaceColor', clrs(g,:), ...
                'MarkerStyle', '.', 'MarkerColor', clrs(g,:));
        end
        hold(ax9,'off');
        set(ax9, 'XTickLabel', lbls, 'XTick', 1:2, ...
            'Box','off','TickDir','out','FontSize',10);
        ylabel(ax9, 'Channel capacity (bits)', 'FontSize', 11);
        title(ax9, 'B   Capacity by causal status', 'FontSize', 12, 'FontWeight', 'bold');

        % Wilcoxon test annotation
        if numel(causal_cap)>3 && numel(noncausal_cap)>3
            p_wrs = ranksum(causal_cap, noncausal_cap);
            text(ax9, 1.5, max([causal_cap; noncausal_cap])*1.05, ...
                sprintf('p = %.3g', p_wrs), 'HorizontalAlignment','center','FontSize',10);
        end
    end

    savefig_both(fig4, out_dir, 'Fig4_causal_validation');
end

%% =========================================================================
%% Figure 5 — NMF communication programs
%% =========================================================================

if isstruct(comm_programs) && isfield(comm_programs, 'rank') && ...
   comm_programs.rank > 0 && ...
   isfield(comm_programs, 'W') && ~isempty(comm_programs.W)

    fig5 = figure('Name','Fig5_NMF_Programs', ...
        'Position',[100 100 1000 420], 'Color','w');

    %--- 5A: Program loadings heatmap (L-R × programs) ---
    ax10 = subplot(1,2,1, 'Parent', fig5);
    W = comm_programs.W;
    % Normalise each program to [0,1]
    W_norm = W ./ (max(W,[],1) + 1e-10);
    % Show top 20 L-R pairs by max loading
    [~, top_idx] = sort(max(W_norm,[],2), 'descend');
    top_idx = top_idx(1:min(20, length(top_idx)));
    imagesc(ax10, W_norm(top_idx,:)');
    colormap(ax10, flipud(gray));
    colorbar(ax10);
    ylabel(ax10, 'Communication program', 'FontSize', 11);
    xlabel(ax10, 'L-R pair (top 20)',     'FontSize', 11);
    title(ax10, 'A   Program loadings (NMF)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax10, 'YTick', 1:comm_programs.rank, 'FontSize', 9, 'TickDir','out');

    if isfield(comm_programs, 'program_table') && istable(comm_programs.program_table)
        pt_prog = comm_programs.program_table;
        if ismember('Ligand', pt_prog.Properties.VariableNames) && ...
           ismember('Receptor', pt_prog.Properties.VariableNames)
            lr_labels = cellfun(@(l,r) sprintf('%s→%s',l,r), ...
                pt_prog.Ligand(top_idx), pt_prog.Receptor(top_idx), ...
                'UniformOutput', false);
            set(ax10, 'XTickLabel', lr_labels, 'XTick', 1:length(top_idx));
            xtickangle(ax10, 45);
        end
    end

    %--- 5B: Program purity bar ---
    ax11 = subplot(1,2,2, 'Parent', fig5);
    if isfield(comm_programs, 'program_purity') && ...
       ~isempty(comm_programs.program_purity)
        pur = comm_programs.program_purity;
        bar(ax11, 1:length(pur), pur, 'FaceColor',[0.18 0.45 0.70],'EdgeColor','none');
        ylim(ax11,[0 1]);
        xlabel(ax11, 'Communication program',  'FontSize', 11);
        ylabel(ax11, 'Program purity',          'FontSize', 11);
        title(ax11, 'B   Program purity', 'FontSize', 12, 'FontWeight', 'bold');
        set(ax11, 'Box','off','TickDir','out','FontSize',10);
        yline(ax11, 0.7, 'r--');
        text(ax11, length(pur)*0.8, 0.73, 'Purity threshold', ...
            'FontSize', 9, 'Color', 'r');
    end

    savefig_both(fig5, out_dir, 'Fig5_nmf_programs');
end

fprintf('Paper figures saved to: %s\n\n', out_dir);

end


%% =========================================================================
%% Local helpers
%% =========================================================================

function savefig_both(fig, out_dir, name)
% Save as PDF (vector) and PNG (300 dpi).
pdf_path = fullfile(out_dir, [name '.pdf']);
png_path = fullfile(out_dir, [name '.png']);
try
    exportgraphics(fig, pdf_path, 'ContentType','vector');
    fprintf('  Saved: %s.pdf\n', name);
catch
    try
        print(fig, pdf_path, '-dpdf', '-r0');
    catch
    end
end
try
    exportgraphics(fig, png_path, 'Resolution', 300);
    fprintf('  Saved: %s.png\n', name);
catch
    try
        print(fig, png_path, '-dpng', '-r300');
    catch
    end
end
close(fig);
end

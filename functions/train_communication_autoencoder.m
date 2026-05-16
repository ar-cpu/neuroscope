function [autoencoder, encoder, latent, perf] = train_communication_autoencoder( ...
        expr_matrix, gene_names, condition_specific, cell_metadata)
% TRAIN_COMMUNICATION_AUTOENCODER  Learn a compact latent representation of
%   cell communication profiles via a sparse autoencoder.
%
% The autoencoder is trained on per-cell expression profiles, with genes
% restricted to those participating in significant communication channels
% (as identified in condition_specific). This focuses the representation on
% communication-relevant variation rather than the full transcriptome.
%
% ARCHITECTURE
% ------------
%   Input layer  :  n_comm_genes features (log1p-normalised expression)
%   Encoder      :  64 → 32 units, ReLU + L1 sparsity penalty
%   Latent space :  16-dimensional bottleneck
%   Decoder      :  32 → 64 → n_comm_genes units, ReLU / linear output
%
% If the Deep Learning Toolbox is not available, the function falls back to
% a PCA-based dimensionality reduction and returns a struct that mirrors the
% autoencoder output interface.
%
% USAGE
%   [ae, enc, latent, perf] = train_communication_autoencoder( ...
%       expr_matrix, gene_names, condition_specific, cell_metadata)
%
% INPUTS
%   expr_matrix        - (n_cells × n_genes) sparse or full numeric matrix
%                        of normalised expression values.
%   gene_names         - Table with variable 'gene_name' listing the gene
%                        corresponding to each column of expr_matrix.
%   condition_specific - Differential channel table from
%                        identify_differential_channels(). Used to select
%                        communication-relevant genes.
%   cell_metadata      - Table with at minimum a 'Cluster' column.
%
% OUTPUTS
%   autoencoder - Trained autoencoder object (or PCA fallback struct). Use
%                 autoencoder.encode(X) to obtain latent codes for new data.
%   encoder     - Encoder sub-network (or PCA projection matrix).
%   latent      - (n_cells × n_latent) matrix of latent codes for the
%                 training data. Can be passed to train_malignancy_classifier.
%   perf        - Struct with performance metrics:
%                   .reconstruction_mse   Mean squared reconstruction error
%                   .variance_explained   Fraction of variance captured
%                   .n_genes_used         Number of communication genes used
%                   .method               'autoencoder' or 'pca_fallback'

if nargin < 3 || isempty(condition_specific)
    condition_specific = table();
end
if nargin < 4 || isempty(cell_metadata)
    cell_metadata = table();
end

fprintf('  [Autoencoder] Selecting communication-relevant genes...\n');

%% Select communication genes

gene_list = gene_names.gene_name;
if istable(gene_names) && any(strcmp(gene_names.Properties.VariableNames, 'gene_name'))
    gene_list = upper(gene_names.gene_name);
end

comm_genes = {};
if ~isempty(condition_specific) && height(condition_specific) > 0
    if any(strcmp(condition_specific.Properties.VariableNames, 'Ligand'))
        comm_genes = unique([upper(condition_specific.Ligand); ...
                             upper(condition_specific.Receptor)]);
    end
end

gene_idx = [];
for g = 1:length(comm_genes)
    idx = find(strcmp(gene_list, comm_genes{g}), 1);
    if ~isempty(idx), gene_idx(end+1) = idx; end  %#ok<AGROW>
end

if isempty(gene_idx)
    % Fallback: use top 500 most variable genes
    if issparse(expr_matrix), X_full = full(expr_matrix);
    else, X_full = expr_matrix; end
    gene_var = var(X_full, 0, 1);
    [~, gene_idx] = maxk(gene_var, min(500, size(expr_matrix, 2)));
end

gene_idx = unique(gene_idx);
if issparse(expr_matrix)
    X = full(expr_matrix(:, gene_idx));
else
    X = expr_matrix(:, gene_idx);
end
X = log1p(X);   % log1p normalisation

n_genes  = size(X, 2);
n_latent = min(16, floor(n_genes / 4));

fprintf('  [Autoencoder] %d cells × %d communication genes → %d-D latent space\n', ...
    size(X,1), n_genes, n_latent);

%% Train autoencoder or fall back to PCA

use_dl = license('test', 'Neural_Network_Toolbox') || ...
         license('test', 'Deep_Learning_Toolbox');

perf = struct();
perf.n_genes_used = n_genes;

if use_dl
    try
        [autoencoder, encoder, latent, perf] = train_with_deep_learning(X, n_latent, perf);
        return;
    catch ME
        fprintf('  [Autoencoder] Deep Learning Toolbox error (%s); using PCA fallback.\n', ...
            ME.message);
    end
end

% PCA fallback
[coeff, score, ~, ~, explained] = pca(X, 'NumComponents', n_latent, ...
    'Economy', true);

latent = score(:, 1:min(n_latent, size(score,2)));

% Reconstruction error
X_recon = latent * coeff(:, 1:size(latent,2))';
mse = mean((X - X_recon).^2, 'all');

perf.reconstruction_mse = mse;
perf.variance_explained  = sum(explained(1:n_latent)) / 100;
perf.method              = 'pca_fallback';

% Package as a struct that mimics the autoencoder interface
autoencoder = struct();
autoencoder.coeff   = coeff;
autoencoder.n_latent = n_latent;
autoencoder.gene_idx = gene_idx;
autoencoder.encode  = @(Xnew) log1p(full(Xnew(:, gene_idx))) * coeff(:, 1:n_latent);

encoder = coeff(:, 1:n_latent);

fprintf('  [Autoencoder] PCA: %.1f%% variance in %d PCs, MSE=%.4f\n', ...
    perf.variance_explained * 100, n_latent, mse);

end


%% ── Deep Learning Toolbox implementation ─────────────────────────────────────

function [autoencoder, encoder, latent, perf] = train_with_deep_learning(X, n_latent, perf)

n_genes = size(X, 2);

% Encoder: input → 64 → 32 → n_latent
enc_layers = [
    featureInputLayer(n_genes, 'Normalization', 'zscore', 'Name', 'input')
    fullyConnectedLayer(64,     'Name', 'enc1')
    reluLayer(                  'Name', 'relu1')
    fullyConnectedLayer(32,     'Name', 'enc2')
    reluLayer(                  'Name', 'relu2')
    fullyConnectedLayer(n_latent,'Name','latent')
];

% Decoder: n_latent → 32 → 64 → output
dec_layers = [
    fullyConnectedLayer(32,     'Name', 'dec1')
    reluLayer(                  'Name', 'relu3')
    fullyConnectedLayer(64,     'Name', 'dec2')
    reluLayer(                  'Name', 'relu4')
    fullyConnectedLayer(n_genes,'Name', 'output')
    regressionLayer(            'Name', 'reglayer')
];

layers = [enc_layers; dec_layers];

opts = trainingOptions('adam', ...
    'MaxEpochs',         50,  ...
    'MiniBatchSize',     min(128, floor(size(X,1)/4)), ...
    'InitialLearnRate',  1e-3, ...
    'L2Regularization',  1e-4, ...
    'Verbose',           false, ...
    'Plots',             'none');

autoencoder = trainNetwork(X, X, layers, opts);
encoder     = autoencoder;  % full network; latent = activations at 'latent' layer

latent = activations(autoencoder, X, 'latent', 'OutputAs', 'rows');

X_recon = predict(autoencoder, X);
mse = mean((X - X_recon).^2, 'all');

perf.reconstruction_mse = mse;
perf.variance_explained  = 1 - mse / var(X(:));
perf.method              = 'autoencoder';

fprintf('  [Autoencoder] Deep Learning: MSE=%.4f, VarExpl=%.1f%%\n', ...
    mse, perf.variance_explained * 100);

end

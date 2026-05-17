function [autoencoder, encoder, latent, perf] = train_communication_autoencoder( ...
        expr_matrix, gene_names, condition_specific, cell_metadata)
% TRAIN_COMMUNICATION_AUTOENCODER  Sparse autoencoder on communication-gene profiles.
%
% Architecture: n_comm_genes → 64 → 32 → 16-D latent → 32 → 64 → n_comm_genes
% PCA fallback when toolbox unavailable.
%
% Usage
%   [ae, enc, latent, perf] = train_communication_autoencoder( ...
%       expr_matrix, gene_names, condition_specific, cell_metadata)
%
% Inputs
%   expr_matrix        - (n_cells × n_genes) sparse or full matrix of normalised expression
%   gene_names         - Table with variable 'gene_name' per column of expr_matrix
%   condition_specific - Differential channel table from identify_differential_channels()
%   cell_metadata      - Table with at minimum a 'Cluster' column
%
% Outputs
%   autoencoder - Trained network or PCA fallback struct; autoencoder.encode(X) gives codes
%   encoder     - Encoder sub-network or PCA projection matrix
%   latent      - (n_cells × n_latent) latent codes for training data
%   perf        - Struct: .reconstruction_mse, .variance_explained, .n_genes_used, .method

if nargin < 3 || isempty(condition_specific)
    condition_specific = table();
end
if nargin < 4 || isempty(cell_metadata)
    cell_metadata = table();
end

fprintf('  [Autoencoder] Selecting communication-relevant genes...\n');

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
X = log1p(X);

n_genes  = size(X, 2);
n_latent = min(16, floor(n_genes / 4));

fprintf('  [Autoencoder] %d cells × %d communication genes → %d-D latent space\n', ...
    size(X,1), n_genes, n_latent);

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

X_recon = latent * coeff(:, 1:size(latent,2))';
mse = mean((X - X_recon).^2, 'all');

perf.reconstruction_mse = mse;
perf.variance_explained  = sum(explained(1:n_latent)) / 100;
perf.method              = 'pca_fallback';

% Struct mirroring the autoencoder interface
autoencoder = struct();
autoencoder.coeff    = coeff;
autoencoder.n_latent = n_latent;
autoencoder.gene_idx = gene_idx;
autoencoder.encode   = @(Xnew) log1p(full(Xnew(:, gene_idx))) * coeff(:, 1:n_latent);

encoder = coeff(:, 1:n_latent);

fprintf('  [Autoencoder] PCA: %.1f%% variance in %d PCs, MSE=%.4f\n', ...
    perf.variance_explained * 100, n_latent, mse);

end


% Deep Learning Toolbox implementation

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

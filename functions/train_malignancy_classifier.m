function [classifier, perf] = train_malignancy_classifier( ...
        condition_latent, ref_latent, test_fraction)
% TRAIN_MALIGNANCY_CLASSIFIER  Train a binary classifier to distinguish
%   condition cells from reference cells in latent communication space.
%
% Uses the latent representations produced by train_communication_autoencoder
% to learn a compact separator between two biological states. This reveals
% which communication-latent dimensions best discriminate the conditions.
%
% When the Statistics and Machine Learning Toolbox is available, trains a
% regularised logistic regression (lasso penalty) or linear SVM. Falls back
% to a nearest-centroid (LDA-style) classifier otherwise.
%
% USAGE
%   [classifier, perf] = train_malignancy_classifier(condition_latent, ref_latent)
%   [classifier, perf] = train_malignancy_classifier(condition_latent, ref_latent, test_fraction)
%
% INPUTS
%   condition_latent - (n_cells_cond × n_latent) matrix of latent codes for
%                      the condition of interest (from train_communication_autoencoder).
%   ref_latent       - (n_cells_ref × n_latent) matrix of latent codes for
%                      the reference / control cells.
%   test_fraction    - Fraction of cells held out for evaluation (default: 0.2).
%
% OUTPUTS
%   classifier - Trained classifier object, or a struct with a .predict
%                function handle for the fallback nearest-centroid case.
%   perf       - Struct with performance metrics:
%                  .accuracy          Overall accuracy on held-out test set
%                  .auc               Area under ROC curve (condition = positive)
%                  .n_condition       Number of condition cells used
%                  .n_reference       Number of reference cells used
%                  .method            'svm', 'logistic', or 'nearest_centroid'

if nargin < 3 || isempty(test_fraction)
    test_fraction = 0.2;
end

fprintf('  [Classifier] %d condition vs %d reference cells in %d-D latent space\n', ...
    size(condition_latent,1), size(ref_latent,1), size(condition_latent,2));

%% Assemble labelled dataset

X = [condition_latent; ref_latent];
y = [ones(size(condition_latent,1),1); zeros(size(ref_latent,1),1)];

n = size(X, 1);
rng(42);
idx_perm  = randperm(n);
n_test    = max(1, round(n * test_fraction));
test_idx  = idx_perm(1:n_test);
train_idx = idx_perm(n_test+1:end);

X_train = X(train_idx, :);  y_train = y(train_idx);
X_test  = X(test_idx,  :);  y_test  = y(test_idx);

perf = struct();
perf.n_condition = size(condition_latent, 1);
perf.n_reference = size(ref_latent, 1);

%% Attempt toolbox classifiers

has_stats = license('test', 'Statistics_Toolbox') || ...
            license('test', 'Statistics_and_Machine_Learning_Toolbox');

if has_stats
    try
        [classifier, perf] = train_with_svm(X_train, y_train, X_test, y_test, perf);
        return;
    catch ME_svm
        try
            [classifier, perf] = train_with_logistic(X_train, y_train, X_test, y_test, perf);
            return;
        catch ME_lr
            fprintf('  [Classifier] Toolbox error (%s); using nearest-centroid fallback.\n', ...
                ME_lr.message);
        end
    end
end

%% Nearest-centroid fallback

[classifier, perf] = train_nearest_centroid(X_train, y_train, X_test, y_test, perf);

end


%% ── SVM ──────────────────────────────────────────────────────────────────────

function [classifier, perf] = train_with_svm(X_train, y_train, X_test, y_test, perf)

classifier = fitcsvm(X_train, y_train, ...
    'KernelFunction', 'linear', ...
    'Standardize',    true,     ...
    'BoxConstraint',  1.0);

[y_pred, scores] = predict(classifier, X_test);
perf.accuracy = mean(y_pred == y_test);

[~, ~, ~, perf.auc] = perfcurve(y_test, scores(:,2), 1);
perf.method = 'svm';

fprintf('  [Classifier] SVM: accuracy=%.3f, AUC=%.3f\n', perf.accuracy, perf.auc);
end


%% ── Logistic regression ──────────────────────────────────────────────────────

function [classifier, perf] = train_with_logistic(X_train, y_train, X_test, y_test, perf)

classifier = fitclinear(X_train, y_train, ...
    'Learner',        'logistic', ...
    'Regularization', 'lasso',   ...
    'Lambda',         1e-3);

[y_pred, scores] = predict(classifier, X_test);
perf.accuracy = mean(y_pred == y_test);

[~, ~, ~, perf.auc] = perfcurve(y_test, scores(:,2), 1);
perf.method = 'logistic';

fprintf('  [Classifier] Logistic: accuracy=%.3f, AUC=%.3f\n', perf.accuracy, perf.auc);
end


%% ── Nearest-centroid fallback ────────────────────────────────────────────────

function [classifier, perf] = train_nearest_centroid(X_train, y_train, X_test, y_test, perf)

centroid_pos = mean(X_train(y_train == 1, :), 1);
centroid_neg = mean(X_train(y_train == 0, :), 1);

predict_fn = @(X) double( ...
    sum((X - centroid_pos).^2, 2) < sum((X - centroid_neg).^2, 2));

y_pred = predict_fn(X_test);
perf.accuracy = mean(y_pred == y_test);

% Signed distance to decision boundary as soft score for AUC
dist_pos = sqrt(sum((X_test - centroid_pos).^2, 2));
dist_neg = sqrt(sum((X_test - centroid_neg).^2, 2));
score    = dist_neg - dist_pos;  % larger → closer to positive centroid

[~, ~, ~, perf.auc] = perfcurve(y_test, score, 1);
perf.method = 'nearest_centroid';

classifier = struct();
classifier.centroid_pos = centroid_pos;
classifier.centroid_neg = centroid_neg;
classifier.predict      = predict_fn;

fprintf('  [Classifier] Nearest-centroid: accuracy=%.3f, AUC=%.3f\n', ...
    perf.accuracy, perf.auc);
end

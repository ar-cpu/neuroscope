function [classifier, perf] = train_malignancy_classifier( ...
        condition_latent, ref_latent, test_fraction)
% TRAIN_MALIGNANCY_CLASSIFIER  Binary classifier on latent communication space.
%
% Usage
%   [classifier, perf] = train_malignancy_classifier(condition_latent, ref_latent)
%   [classifier, perf] = train_malignancy_classifier(condition_latent, ref_latent, test_fraction)
%
% Inputs
%   condition_latent - (n_cells_cond × n_latent) latent codes from train_communication_autoencoder
%   ref_latent       - (n_cells_ref × n_latent) latent codes for reference/control cells
%   test_fraction    - Fraction held out for evaluation (default: 0.2)
%
% Outputs
%   classifier - Trained model, or struct with .predict handle (nearest-centroid fallback)
%   perf       - Struct: .accuracy, .auc, .n_condition, .n_reference, .method

if nargin < 3 || isempty(test_fraction)
    test_fraction = 0.2;
end

fprintf('  [Classifier] %d condition vs %d reference cells in %d-D latent space\n', ...
    size(condition_latent,1), size(ref_latent,1), size(condition_latent,2));

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

[classifier, perf] = train_nearest_centroid(X_train, y_train, X_test, y_test, perf);

end


% SVM classifier

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


% Logistic regression classifier

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


% Nearest-centroid fallback

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

% Manual trapezoidal AUROC — no toolbox required
[~, ord]     = sort(score, 'descend');
labs_sorted  = y_test(ord);
n_pos        = max(sum(y_test == 1), 1);
n_neg        = max(sum(y_test == 0), 1);
tp           = cumsum(labs_sorted == 1) / n_pos;
fp           = cumsum(labs_sorted == 0) / n_neg;
perf.auc     = trapz([0; fp], [0; tp]);
perf.method  = 'nearest_centroid';

classifier = struct();
classifier.centroid_pos = centroid_pos;
classifier.centroid_neg = centroid_neg;
classifier.predict      = predict_fn;

fprintf('  [Classifier] Nearest-centroid: accuracy=%.3f, AUC=%.3f\n', ...
    perf.accuracy, perf.auc);
end

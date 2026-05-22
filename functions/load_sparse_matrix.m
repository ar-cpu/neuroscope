function sparse_mat = load_sparse_matrix(filename)
% LOAD_SPARSE_MATRIX  Load a sparse gene-by-cell (or cell-by-gene) matrix
%                     from either a native MATLAB .mat file or an HDF5-backed
%                     scipy sparse CSR archive produced by Python preprocessing
%                     pipelines (e.g., Scanpy, scVI).
%
% The function first attempts to load a standard MATLAB .mat file that contains
% a variable named 'sparse_matrix'. If that fails — most commonly because the
% file is in HDF5 format — it falls back to reading the three CSR arrays
% (/data, /indices, /indptr) and the /shape attribute directly from the HDF5
% store and reconstructs a native MATLAB sparse matrix.
%
% Index conversion: scipy stores CSR indices with 0-based offsets. These are
% converted to 1-based MATLAB indices before the sparse() call.
%
% Usage:
%   sparse_mat = load_sparse_matrix(filename)
%
% Inputs:
%   filename   - Path to a .mat or HDF5 file containing the sparse matrix.
%                For HDF5 files the following datasets must exist at the root:
%                  /data    - non-zero values (float32 or float64)
%                  /indices - column indices, 0-based (int32)
%                  /indptr  - row pointer array, 0-based (int32)
%                  /shape   - [n_rows, n_cols] (int32 or int64)
%
% Outputs:
%   sparse_mat - MATLAB sparse double matrix of size [n_rows x n_cols].
%
% Errors:
%   Raises an error if neither format can be parsed or if HDF5 datasets are
%   missing. The caller should wrap this in a try/catch for graceful fallback.

    try
        data = load(filename);
        fields = fieldnames(data);
        mat_field = '';
        for preferred = {'sparse_matrix','X','expr_matrix'}
            if isfield(data, preferred{1})
                mat_field = preferred{1};
                break;
            end
        end
        if ~isempty(mat_field)
            sparse_mat = data.(mat_field);
            fprintf('     Loaded standard format (%s field): %dx%d\n', mat_field, size(sparse_mat));
        elseif ~isempty(fields)
            sparse_mat = data.(fields{1});
            fprintf('     Loaded standard format (%s field): %dx%d\n', fields{1}, size(sparse_mat));
        else
            error('No sparse matrix field found in MAT file');
        end
    catch
        fprintf('     Loading HDF5 format...\n');
        data_vals = h5read(filename, '/data');
        indices = h5read(filename, '/indices');
        indptr = h5read(filename, '/indptr');
        shape_val = h5read(filename, '/shape');
        
        data_vals = double(data_vals);
        indices = double(indices) + 1;
        indptr = double(indptr) + 1;
        
        total_nnz = length(data_vals);
        row_ind = zeros(total_nnz, 1);
        col_ind = zeros(total_nnz, 1);
        
        pos = 1;
        for row_idx = 1:(length(indptr)-1)
            start_idx = indptr(row_idx);
            end_idx = indptr(row_idx+1) - 1;
            if start_idx <= end_idx
                n_elem = end_idx - start_idx + 1;
                row_ind(pos:(pos+n_elem-1)) = row_idx;
                col_ind(pos:(pos+n_elem-1)) = indices(start_idx:end_idx);
                pos = pos + n_elem;
            end
        end
        
        sparse_mat = sparse(row_ind, col_ind, data_vals, shape_val(1), shape_val(2));
        fprintf('     Loaded HDF5: %dx%d with %d non-zeros\n', ...
                size(sparse_mat,1), size(sparse_mat,2), nnz(sparse_mat));
    end
end
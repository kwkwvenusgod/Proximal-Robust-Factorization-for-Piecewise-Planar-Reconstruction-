K = 250;
root = 'C:\MyProject\Proximal-Robust-Factorization-for-Piecewise-Planar-Reconstruction-\Proximal-Robust-Factorization-for-Piecewise-Planar-Reconstruction\data\chessboardforward\';
files = dir([root '*.ppm']);
dispOn = true;

% infer the TSPs
[sp_labels] = TSP(K, root, files, dispOn);

% save the results
save('results/labels_chessboardforward_K250.mat', 'sp_labels');

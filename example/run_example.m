%% hpDPC example: A549 cells
% Reconstructs an hpDPC phase image from supplied pDPC and TIE phase maps
% and compares the result with the supplied reference hpDPC reconstruction.
%
% The supplied reference hpDPC TIFF is stored as 32-bit floating point.
% Fusion is performed in double precision, after which the recreated hpDPC
% image is converted to single precision for like-for-like comparison with
% the stored reference.
%
% Accompanying software for:
% "Hybrid Transport-of-Intensity and Polarization Differential Phase
% Contrast for Extended Spatial-Frequency Phase Imaging"
%
% Requirements:
%   MATLAB R2020a or later
%   hpDPC_fusion.m from the repository fusion directory
%
% Copyright (c) 2026 Fraser Montandon and contributors
% SPDX-License-Identifier: MIT

clearvars;
close all;
clc;

%% hpDPC fusion parameters

params.dx      = 3.45e-6 / 20;  % sample-plane sampling [m/pixel]
params.k0      = 0.05;           % normalized radial transition frequency
params.alpha   = 2;              % transition exponent
params.reg_eps = 1e-3;           % spectral regularization parameter

%% Locate repository files

scriptFolder = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(scriptFolder);
fusionFolder = fullfile(repositoryRoot, 'fusion');
dataFolder = fullfile(scriptFolder, 'A549');

addpath(fusionFolder);

assert(exist('hpDPC_fusion', 'file') == 2, ...
    ['hpDPC_fusion.m was not found. Expected location: ' ...
     fullfile(fusionFolder, 'hpDPC_fusion.m')]);

%% Load example phase maps

fprintf('Loading A549 example phase maps...\n');

% Use double precision for fusion inputs.
phiDPC = double(readPhase(fullfile(dataFolder, 'A549_pDPC')));
phiTIE = double(readPhase(fullfile(dataFolder, 'A549_TIE')));

% Preserve the native precision of the stored reference TIFF.
phiExpected = readPhase(fullfile(dataFolder, 'A549_hpDPC_reference'));

assert(isequal(size(phiDPC), size(phiTIE), size(phiExpected)), ...
    'pDPC, TIE, and reference hpDPC images must have identical dimensions.');

assert(isa(phiExpected, 'single'), ...
    ['A549_hpDPC_reference.tif must be stored as a 32-bit floating-point ' ...
     'TIFF (single precision).']);

fprintf('Image size: %d x %d pixels\n', size(phiTIE, 1), size(phiTIE, 2));

%% Recreate hpDPC reconstruction

fprintf('Running hpDPC fusion... ');
tFusion = tic;

phiRecreatedDouble = hpDPC_fusion( ...
    'phiTIE', phiTIE, ...
    'phiDPC', phiDPC, ...
    'dx', params.dx, ...
    'k0', params.k0, ...
    'alpha', params.alpha, ...
    'reg_eps', params.reg_eps);

fprintf('done (%.2f s).\n', toc(tFusion));

% Match the numerical representation of the stored reference.
phiRecreated = single(phiRecreatedDouble);

%% Compare with the supplied reference

errorDouble = double(phiRecreated) - double(phiExpected);
absoluteErrorDouble = abs(errorDouble);

rmse = sqrt(mean(errorDouble(:).^2));
mae = mean(absoluteErrorDouble(:));
maxAbsError = max(absoluteErrorDouble(:));

referenceRange = ...
    double(max(phiExpected(:))) - double(min(phiExpected(:)));

if referenceRange > 0
    relativeRMSE = rmse / referenceRange;
else
    relativeRMSE = NaN;
end

exactMatch = isequal(phiRecreated, phiExpected);

fprintf('\nA549 hpDPC reproduction\n');
fprintf('-----------------------\n');
fprintf('Exact match        : %s\n', char(string(exactMatch)));
fprintf('RMSE               : %.6e rad\n', rmse);
fprintf('MAE                : %.6e rad\n', mae);
fprintf('Maximum error      : %.6e rad\n', maxAbsError);
fprintf('Relative RMSE      : %.6e\n\n', relativeRMSE);

%% Display comparison

phaseImages = { ...
    phiDPC, ...
    phiTIE, ...
    double(phiRecreated), ...
    double(phiExpected)};

phaseMinimum = min(cellfun(@(x) min(x(:)), phaseImages));
phaseMaximum = max(cellfun(@(x) max(x(:)), phaseImages));
phaseLimits = [phaseMinimum, phaseMaximum];

if phaseLimits(1) == phaseLimits(2)
    phaseLimits = phaseLimits + [-0.5, 0.5];
end

figureHandle = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [100, 100, 1750, 430], ...
    'Name', 'A549 hpDPC example');

layout = tiledlayout(figureHandle, 1, 5, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

panelTitles = { ...
    'pDPC', ...
    'TIE', ...
    'Recreated hpDPC', ...
    'Reference hpDPC'};

phaseAxes = gobjects(1, 4);

for panelIndex = 1:4
    phaseAxes(panelIndex) = nexttile(layout);
    imagesc(phaseAxes(panelIndex), phaseImages{panelIndex}, phaseLimits);
    axis(phaseAxes(panelIndex), 'image');
    axis(phaseAxes(panelIndex), 'off');
    colormap(phaseAxes(panelIndex), 'gray');
    title(phaseAxes(panelIndex), panelTitles{panelIndex}, ...
        'FontWeight', 'normal');
end

phaseColorbar = colorbar(phaseAxes(4), 'eastoutside');
phaseColorbar.Label.String = 'Phase (rad)';
phaseColorbar.Label.FontWeight = 'normal';

axError = nexttile(layout);
imagesc(axError, absoluteErrorDouble);
axis(axError, 'image');
axis(axError, 'off');
colormap(axError, 'parula');
title(axError, '|Recreated - reference|', 'FontWeight', 'normal');

errorColorbar = colorbar(axError, 'eastoutside');
errorColorbar.Label.String = 'Absolute phase error (rad)';
errorColorbar.Label.FontWeight = 'normal';

if exactMatch
    figureTitle = 'A549 hpDPC reproduction | Exact match';
else
    figureTitle = sprintf( ...
        ['A549 hpDPC reproduction | RMSE = %.3e rad | ' ...
         'max error = %.3e rad'], ...
        rmse, maxAbsError);
end

title(layout, figureTitle, 'FontWeight', 'normal');

%% Local function

function phi = readPhase(stem)
%READPHASE Read a floating-point phase TIFF or MAT file.
%
% Accepted files:
%   .tif
%   .tiff
%   .mat
%
% A MAT file must contain a variable named "phi".

    candidates = { ...
        [stem '.tif'], ...
        [stem '.tiff'], ...
        [stem '.mat']};

    present = cellfun(@isfile, candidates);

    assert(nnz(present) == 1, ...
        'Expected exactly one .tif, .tiff, or .mat file for:\n%s', stem);

    filename = candidates{find(present, 1)};
    [~, ~, extension] = fileparts(filename);

    if strcmpi(extension, '.mat')
        data = load(filename, 'phi');
        assert(isfield(data, 'phi'), ...
            '%s must contain a variable named "phi".', filename);
        raw = data.phi;
    else
        info = imfinfo(filename);
        assert(numel(info) == 1, ...
            '%s must contain a single phase image.', filename);

        raw = imread(filename);

        assert(isfloat(raw), ...
            ['%s is stored as an integer TIFF. The example requires a ' ...
             'floating-point quantitative phase image.'], ...
            filename);
    end

    assert(isnumeric(raw) && isreal(raw) && ismatrix(raw) && ...
        ~isempty(raw) && all(isfinite(raw(:))), ...
        '%s must contain a finite, real, two-dimensional phase map.', ...
        filename);

    phi = raw;
end

%% hpDPC simulation
% Simulation accompanying:
% "Hybrid Transport-of-Intensity and Polarization Differential Phase
% Contrast for Extended Spatial-Frequency Phase Imaging"
%
% The script compares ground-truth phase with TIE, matched and unmatched
% pDPC, conventional hybrid TIE-pDPC, and hpDPC reconstructions.
%
% The partially coherent forward model uses discrete source-point sampling.
% targetSourceSamples = 80 is used for the manuscript configuration.
%
% Requirements:
%   MATLAB R2024b
%
% Copyright (c) 2026 Fraser Montandon and contributors
% SPDX-License-Identifier: MIT

clearvars;
close all;
clc;

%% Simulation parameters

lambdaUm = 0.523;              % Wavelength [um]
objectiveNA = 0.40;            % Objective numerical aperture
tieIlluminationNA = 0.15;      % Symmetric TIE illumination NA

etaMatched = 1.00;
etaUnmatched = 0.75;
dpcMatchedNA = etaMatched * objectiveNA;
dpcUnmatchedNA = etaUnmatched * objectiveNA;
annularInnerNA = tieIlluminationNA;

pixelSizeUm = 0.1725;          % Simulation-plane sampling [um]
nPixels = 512;                 % Square simulation grid size
defocusUm = 4.0;               % TIE defocus distance [um]

dpcRegularization = [1e-2, 2e-3];  % Absorption and phase regularization
tieIntensityOffset = 1e-2;          % Offset for low-intensity regions in TIE inversion

fusionPixelSize = 1;           % Only normalized spatial frequency is used
fusionCutoff = 0.05;
fusionOrder = 2;
fusionRegularization = 1e-3;

targetSourceSamples = 80;      % Source-point sampling target
runSourceConvergence = false;  % Optional convergence analysis

validateattributes(lambdaUm, {'numeric'}, {'real','finite','positive','scalar'});
validateattributes(objectiveNA, {'numeric'}, {'real','finite','positive','scalar','<=',1});
validateattributes(tieIlluminationNA, {'numeric'}, ...
    {'real','finite','nonnegative','scalar','<=',objectiveNA});
validateattributes(pixelSizeUm, {'numeric'}, {'real','finite','positive','scalar'});
validateattributes(nPixels, {'numeric'}, {'real','finite','integer','positive','scalar'});
validateattributes(defocusUm, {'numeric'}, {'real','finite','nonzero','scalar'});
validateattributes(targetSourceSamples, {'numeric'}, ...
    {'real','finite','integer','positive','scalar'});

if annularInnerNA >= dpcUnmatchedNA
    error('hpDPC:InvalidAnnulus', ...
        'The annular inner NA must be smaller than its outer NA.');
end

fprintf('\n============================================================\n');
fprintf(' hpDPC SIMULATION\n');
fprintf('============================================================\n');
fprintf('Objective NA             : %.2f\n', objectiveNA);
fprintf('TIE illumination NA      : %.2f\n', tieIlluminationNA);
fprintf('Matched pDPC NAill       : %.2f  (eta = %.2f)\n', ...
    dpcMatchedNA, etaMatched);
fprintf('Unmatched pDPC NAill     : %.2f  (eta = %.2f)\n', ...
    dpcUnmatchedNA, etaUnmatched);
fprintf('hpDPC central-hole NA    : %.2f\n', annularInnerNA);
fprintf('Target source samples    : %d\n', targetSourceSamples);
fprintf('============================================================\n\n');

%% Spatial-frequency grid and objective pupil

fieldWidthUm = nPixels * pixelSizeUm;
frequencySpacing = 1 / fieldWidthUm;
frequencyAxis = (-nPixels/2:nPixels/2-1) * frequencySpacing;
[frequencyX, frequencyY] = meshgrid(frequencyAxis, frequencyAxis);

objectivePupil = double( ...
    frequencyX.^2 + frequencyY.^2 <= (objectiveNA / lambdaUm)^2);

%% Ground-truth phase object

phaseTruth = createGroundTruthPhase(nPixels, fieldWidthUm);
objectSpectrum = fftshift(fft2(exp(1i * phaseTruth)));

%% Illumination distributions

dpcAnglesDeg = [0, 180, 90, 270];

sourceMatched = makeDPCSources( ...
    dpcAnglesDeg, dpcMatchedNA, 0, lambdaUm, frequencyX, frequencyY);

sourceUnmatched = makeDPCSources( ...
    dpcAnglesDeg, dpcUnmatchedNA, 0, lambdaUm, frequencyX, frequencyY);

annularInnerRatio = annularInnerNA / dpcUnmatchedNA;
sourceAnnular = makeDPCSources( ...
    dpcAnglesDeg, dpcUnmatchedNA, annularInnerRatio, ...
    lambdaUm, frequencyX, frequencyY);

radialNA = hypot(frequencyX, frequencyY) * lambdaUm;
sourceTIE = double(radialNA <= tieIlluminationNA);

%% Defocus propagation transfer function

defocusTransfer = exp(1i * 2*pi * defocusUm .* sqrt(max( ...
    0, (1/lambdaUm)^2 - frequencyX.^2 - frequencyY.^2)));

%% pDPC simulations and reconstructions

fprintf('Reconstructing matched pDPC...\n');
phaseDPCMatched = reconstructDPC( ...
    objectSpectrum, sourceMatched, objectivePupil, defocusTransfer, ...
    frequencyAxis, frequencySpacing, dpcRegularization, ...
    targetSourceSamples);

fprintf('Reconstructing unmatched pDPC...\n');
phaseDPCUnmatched = reconstructDPC( ...
    objectSpectrum, sourceUnmatched, objectivePupil, defocusTransfer, ...
    frequencyAxis, frequencySpacing, dpcRegularization, ...
    targetSourceSamples);

fprintf('Reconstructing annular pDPC...\n');
phaseDPCAnnular = reconstructDPC( ...
    objectSpectrum, sourceAnnular, objectivePupil, defocusTransfer, ...
    frequencyAxis, frequencySpacing, dpcRegularization, ...
    targetSourceSamples);

%% TIE simulation and reconstruction

fprintf('Reconstructing TIE...\n');
[~, intensityFocus, intensityDefocus] = runForwardModel( ...
    objectSpectrum, sourceTIE, objectivePupil, defocusTransfer, ...
    frequencyAxis, frequencySpacing, targetSourceSamples);

phaseTIE = reconstructTIE( ...
    intensityFocus, intensityDefocus, defocusUm, lambdaUm, ...
    frequencyX, frequencyY, tieIntensityOffset);

%% Spectral fusion

fprintf('Performing spectral fusion...\n');

% Hybrid reference: TIE + unmatched conventional pDPC.
phaseHybrid = fuseTIEDPC( ...
    'PhaseTIE', phaseTIE, ...
    'PhaseDPC', phaseDPCUnmatched, ...
    'PixelSize', fusionPixelSize, ...
    'Cutoff', fusionCutoff, ...
    'Order', fusionOrder, ...
    'Regularization', fusionRegularization);

% hpDPC: TIE + unmatched annular pDPC.
phaseHpDPC = fuseTIEDPC( ...
    'PhaseTIE', phaseTIE, ...
    'PhaseDPC', phaseDPCAnnular, ...
    'PixelSize', fusionPixelSize, ...
    'Cutoff', fusionCutoff, ...
    'Order', fusionOrder, ...
    'Regularization', fusionRegularization);

%% Comparison figure
% Shift each displayed phase map to zero minimum for visualization only.

zeroMinimum = @(array) array - min(array(:));
displayData = {
    zeroMinimum(phaseTruth)
    zeroMinimum(phaseTIE)
    zeroMinimum(phaseDPCMatched)
    zeroMinimum(phaseDPCUnmatched)
    zeroMinimum(phaseHybrid)
    zeroMinimum(phaseHpDPC)
};

panelTitles = {
    '(a) Ground-truth phase'
    '(b) TIE (\Deltaz = 4 \mum)'
    '(c) pDPC, \eta = 1'
    '(d) pDPC, \eta = 0.75'
    '(e) Hybrid TIE-pDPC'
    '(f) hpDPC'
};

plotComparisonFigure(displayData, panelTitles);

%% Optional source-sampling convergence test

if runSourceConvergence
    sourceSampleTargets = [40, 80, 160, 320];
    convergenceTable = runConvergenceTest( ...
        sourceSampleTargets, phaseHpDPC, objectSpectrum, sourceAnnular, ...
        sourceTIE, objectivePupil, defocusTransfer, frequencyAxis, ...
        frequencySpacing, dpcRegularization, defocusUm, lambdaUm, ...
        frequencyX, frequencyY, tieIntensityOffset, fusionPixelSize, ...
        fusionCutoff, fusionOrder, fusionRegularization);

    fprintf('\n');
    disp(convergenceTable);
end

fprintf('\nSimulation complete.\n');

%% Local functions

function phase = createGroundTruthPhase(nPixels, fieldWidthUm)
% Create the text-plus-microlens phase target.

    textFigure = figure( ...
        'Visible', 'off', ...
        'Position', [0, 0, nPixels, nPixels], ...
        'Color', 'k');
    figureCleanup = onCleanup(@() closeIfValid(textFigure)); 

    text(0.5, 0.5, 'UCT', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 140, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'Color', 'w');
    axis off;
    set(gca, 'Position', [0, 0, 1, 1]);

    drawnow;
    frame = getframe(gca);
    textImage = frame.cdata;

    % The rendered target is white on black, so one RGB channel is
    % sufficient for conversion to a two-dimensional intensity image.
    if ndims(textImage) == 3
        textImage = textImage(:, :, 1);
    end

    if ~isequal(size(textImage), [nPixels, nPixels])
        textImage = resizeImage(textImage, nPixels, nPixels);
    end

    textMask = double(textImage > 100);

    spatialAxis = linspace(-fieldWidthUm/2, fieldWidthUm/2, nPixels);
    [xGrid, yGrid] = meshgrid(spatialAxis, spatialAxis);
    lensRadius = fieldWidthUm / 2.2;
    radialPosition = hypot(xGrid, yGrid);

    lensPhase = 3.5 * (1 - (radialPosition/lensRadius).^2);
    lensPhase(radialPosition > lensRadius) = 0;

    phase = lensPhase + textMask;
    phase(phase < 0) = 0;
end

function closeIfValid(figureHandle)


    if isgraphics(figureHandle)
        close(figureHandle);
    end
end

function resized = resizeImage(image, nRows, nColumns)


    image = double(image);
    sourceRows = linspace(1, size(image, 1), nRows);
    sourceColumns = linspace(1, size(image, 2), nColumns);
    [queryColumns, queryRows] = meshgrid(sourceColumns, sourceRows);

    resized = interp2(image, queryColumns, queryRows, 'linear');
end

function sources = makeDPCSources( ...
        anglesDeg, outerNA, innerRatio, lambdaUm, frequencyX, frequencyY)
% Generate asymmetric DPC source distributions.

    sources = zeros([size(frequencyX), numel(anglesDeg)]);
    for sourceIndex = 1:numel(anglesDeg)
        sources(:, :, sourceIndex) = generateDPCSource( ...
            anglesDeg(sourceIndex), outerNA, innerRatio, lambdaUm, ...
            frequencyX, frequencyY);
    end
end

function source = generateDPCSource( ...
        angleDeg, outerNA, innerRatio, lambdaUm, frequencyX, frequencyY)
% Generate one semicircular or semiannular source.

    radialNA = hypot(frequencyX, frequencyY) * lambdaUm;
    aperture = radialNA <= outerNA & radialNA >= innerRatio * outerNA;

    switch angleDeg
        case 0
            halfPlane = frequencyX >= 0;
        case 90
            halfPlane = frequencyY >= 0;
        case 180
            halfPlane = frequencyX <= 0;
        case 270
            halfPlane = frequencyY <= 0;
        otherwise
            error('hpDPC:UnsupportedSourceAngle', ...
                'Only source angles of 0, 90, 180, and 270 degrees are supported.');
    end

    source = double(aperture & halfPlane);
end

function phase = reconstructDPC( ...
        objectSpectrum, source, pupil, defocusTransfer, frequencyAxis, ...
        frequencySpacing, regularization, targetSourceSamples)
% Simulate DPC measurements and reconstruct phase.

    intensity = runForwardModel( ...
        objectSpectrum, source, pupil, defocusTransfer, frequencyAxis, ...
        frequencySpacing, targetSourceSamples);
    normalizedData = prepareDPCInput(intensity);
    [~, phase] = solveDPC(normalizedData, source, pupil, regularization);
end

function [measuredIntensity, symmetricIntensity, defocusedIntensity] = ...
        runForwardModel(objectSpectrum, source, pupil, defocusTransfer, ...
        frequencyAxis, frequencySpacing, targetSourceSamples)
% Partially coherent, source-sampled image formation.
%
% Source-point subsampling uses:
%   stride = max(1, floor(numberOfSourcePixels/targetSourceSamples)).

    [nRows, nColumns, nSources] = size(source);
    measuredIntensity = zeros(nRows, nColumns, nSources);
    symmetricIntensity = zeros(nRows, nColumns);
    defocusedIntensity = zeros(nRows, nColumns);

    for sourceIndex = 1:nSources
        [sourceRows, sourceColumns] = find(source(:, :, sourceIndex));
        currentIntensity = zeros(nRows, nColumns);

        sourceStride = max(1, floor( ...
            numel(sourceColumns) / targetSourceSamples));

        for pointIndex = 1:sourceStride:numel(sourceColumns)
            sourceFrequencyX = frequencyAxis(sourceColumns(pointIndex));
            sourceFrequencyY = frequencyAxis(sourceRows(pointIndex));
            shiftX = round(sourceFrequencyX / frequencySpacing);
            shiftY = round(sourceFrequencyY / frequencySpacing);

            shiftedSpectrum = circshift( ...
                objectSpectrum, [-shiftY, -shiftX]) .* pupil;

            focusedField = ifft2(ifftshift(shiftedSpectrum));
            currentIntensity = currentIntensity + abs(focusedField).^2;

            defocusedField = ifft2(ifftshift( ...
                shiftedSpectrum .* defocusTransfer));
            defocusedIntensity = defocusedIntensity + abs(defocusedField).^2;
        end

        measuredIntensity(:, :, sourceIndex) = currentIntensity;
        symmetricIntensity = symmetricIntensity + currentIntensity;
    end

    backgroundIntensity = mean(symmetricIntensity(:));
    if ~isfinite(backgroundIntensity) || backgroundIntensity <= 0
        error('hpDPC:InvalidIntensityNormalization', ...
            'The simulated mean background intensity must be positive and finite.');
    end

    if nSources > 1
        measuredIntensity = measuredIntensity / (backgroundIntensity/nSources);
    else
        measuredIntensity = measuredIntensity / backgroundIntensity;
    end
    symmetricIntensity = symmetricIntensity / backgroundIntensity;
    defocusedIntensity = defocusedIntensity / backgroundIntensity;
end

function normalizedData = prepareDPCInput(intensity)
% Normalize DPC images and transform them to Fourier space.

    normalizedData = zeros(size(intensity));
    for imageIndex = 1:size(intensity, 3)
        image = intensity(:, :, imageIndex);
        normalizedData(:, :, imageIndex) = fft2( ...
            image / mean(image(:)) - 1);
    end
end

function [absorption, phase] = solveDPC( ...
        normalizedData, source, pupil, regularization)
% Perform Tikhonov-regularized pDPC inversion.

    validateattributes(regularization, {'numeric'}, ...
        {'real','finite','nonnegative','numel',2});

    phaseTransfer = zeros(size(source));
    absorptionTransfer = zeros(size(source));
    for sourceIndex = 1:size(source, 3)
        [phaseTransfer(:, :, sourceIndex), ...
            absorptionTransfer(:, :, sourceIndex)] = calculateDPCTransfer( ...
            source(:, :, sourceIndex), pupil);
    end

    matrix11 = sum(abs(absorptionTransfer).^2, 3) + regularization(1);
    matrix12 = sum(conj(absorptionTransfer) .* phaseTransfer, 3);
    matrix21 = sum(conj(phaseTransfer) .* absorptionTransfer, 3);
    matrix22 = sum(abs(phaseTransfer).^2, 3) + regularization(2);
    determinant = matrix11 .* matrix22 - matrix12 .* matrix21;

    data1 = sum(normalizedData .* conj(absorptionTransfer), 3);
    data2 = sum(normalizedData .* conj(phaseTransfer), 3);

    absorption = real(ifft2( ...
        (data1 .* matrix22 - data2 .* matrix12) ./ determinant));
    phase = real(ifft2( ...
        (data2 .* matrix11 - data1 .* matrix21) ./ determinant));
end

function [phaseTransfer, absorptionTransfer] = ...
        calculateDPCTransfer(source, pupil)
% Calculate phase and absorption transfer functions.

    shiftedSource = fftshift(source);
    paddedSource = zeros(size(shiftedSource) + 1, 'like', shiftedSource);
    paddedSource(1:end-1, 1:end-1) = shiftedSource;

    rotatedSource = rot90(paddedSource, 2);
    rotatedSource = ifftshift(rotatedSource(1:end-1, 1:end-1));

    crossSpectrum = conj(fft2(rotatedSource .* pupil)) .* fft2(pupil);
    absorptionTransfer = 2 * ifft2(real(crossSpectrum));
    phaseTransfer = 2 * ifft2(1i * imag(crossSpectrum));

    dcNormalization = sum(rotatedSource .* abs(pupil).^2, 'all');
    if dcNormalization == 0
        error('hpDPC:ZeroTransferNormalization', ...
            'The DPC transfer-function normalization is zero.');
    end

    absorptionTransfer = absorptionTransfer / dcNormalization;
    phaseTransfer = 1i * phaseTransfer / dcNormalization;
end

function phase = reconstructTIE( ...
        intensityFocus, intensityDefocus, defocusUm, lambdaUm, ...
        frequencyX, frequencyY, intensityOffset)
% Reconstruct phase using a two-plane forward difference.

    waveNumber = 2*pi / lambdaUm;
    axialDerivative = (intensityDefocus - intensityFocus) / defocusUm;

    frequencyXUnshifted = ifftshift(frequencyX);
    frequencyYUnshifted = ifftshift(frequencyY);
    laplacian = -4*pi^2 * ( ...
        frequencyXUnshifted.^2 + frequencyYUnshifted.^2);

    laplacian(1, 1) = 1;
    inverseLaplacian = 1 ./ laplacian;
    inverseLaplacian(1, 1) = 0;

    auxiliaryPotential = ifft2( ...
        fft2(-waveNumber * axialDerivative) .* inverseLaplacian);
    derivativeX = real(ifft2(fft2(auxiliaryPotential) .* ...
        (1i * 2*pi * frequencyXUnshifted)));
    derivativeY = real(ifft2(fft2(auxiliaryPotential) .* ...
        (1i * 2*pi * frequencyYUnshifted)));

    rightHandSide = ...
        fft2(derivativeX ./ (intensityFocus + intensityOffset)) .* ...
            (1i * 2*pi * frequencyXUnshifted) + ...
        fft2(derivativeY ./ (intensityFocus + intensityOffset)) .* ...
            (1i * 2*pi * frequencyYUnshifted);

    phase = real(ifft2(rightHandSide .* inverseLaplacian));
    phase = phase - mean(phase(:));
end

function [fusedPhase, diagnostics] = fuseTIEDPC(varargin)
% Fuse low-frequency TIE and higher-frequency pDPC phase.

    parser = inputParser;
    parser.FunctionName = mfilename;
    addParameter(parser, 'PhaseTIE', [], ...
        @(x) isnumeric(x) && ismatrix(x) && all(isfinite(x(:))));
    addParameter(parser, 'PhaseDPC', [], ...
        @(x) isnumeric(x) && ismatrix(x) && all(isfinite(x(:))));
    addParameter(parser, 'PixelSize', 1, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'Cutoff', 0.05, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'Order', 2, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'Regularization', 1e-3, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    parse(parser, varargin{:});

    phaseTIE = double(parser.Results.PhaseTIE);
    phaseDPC = double(parser.Results.PhaseDPC);
    pixelSize = parser.Results.PixelSize;
    cutoff = parser.Results.Cutoff;
    filterOrder = parser.Results.Order;
    regularization = parser.Results.Regularization;

    if isempty(phaseTIE) || isempty(phaseDPC)
        error('hpDPC:MissingFusionInput', ...
            'Both PhaseTIE and PhaseDPC must be supplied.');
    end
    if ~isequal(size(phaseTIE), size(phaseDPC))
        error('hpDPC:FusionSizeMismatch', ...
            'PhaseTIE and PhaseDPC must have the same dimensions.');
    end

    [nRows, nColumns] = size(phaseTIE);
    rowFrequency = ifftshift( ...
        (-floor(nRows/2):ceil(nRows/2)-1) / (nRows * pixelSize));
    columnFrequency = ifftshift( ...
        (-floor(nColumns/2):ceil(nColumns/2)-1) / (nColumns * pixelSize));
    [frequencyColumns, frequencyRows] = meshgrid( ...
        columnFrequency, rowFrequency);

    radialFrequency = hypot(frequencyColumns, frequencyRows);
    normalizedFrequency = radialFrequency / (max(radialFrequency(:)) + eps);

    tieWeight = exp(-(normalizedFrequency / (cutoff + eps)).^filterOrder);
    tieWeight = tieWeight / (max(tieWeight(:)) + eps);
    dpcWeight = 1 - tieWeight;
    dpcWeight(normalizedFrequency < 1e-6) = 0;

    fusedSpectrum = tieWeight .* fft2(phaseTIE) + ...
        dpcWeight .* fft2(phaseDPC);
    fusedSpectrum = fusedSpectrum ./ (1 + regularization .* ...
        (normalizedFrequency / 0.99).^2);

    fusedPhase = real(ifft2(fusedSpectrum));
    fusedPhase = fusedPhase - mean(fusedPhase(:)) + mean(phaseTIE(:));

    if nargout > 1
        diagnostics = struct( ...
            'RadialFrequency', radialFrequency, ...
            'NormalizedFrequency', normalizedFrequency, ...
            'TIEWeight', tieWeight, ...
            'DPCWeight', dpcWeight, ...
            'Cutoff', cutoff, ...
            'Order', filterOrder, ...
            'Regularization', regularization, ...
            'PixelSize', pixelSize);
    end
end

function figureHandle = plotComparisonFigure(displayData, panelTitles)
% Plot six comparison panels.

    allData = cat(3, displayData{:});
    commonLimits = [min(allData(:)), max(allData(:))];

    figureHandle = figure( ...
        'Color', 'w', ...
        'Units', 'centimeters', ...
        'Position', [2, 2, 17.8, 11.2]);
    layout = tiledlayout(figureHandle, 2, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for panelIndex = 1:numel(displayData)
        axesHandle = nexttile(layout);
        imagesc(axesHandle, displayData{panelIndex});
        axis(axesHandle, 'image');
        axis(axesHandle, 'off');
        clim(axesHandle, commonLimits);
        title(axesHandle, panelTitles{panelIndex}, ...
            'Interpreter', 'tex', ...
            'FontName', 'Arial', ...
            'FontSize', 9, ...
            'FontWeight', 'bold');
        set(axesHandle, 'FontName', 'Arial', 'LineWidth', 0.8);
    end

    colormap(figureHandle, 'parula');
    colorbarHandle = colorbar;
    colorbarHandle.Layout.Tile = 'east';
    ylabel(colorbarHandle, 'Phase (rad)', ...
        'FontName', 'Arial', 'FontSize', 10);
    colorbarHandle.FontName = 'Arial';
    colorbarHandle.FontSize = 9;
    colorbarHandle.LineWidth = 0.8;
end

function resultsTable = runConvergenceTest( ...
        sampleTargets, manuscriptPhase, objectSpectrum, sourceAnnular, ...
        sourceTIE, pupil, defocusTransfer, frequencyAxis, frequencySpacing, ...
        dpcRegularization, defocusUm, lambdaUm, frequencyX, frequencyY, ...
        tieIntensityOffset, fusionPixelSize, fusionCutoff, fusionOrder, ...
        fusionRegularization)
% Evaluate sensitivity to source-sampling density.

    fprintf('\n============================================================\n');
    fprintf(' SOURCE-SAMPLING CONVERGENCE TEST\n');
    fprintf('============================================================\n');

    nTests = numel(sampleTargets);
    reconstructedPhases = cell(nTests, 1);

    for testIndex = 1:nTests
        targetSamples = sampleTargets(testIndex);
        fprintf('\nRunning targetSourceSamples = %d\n', targetSamples);

        phaseDPC = reconstructDPC( ...
            objectSpectrum, sourceAnnular, pupil, defocusTransfer, ...
            frequencyAxis, frequencySpacing, dpcRegularization, targetSamples);
        [~, intensityFocus, intensityDefocus] = runForwardModel( ...
            objectSpectrum, sourceTIE, pupil, defocusTransfer, ...
            frequencyAxis, frequencySpacing, targetSamples);
        phaseTIE = reconstructTIE( ...
            intensityFocus, intensityDefocus, defocusUm, lambdaUm, ...
            frequencyX, frequencyY, tieIntensityOffset);

        reconstructedPhases{testIndex} = fuseTIEDPC( ...
            'PhaseTIE', phaseTIE, ...
            'PhaseDPC', phaseDPC, ...
            'PixelSize', fusionPixelSize, ...
            'Cutoff', fusionCutoff, ...
            'Order', fusionOrder, ...
            'Regularization', fusionRegularization);
    end

    reference = reconstructedPhases{end};
    reference = reference - mean(reference(:));
    rmse = zeros(nTests, 1);
    maximumAbsoluteDifference = zeros(nTests, 1);
    correlation = zeros(nTests, 1);

    for testIndex = 1:nTests
        currentPhase = reconstructedPhases{testIndex};
        currentPhase = currentPhase - mean(currentPhase(:));
        difference = currentPhase - reference;
        rmse(testIndex) = sqrt(mean(difference(:).^2));
        maximumAbsoluteDifference(testIndex) = max(abs(difference(:)));

        if testIndex == nTests
            correlation(testIndex) = 1;
        else
            correlationMatrix = corrcoef(currentPhase(:), reference(:));
            correlation(testIndex) = correlationMatrix(1, 2);
        end
    end

    resultsTable = table( ...
        sampleTargets(:), rmse, maximumAbsoluteDifference, correlation, ...
        'VariableNames', { ...
            'TargetSourceSamples', ...
            'RMSE_vs_320_rad', ...
            'MaxAbsDifference_vs_320_rad', ...
            'Correlation_vs_320'});

    referenceIndex = find(sampleTargets == 80, 1);
    if ~isempty(referenceIndex)
        referenceDifference = ...
            reconstructedPhases{referenceIndex} - manuscriptPhase;
        fprintf('\nVerification of targetSourceSamples = 80:\n');
        fprintf('Maximum absolute difference = %.16g rad\n', ...
            max(abs(referenceDifference(:))));
        fprintf('RMSE                        = %.16g rad\n', ...
            sqrt(mean(referenceDifference(:).^2)));
    end

    figure('Color', 'w', 'Name', 'Source-sampling convergence');
    yyaxis left;
    semilogy(sampleTargets, rmse, 'o-', 'LineWidth', 1.5);
    ylabel(sprintf('RMSE relative to target = %d (rad)', sampleTargets(end)));
    yyaxis right;
    plot(sampleTargets, correlation, 's-', 'LineWidth', 1.5);
    ylabel(sprintf('Correlation with target = %d', sampleTargets(end)));
    xlabel('Target source samples');
    title('Partially coherent source-sampling convergence');
    grid on;
end

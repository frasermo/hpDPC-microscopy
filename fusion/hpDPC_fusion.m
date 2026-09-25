function [phi_fused, opts_out] = hpDPC_fusion(varargin)
%HPDPC_FUSION Fuse precomputed TIE and pDPC phase maps in the Fourier domain.
%
%   PHI_FUSED = HPDPC_FUSION( ...
%       'phiTIE', PHI_TIE, ...
%       'phiDPC', PHI_DPC, ...
%       'dx', DX, ...
%       'k0', 0.05, ...
%       'alpha', 2, ...
%       'reg_eps', 1e-3)
%
% combines low-spatial-frequency information from a transport-of-intensity
% equation (TIE) phase reconstruction with complementary higher-spatial-
% frequency information from a polarization differential phase-contrast
% (pDPC) phase reconstruction.
%
% The two input phase maps must have identical dimensions and should be
% co-registered, expressed with a consistent phase sign convention, and
% sampled on the same computational grid before fusion.
%
% Name-value inputs
% -----------------
% phiTIE      TIE phase map (required).
% phiDPC      pDPC phase map of the same size (required).
% dx          Equal spatial sampling interval along both image axes.
%             Because the radial frequency is normalized by its maximum,
%             dx does not alter the weighting profile; it is retained for
%             physical-frequency diagnostics. Default: 1.
% k0          Dimensionless radial transition frequency. Default: 0.05.
% alpha       Spectral-transition exponent. Default: 2.
% reg_eps     Weak frequency-dependent regularization. Default: 1e-3.
% return_opts Return diagnostic spectra and weights when true. Default:
%             false.
%
% Outputs
% -------
% phi_fused   Fused phase reconstruction, returned as double precision.
% opts_out    Diagnostic structure when return_opts is true; otherwise an
%             empty structure.
%
% Method
% ------
% For normalized radial spatial frequency rho = K/Kmax,
%
%   W_TIE = exp(-(rho/k0).^alpha),    W_DPC = 1 - W_TIE.
%
% The weighted spectrum is weakly regularized before inverse transformation.
% The arbitrary global phase offset is matched to the mean TIE phase.
%
% Accompanying software for:
% "Hybrid Transport-of-Intensity and Polarization Differential Phase
% Contrast for Extended Spatial-Frequency Phase Imaging"
%
% Copyright (c) 2026 Fraser Montandon and contributors
% SPDX-License-Identifier: MIT

%% Parse inputs

parser = inputParser;
parser.FunctionName = mfilename;

addParameter(parser, 'phiTIE', [], ...
    @(x) isnumeric(x) && isreal(x) && ismatrix(x));
addParameter(parser, 'phiDPC', [], ...
    @(x) isnumeric(x) && isreal(x) && ismatrix(x));
addParameter(parser, 'dx', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'k0', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'alpha', 2, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'reg_eps', 1e-3, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(parser, 'return_opts', false, ...
    @(x) islogical(x) && isscalar(x));

parse(parser, varargin{:});
options = parser.Results;

%% Validate phase inputs

if isempty(options.phiTIE)
    error('hpDPC_fusion:MissingTIE', 'phiTIE must be provided.');
end
if isempty(options.phiDPC)
    error('hpDPC_fusion:MissingDPC', 'phiDPC must be provided.');
end
if ~isequal(size(options.phiTIE), size(options.phiDPC))
    error('hpDPC_fusion:SizeMismatch', ...
        'phiTIE and phiDPC must have identical dimensions.');
end
if any(~isfinite(options.phiTIE(:))) || ...
        any(~isfinite(options.phiDPC(:)))
    error('hpDPC_fusion:NonFiniteInput', ...
        'phiTIE and phiDPC must contain finite values only.');
end

phiTIE = double(options.phiTIE);
phiDPC = double(options.phiDPC);
[nRows, nColumns] = size(phiTIE);

if nRows < 2 || nColumns < 2
    error('hpDPC_fusion:InputTooSmall', ...
        'Phase maps must contain at least two rows and two columns.');
end

%% Construct the Fourier-frequency grid

rowFrequency = ifftshift( ...
    (-floor(nRows/2):ceil(nRows/2)-1).' / (nRows * options.dx));
columnFrequency = ifftshift( ...
    (-floor(nColumns/2):ceil(nColumns/2)-1) / ...
    (nColumns * options.dx));

% Implicit expansion avoids storing separate full 2-D X/Y frequency grids.
radialFrequency = hypot(rowFrequency, columnFrequency);
maximumFrequency = max(radialFrequency(:));

if ~isfinite(maximumFrequency) || maximumFrequency <= 0
    error('hpDPC_fusion:InvalidFrequencyGrid', ...
        'Unable to construct a valid Fourier-frequency grid.');
end

normalizedFrequency = radialFrequency / maximumFrequency;

%% Construct complementary spectral weights

tieWeight = exp(-(normalizedFrequency / options.k0).^options.alpha);
tieWeight = tieWeight / max(tieWeight(:));
dpcWeight = 1 - tieWeight;


dpcWeight(normalizedFrequency < 1e-6) = 0;

%% Blend and regularize the phase spectra

if options.return_opts
    tieSpectrum = fft2(phiTIE);
    dpcSpectrum = fft2(phiDPC);
    fusedSpectrum = tieWeight .* tieSpectrum + dpcWeight .* dpcSpectrum;
else

    fusedSpectrum = tieWeight .* fft2(phiTIE);
    fusedSpectrum = fusedSpectrum + dpcWeight .* fft2(phiDPC);
end

fusedSpectrum = fusedSpectrum ./ (1 + options.reg_eps .* ...
    (normalizedFrequency / 0.99).^2);

%% Return to the spatial domain and set the phase-offset convention

phi_fused = real(ifft2(fusedSpectrum));
phi_fused = phi_fused - mean(phi_fused(:)) + mean(phiTIE(:));

%% Optional diagnostic output

opts_out = struct();
if options.return_opts
    opts_out.W_TIE = tieWeight;
    opts_out.W_DPC = dpcWeight;
    opts_out.FT_TIE = tieSpectrum;
    opts_out.FT_DPC = dpcSpectrum;
    opts_out.PhiF = fusedSpectrum;
    opts_out.K = radialFrequency;
    opts_out.Kmax = maximumFrequency;
    opts_out.Knorm = normalizedFrequency;
    opts_out.k0 = options.k0;
    opts_out.alpha = options.alpha;
    opts_out.reg_eps = options.reg_eps;
    opts_out.dx = options.dx;
end
end

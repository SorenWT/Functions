function correlationMats = single_subject(dataFile,parcelFile,outDir,sessionName)
    % Example analysis for a single subject
    % 
    % INPUTS
    % - dataFile: A single meeg objects or filenames for an meeg object
    % - parcelFile: choose a binary ROI map. Take care that the resolution of the 
    %               nifti file matches that of the source reconstruction.
    % - outdir: choose a results directory
    % - sessionName: Optionally specify name of sessions (one for each entry in Dlist)

    if nargin < 4 || isempty(sessionName) 
        sessionName = 'single_subject';
    end
    
    if isa(dataFile,'meeg')
        assert(length(dataFile) == 1,'Only one input file is supported');
    end

    % set a save file name
    resultsName = fullfile(outDir, 'myResults');

    % setup the ROI network settings
    Settings = struct();
    Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
    Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
    Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
    Settings.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
    Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
    Settings.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
    Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
    Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
    Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
    Settings.EnvelopeParams.windowLength = 3; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
    Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
    Settings.EnvelopeParams.takeLogs  = true;                           % perform analysis on logarithm of envelope. This improves normality assumption
    Settings.frequencyBands           = {[1.3 4],[4 8],[8 13],[13 30],[30 85],[]};          % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
    Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA',  'peakVoxel' or 'spatialBasis'
    Settings.outputDirectory          = outDir;                         % Set a directory for the results output
    Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
    Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
    Settings.sessionName              = sessionName; 
    Settings.SaveCorrected            = struct('timeCourses',   true, ...  % save corrected timecourses
                                               'envelopes',     true,  ...  % save corrected power envelopes
                                               'variances',     true);     % save mean power in each ROI before correction
                              
    % Run the ROI network analysis
    Settings        = ROInets.check_inputs(Settings);
    correlationMats = ROInets_IRASA(dataFile, Settings, resultsName);

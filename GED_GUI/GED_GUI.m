% The authors of this code are Adrija Chatterjee and Prof. Pragathi P.
% Balasubramani, Translational Neuroscience and Technology Lab(Transit),
% Department of Cognitive Science, Indian Institute of Technology, Kanpur.

% References: 
% 1. https://github.com/mikexcohen/GED_tutorial
% 2. https://doi.org/10.1016/j.neuroimage.2021.118809              

% Please add the topoplotindie.m file in your path- which is shared in the
% GED_GUI folder (LINE 27).

% It will open a GUI-interface where the two set files(2 conditions) can be uploaded-
% Here, they are coded as baseline and experimental. 
% Generalised Eigen Decomposition (GED) spatially localise the EEG signals
% that orthogonally differentiates the two conditions. 

% contact: adrijac23@iitk.ac.in, pbalasub@iitk.ac.in 

%%


function GED_GUI()

    if isdeployed
        addpath(fullfile(ctfroot, 'GED_tutorial-main'));
    else
        % addpath('ADD PATH'); 
        %the code will work after adding the topoplotindie.m file.
    end

    % Create the main GUI window
    fig = uifigure('Name', 'GED Analysis', 'Position', [100, 100, 900, 600]);
    layout = uigridlayout(fig, [4, 3], 'ColumnWidth', {'1x', '1x', '1x'}, 'RowHeight', {'fit', 'fit', '2x', '2x'});

    % Baseline file selection
    uilabel(layout, 'Text', 'Baseline File(.set):', 'HorizontalAlignment', 'right');
    baselineLabel = uilabel(layout, 'Text', 'Not Selected', 'HorizontalAlignment', 'left');
    baselineButton = uibutton(layout, 'Text', 'Select File', 'ButtonPushedFcn', @(btn, event) selectFile('baseline'));

    % Experiment file selection
    uilabel(layout, 'Text', 'Experiment File(.set):', 'HorizontalAlignment', 'right');
    experimentLabel = uilabel(layout, 'Text', 'Not Selected', 'HorizontalAlignment', 'left');
    experimentButton = uibutton(layout, 'Text', 'Select File', 'ButtonPushedFcn', @(btn, event) selectFile('experiment'));

    % Run analysis button
    layoutbutton =  
    runButton = uibutton(layout, ...
    'Text', 'Run Analysis', ...
    'FontSize', 12, ...              % Smaller font size
    'BackgroundColor', [0, 0.447, 0.741], ... % Blue color (RGB for MATLAB blue)
    'FontColor', 'white', ...        % White text for contrast
    'ButtonPushedFcn', @(btn, event) runAnalysis()); 

    % Axes for displaying results
    axesR = uiaxes(layout);
    axesR.Title.String = 'R^-1 Matrix';
    
    axesS = uiaxes(layout);
    axesS.Title.String = 'S Matrix';
    
    axesDiff = uiaxes(layout);
    axesDiff.Title.String = 'R^-1 S-DIFF Matrix';
    
    topoplotAxes = uiaxes(layout);
    topoplotAxes.Title.String = 'Topography';
    % colormap(topoplotAxes, parula); 


    % File paths (stored as persistent variables)
    persistent baselineFile experimentFile;
    baselineFile = '';
    experimentFile = '';

    % Select File Callback
    function selectFile(type)
        [file, path] = uigetfile('*.set', 'Select .set File');
        if isequal(file, 0)
            return;
        end
        if strcmp(type, 'baseline')
            baselineFile = fullfile(path, file);
            baselineLabel.Text = baselineFile;
        elseif strcmp(type, 'experiment')
            experimentFile = fullfile(path, file);
            experimentLabel.Text = experimentFile;
        end
    end

    % Run Analysis Callback
    function runAnalysis()
        % Validate inputs
        if isempty(baselineFile)
            uialert(fig, 'Please select a baseline .set file.', 'Error');
            return;
        end
        if isempty(experimentFile)
            uialert(fig, 'Please select an experiment .set file.', 'Error');
            return;
        end

        % Load data
        baseline = pop_loadset(baselineFile);
        diff_trial = pop_loadset(experimentFile);

        % Process covariance matrices
        tmpd_R = preprocessData(baseline.data);
        covR = (tmpd_R * tmpd_R') / size(tmpd_R, 2);

        tmpd_S = preprocessData(diff_trial.data);
        covS = (tmpd_S * tmpd_S') / size(tmpd_S, 2);

        % Perform GED
        [evecs, evals] = eig(covS, covR);
        [~, sidx] = sort(diag(evals), 'descend');
        evecs = evecs(:, sidx);

        % Source separation
        compmap = evecs(:, 1)' * covS;
        [~, se] = max(abs(compmap));
        compmap = compmap * sign(compmap(se));
        chanloc = baseline.chanlocs;

        % Display Covariance Matrices
        imagesc(axesR, inv(covR));
        colorbar(axesR);

        imagesc(axesS, covS);
        colorbar(axesS);

        imagesc(axesDiff, inv(covR) * covS);
        colorbar(axesDiff);
        %%

        % Display Topography
        fig.HandleVisibility = 'on'; %this is the additional change required

        axes(topoplotAxes);
        topoplotIndie(compmap, chanloc, 'numcontour', 0);
        % colorbar(topoplotAxes);
        % scatter(topoplotAxes,rand(1,10),rand(1,10)) % this is working 
        colorbar
        % Completion Message
        % uialert(fig, 'Analysis Complete! Results displayed.', 'Success');
    end

    % Preprocess Data Function
    function data = preprocessData(rawData)
        mean_vals = mean(rawData, 2);
        std_vals = std(rawData, 0, 2);
        data = (rawData - mean_vals) ./ std_vals;
    end
end

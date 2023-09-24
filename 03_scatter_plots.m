% In this script, six scatter plots are generated for the following pairs
%  - bp and R_{-1.9482}
%  - bp and SCI_{-3.6624}
%  - ΔH and R_{-1.2383}
%  - ΔH and SCI_{-2.3554}
% ... with their respective regression lines.

close all; % Before drawing, close any figures already opened
clear;     % Clear all variables

% CONFIG: Line width and font fize for each curve in drawn figures
lineWidth = 2;
fontSize = 26;
% Save plots to images? Set to true if yes
saveToFile = false;
% Note: Dimensions of resulting images are scaled according to each window size.
%       To control the dimensions, after running the whole script, resize each
%       ... figure window and then run only the saveas functions
%       ... manually, by selection, at the end of this script

% Utility: Below is a function to round-off to 4 decimal places | returns string
%          Need to use this function as round(X,4,Type) does not exist in Octave
%          ... and sprintf("%.04f",X) does not round properly for some numbers.
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing Entropy and Heat Capacity of 22 lower benzenoids
expData = {reshape([% Boiling point
   80.1 218 338 340 431 425 429 440 496 493 497
   547  542 535 535 531 519 590 596 594 595 393
   ]', 22, 1), reshape([ % Enthalpy of Formation
   75.2 141   202.7 222.6 271.1 277.1 275.1 310.5 296 289.9 319.2
   323  301.2 348   335   336.3 336.9 296.7 375.6 366 393.3 221.3
]', 22, 1),
   "bp", "Δ𝐻^0_f" % Their labels
};

% 22 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 7 6 8 7  9 6  7  8 8 6  7  9  8  8  9  6  8  8  9  6
  0 4 6 8 8 10 6 12 10 8 8 12 10 10 12 12 10 12 12 12 10 8
  0 1 3 2 5 4  6 3  7  8 8 9  10  7 6  6  7  12  9  9 10 5
]'; % Used for computing indices based on edge endpoint degree partitions

% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(a) | n=1:R_a, n=2:SCI_a, a = alpha
getIndexFns = { % obtains a 22 by 1 matrix containing indices of the benzenoids
  @(a) (sum(d_f.*[4,6,9].^a,2)); % General Randic index
  @(a) (sum(d_f.*[4,5,6].^a,2)); % General SCI
}';
% Cell containing their respective labels
indexName = {"R" "SCI"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two
numIndices = size(getIndexFns,2); % two
numCurves = numData*numIndices;   % four

for edn = 1:numData % edn = experimental data number | 1=bp, 2=ΔH
  for fnn = 1:numIndices % fnn = function number | 1=R_a, 2=SCI_a
    ccFn = @(alpha) corrcoef( % Gets corrcoef between benzenoid prop and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);

    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));

    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);         % For the regression line
    x = [getIndexFns{fnn}(peakAlpha)(1) max(getIndexFns{fnn}(peakAlpha))];
    y = m*x + b;

    % Scatter plot
    this_figure = figure(3*(edn-1)+fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth', lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s_{−%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});

    % To find the standard error of fit
    % Calculate residuals
    residuals = expData{1, edn} - (m * getIndexFns{fnn}(peakAlpha) + b);
    % Calculate the sum of squared residuals (SSE)
    sse = sum(residuals.^2);
    % Calculate degrees of freedom (df)
    n = length(expData{1, edn});
    df = n - 2;  % Two parameters estimated: slope and intercept
    % Calculate the standard error of fit (SE)
    se = sqrt(sse / df);
    fprintf('Standard Error of Fit for %s and %s: %.4f\n', ...
    bestIndexLabel, expData{2, edn}, se);


    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");

    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)

    drawnow;
  end
end

if saveToFile
  % Also enforce more suitable axes for figures 2 and 5
  saveas(figure(1), "03_scatter_E_R.png");
  figure(2);
  axis([0.035 0.1125 250 600]);
  saveas(figure(2), "03_scatter_E_SCI.png");
  saveas(figure(3), "03_scatter_DH_R.png");
  figure(4);
  axis([0.2 0.9625 50 400]);
  saveas(figure(4), "03_scatter_DH_SCI.png");
end

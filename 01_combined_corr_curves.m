% In this script, six values are closely approximated via golden section search
%  - α value for which correlation coefficient ρ is strongest between E and R_a
%  - α value for which ρ is strongest between bp and SCI_α
%  - α value for which ρ is strongest between ΔH and R_α
%  - α value for which ρ is strongest between ΔH and SCI_α
% Additionally, curves for ρ against α near these values are plotted in 4 figs

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
pkg load statistics;

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2; % the line thickness in the graph
fontSize = 16; % the font size in the graph
saveToFile = true; % Set to true to auto-save plots
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
  @(a) (sum(d_f.*[4,6,9].^a,2)); % General Randic index(R)
  @(a) (sum(d_f.*[4,5,6].^a,2)); % General Sum connectivity index(SCI)
}';
% Cell containing their respective labels
indexName = {"R" "SCI"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two (bp & ΔH)
numIndices = size(getIndexFns,2); % two (R & SCI)
numCurves = numData*numIndices;   % four (4 curves)

% All x in visible ranges (both plots - near and far)
xnear = [linspace(-0.63, 0.2, 800); linspace(-3.30, 0.25, 800)];
xfar = linspace(-20,20,800); % xfar is same for both bp and ΔH

% Do the same procedure for each experimental data i.e., bp and ΔH
for ii = 1:numData
  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;

  % Indicate inverval for which R_a is the better indicator (before corr-curves)

  % WARNING: these xmeet1-ymeet2 values are hardcoded, computed separately
  %          ... These are coordinates where ρ-α of {ΔH or E}-R_α
  %          ... intersects ρ-α of {ΔH or E}-SCI_α
  xmeet1 = [
    -0.46299; % for bp
    -2.2887   % for ΔH
  ](ii); % of these 2, either the first value or second is used, depending on ii
  ymeet1 = [
    0.99641; % for bp
    0.97216  % for ΔH
  ](ii);
  xmeet2 = [
    -0.00045;    % for bp
    -0.00010081  % for ΔH
  ](ii);
  ymeet2 = [
    0.99574; % for bp
    0.94512  % for ΔH
  ](ii);

  ybox = [
    0.9968; % for figure 3 (blue dotted line)
    0.98  % for figure 4 (blue dotted line)
  ](ii);
  % Plot the blue dashed box (before drawing the curves so it appear beneath)
  plot([xmeet1 xmeet1 0 0], [0 ybox ybox 0], '--b', 'LineWidth', lineWidth);

  yend = 0; % <-- to be assigned some value later for adjusting visible range

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between bp/ΔH (depending on ii) with specified α
    %                                and R_a/SCI_a (depending on n)
    get_indices_vals = @(alpha) getIndexFns{n}(alpha)(!isnan(expData{1,ii}));
    ccFn = @(alpha) corrcoef(
      get_indices_vals(alpha), % Either R(alpha) or SCI(alpha)
      expData{1,ii}(!isnan(expData{1,ii})) % bp/ΔH
    )(1,2);

    % generate corresponding y values
    ynear = arrayfun(ccFn, xnear(ii,:));
    yfar = arrayfun(ccFn, xfar);

    % Compute peak values via. golden section search, and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(GoldenSectionSearch_Maximum(ccFn, xnear(ii,1), xnear(ii,end), 1e-15))
    peakCorrCoeff = ccFn(peakAlpha)

    X_lm = getIndexFns{n}(peakAlpha);
    Y_lm = expData{1,ii};
    mdl = fitlm(X_lm, Y_lm);
    close

    % Calc confidence intervals
    grad_coef = mdl{3,2}
    grad_lower = mdl{3,4};
    CI_grad = grad_coef - grad_lower
    intercept_coef = mdl{2,2}
    intercept_lower_CI = mdl{2,4};
    CI_intercept = intercept_coef - intercept_lower_CI
    disp("")

    % Generate curve label                  [bp,ΔH]         [R,SCI]
    curveLabels{n} = sprintf("%s and %s_a", expData{2,ii}, indexName{n});

    figure(ii); % One zoomed-out plot for each expData
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curves(n) = plot(xnear(ii,:), ynear, '-', 'LineWidth', lineWidth);

    % Show the peak in the plot: draw indicator lines & display coordinates
    plot([peakAlpha peakAlpha xnear(ii,1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2); % Black dashed indicator lines
    text(peakAlpha, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakAlpha)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');
        % Negative sign entered manually here to bypass the default
        % ... usage of "hypen-minus" instead of "minus" (− vs -)

    yend = max(yend, ynear(end)); % y value to be used as visible y lower bound
  end

  % Mark and write on the plot the limits of alpha where R_a is better than SCI_a
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b',
       'MarkerSize', 16, 'LineWidth', lineWidth/1.5); % Mark with blue asterisks
  text(xmeet1, ymeet1, {'', sprintf(" (−%s, %s)", as_4_dp_str(abs(xmeet1)), as_4_dp_str(ymeet1))},
       'VerticalAlignment', 'top', 'Color', [0, 0, 0.8]); % Write blue text
  text(xmeet2, ymeet2, {'', sprintf("(0, %s) ", as_4_dp_str(ymeet2))},
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);

  % Label this expData's zoomed-in plot
  xlabel('α');
  ylabel('ρ');
  leg = legend(curves, curveLabels); % curves contains all drawn "xnear" curves
  set(leg, 'location', "southwest"); % the location of the legend box

  ybox_space = [
    0.000005; % spacing between upper y-value with the blue dotted lines (figure 3)
    0.0005  % spacing between upper y-value with the blue dotted lines (figure 4)
  ](ii);
  axis([xnear(ii,1) xnear(ii,end) yend ybox+ybox_space]); % Enforce figure's visible range
  drawnow;

  % Label the zoomed-out plot
  figure(ii);
  xlabel('α');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels); % curvesFar contains all drawn "xfar" curves
  set(leg, 'location', "southeast");

  if ii==2
    set(leg, 'location', "southeast");

  end

  hold off;
end

for ii = 1:4
  % Replace hyphens with minuses on negative axes
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));

  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_bp_indices_FAR.png");
  saveas(figure(2), "01_comb_ccurves_DH_indices_FAR.png");
  saveas(figure(3), "01_comb_ccurves_bp_indices_NEAR.png");
  saveas(figure(4), "01_comb_ccurves_DH_indices_NEAR.png");
end

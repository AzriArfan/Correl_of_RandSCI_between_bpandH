# Correl_of_RandSCI_between_bpandH
<h3>Files</h3>
<table>
  <tr><th>File Name</th><th>Description</th></tr>
  <tr><td><a href=01_combined_corr_curves.m>01_combined_corr_curves.m</a></td><td>Generates combined correlation curves plots, grouped by experimental data, for the 4 pairs of data, highlighting the curves' peaks and intersection between two curves with highest correlation in some given alpha</td></tr>
  <tr><td><a href=02_good_alpha_intervals.m>02_good_alpha_intervals.m</a></td><td>Generates 4 correlation curve plots, one for each pair of data, highlighting the alpha intervals which results in good correlation</td></tr>
  <tr><td><a href=03_scatter_plots.m>03_scatter_plots.m</a></td><td>Generates 4 scatter plots, one for each pair of data</td></tr>
  <tr><td><a href=GoldenSectionSearch_Maximum.m>GoldenSectionSearch_Maximum.m</a></td><td><b>IMPORTANT: This must be placed in the working directory as it is referenced in the three other scripts.</b> This is the Octave 8.1 implementation of the <a href="https://en.wikipedia.org/wiki/Golden-section_search">Golden Section Search</a> algorithm. I manually translated this from the <a href="https://en.wikipedia.org/wiki/Golden-section_search">Python code</a>.</td></tr>
</table>
<h3>Syntax</h3>
<p>The above scripts are to be run with the <a href=https://octave.org/>GNU Octave</a> software, not MATLAB. Due to some syntax differences between MATLAB and GNU Octave, e.g. in subsetting matrices/arrays, the scripts will not run (without ammendments) in MATLAB.</p>
<h3>Expected Output</h3>
<h4 align=center>Command Window output (<a href=01_combined_corr_curves.m>1<sup>st</sup> script</a>)</h4>

```
bp against general_R
peakAlpha = -0.330300790051691
peakCorrCoeff = 0.996557927931623
bp against general_SCI
peakAlpha = -0.623360326403817
peakCorrCoeff = 0.996463128711537
ΔH against general_R
peakAlpha = -1.720555062022164
peakCorrCoeff = 0.977343521423170
ΔH against general_SCI
peakAlpha = -3.249655284970219
peakCorrCoeff = 0.975715487711226

```

<h4 align=center>Command Window output (<a href=01_combined_corr_curves.m>2<sup>nd</sup> script</a>)</h4>

```
ρ(bp,R_α) ≥ 0.0.99656 when α ∈ [1.11542, 0.6813]
ρ(bp,SCI_α) ≥ 0.99646 when α ∈ [2.2181, 1.4401]
ρ(ΔH,SO_α) ≥ 0.97734 when α ∈ [2.4965, 0.8520]
ρ(ΔH,R_α) ≥ 0.97572 when α ∈ [4.5788, 1.7437]

```

<h3>Combine correlation curves</h3>

![fargraph_bp](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/f72b55d4-da37-49a1-811b-46812572e3e2)

![fargraph_DH](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/8912b6a4-7858-4ad2-be42-392fd18d6951)

![neargraph_bp](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/82e7da79-5f00-4694-99ee-c15c598d61d2)

![neargraph_DH](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/7850fea1-c5c3-4a12-b2c3-27219a9c1a0e)

<h3>Good alpha intervals</h3>

![EG_bp_and_Rnew](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/c82966a1-412a-45cb-a76a-7967ec28c643)

![EG_bp_and_SCInew](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/a39411bb-fc9f-4daf-b44d-cdd3f092cf3c)

![EG_DH_and_Rnew](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/a29aca9d-d103-4bb6-bd3f-bb65fc593a77)

![EG_DH_and_SCInew](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/68a7ce7f-5586-4934-b785-c1649e1fa82c)

<h3>Scatter plots</h3>

![SP_bp_and_R](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/052c3154-c2eb-4aae-ba43-a23aa3c217ec)

![SP_bp_and_SCI](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/a04cc1d8-c839-4acd-b894-3895089ca782)

![SP_DH_and_R](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/9406a628-ca23-4de3-95c3-d76143f3556e)

![SP_DH_and_SCI](https://github.com/AzriArfan/Correl_of_RandSCI_between_bpandH/assets/145908475/0f69d683-71be-4d62-87ff-ef5144d97c50)



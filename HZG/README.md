# H→Zγ Run 3 mass fits

This branch adds the HZG signal worker and native background fitter to the CMS HiggsDNA FinalFit base. Analysis inputs, fitted parameters and job submission live with the campaign. The HZG workers read a RooWorkspace named `CMS_hzg_workspace` with observable `CMS_hzg_mass`.

## Environment

Use the HZG CombinedLimit build that exports the complete erf PDFs and ROOT dictionaries. The public dependency is [MingxuHeseri/HiggsAnalysis-CombinedLimit at 69b5f15](https://github.com/MingxuHeseri/HiggsAnalysis-CombinedLimit/tree/69b5f154d7a2992bbde67ae1f9e62c6cb4717f49). Activate the corresponding ROOT and Combine runtime before starting a worker. The signal worker uses PyROOT, NumPy, SciPy and pandas. The background build uses a C++ compiler, ROOT with RooFit, Boost and this checkout's top-level tdrStyle sources.

## Signal fit

Run from the repository root with the production input directory, category and seed configuration. In this example the input file has an explicit pooled process named `all`, and the frozen fit uses the five named partitions.

```bash
python3 HZG/Signal/scripts/signalFit.py \
  --inputWSDir "$signal_inputs" --outputPath "$output" --ext run3 \
  --proc all --sourceProcs all --cat "$category" --flav ele --year all \
  --massPoints 125 --useDCB --skipSystematics --skipBeamspotReweigh \
  --fitMassLow 110 --fitMassHigh 140 --nBins 340 --MHPolyOrder 0 \
  --fitPartitions 0,1,2,3,4 --fitPartitionVar fit_partition \
  --fitSeedConfig "$seed_config" --fitObjective chi2 \
  --minimizerMethod L-BFGS-B --minimizerTolerance 1e-7 \
  --runMinuitPolish --strictEdmThreshold 1e-6 --minuitPolishTolerance 0.0001
```

Use `--flav mu` for muons. Add `--useDCBplusGaussian` to `--useDCB` to fit the shared-mean DSCB plus Gaussian. Supply actual source process names when merging separate datasets. The signed nominal sample supplies the normalization, and the fit partitions select the events used for shape estimation. For the unit-factor pooled process the workspace uses `1000 × PDF_norm = nominal sum of weights`; preserve the corresponding card rate. Standalone per-process factors are defined in `Signal/tools/XSBRMap.py`. The saved workspace opens at the nominal `MH = 125 GeV`. The worker writes `wsig_13TeV` below `outdir_<ext>/signalFit/output`.

To include signal systematics, omit `--skipSystematics` and pass `--systematicsFile /path/to/category.pkl` with the desired `--scales` and `--smears`. The default table is `<outputPath>/outdir_<ext>/calcPhotonSyst/pkl/<category>.pkl`. The worker also recognizes the older Signal-directory layouts. The DataFrame needs one matching row for the systematics process (including diagonal-process substitution) and, where present, the `cat` and `flav` columns. A table without `flav` uses the same unique row for both channels. For nuisance `N`, provide `N_mean`, `N_sigma` and `N_rate` as fractional shifts. Historical energy/type-qualified columns are supported. In a Combine card, declare the analysis-approved nuisance constraints and verify the saved workspace response at ±1σ. Run the table regression with `python3 HZG/Signal/test/test_systematics.py`.

## Native background fit

In the CMSSW runtime build and run the integral and plot-normalization checks with

```bash
make -C HZG/Background -j4 check
```

For a standalone installation, set `COMBINE_INCLUDE_ROOT` to the directory containing `HiggsAnalysis/CombinedLimit/interface`, and `COMBINE_LIBDIR` to the directory containing its library. Keep the installed Boost and ROOT headers and libraries visible to the compiler. For Apple Clang in a conda ROOT environment use `CXX=/usr/bin/clang++ LDFLAGS=-fuse-ld=/usr/bin/ld` when required by the SDK.

Example for a blinded sideband shape check with an explicitly selected family and order. Run from an empty output directory with writable paths.

```bash
export HZG_BACKGROUND_PARAM_CONFIG="$parameter_json"
"$finalfit/HZG/Background/bin/background_fit" \
  -i "$background_workspace" --isData 1 --flashggCats "$category" --singleCat 0 \
  --year all --blindFit --mgg_low "$low" --mgg_high "$high" \
  --funclist ExpBernsteinErfTurnOn --candidateOrders 2 \
  --skipSSTest --saveMultiPdf envelope.root --datfile fit.dat --outDir plots
```

The low-level CLI requires explicit `--blindFit` for sideband fitting. The default objective is the sideband NLL in that mode. `--dochi2Fit` selects the binned objective. To run the native spurious-signal test, supply `--sstInputPath` and `--sstSignalPath` and omit `--skipSSTest`. Scientific acceptance uses the study's family list, fit ranges and selection criteria.

For frozen erf envelopes, `Background/native_erf_envelope.py::native_envelope` reconstructs the accepted members with their recorded parameters, discrete index and penalties. Load the HZG CombinedLimit library and the built `libHZGBackground`, add `HZG/Background/interface` to the ROOT interpreter include path and declare `#include "PdfFactory.h"` before calling it. The returned `multipdf` workspace can be consumed by a shapes card and standard `text2workspace.py`.

## Source history

The upstream base is [cms-analysis/flashggFinalFit at f715b266](https://github.com/cms-analysis/flashggFinalFit/commit/f715b266dea7776d29353afab40fb7f815c3e7a0). The HZG branch includes both original tips of [FlashggFinalFit_HZGamma](https://github.com/Vvvvvvvictor/FlashggFinalFit_HZGamma) as merge parents, `41ac9c1` and `dc4ef86`. They preserve JackmanHe's and Mingtao's original commits, including signal pooling, turn-on models, and sideband fitting. The production port follows the later CERN14 signal and CERN16 background implementations. The source tree for HZG contains the worker dependencies; historical configurations, generated cards, binaries, trial scripts and outputs remain accessible in the original history.

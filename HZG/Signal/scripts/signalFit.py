print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZG SIGNAL FITTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
import ROOT
import pandas as pd
import pickle
import math
import os, sys
import json
from optparse import OptionParser
import glob
import re
from collections import OrderedDict as od

from pathlib import Path

hzg_dir = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(hzg_dir / "tools"), str(hzg_dir / "Signal" / "tools")]
if ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") < 0:
  raise RuntimeError("Load the CombinedLimit runtime before running the HZG signal fit")

from commonTools import *
from commonObjects import *

from replacementMap import globalReplacementMap
from XSBRMap import *
from simultaneousFit import *
from finalModel import *
from plottingTools import *

# Constant
MHLow, MHHigh = '120', '130'
MHNominal = '125'

def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HZG SIGNAL FITTER (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  sys.exit(1)

def get_options():
  parser = OptionParser()
  parser.add_option("--xvar", dest='xvar', default='CMS_hzg_mass', help="Observable to fit")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--outputPath", dest='outputPath', default='.', help="Output path")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--proc", dest='proc', default='', help="Signal process")
  parser.add_option("--sourceProcs", dest='sourceProcs', default='', help="Comma-separated source processes to merge for shape datasets")
  parser.add_option("--cat", dest='cat', default='', help="RECO category")
  parser.add_option("--flav", dest='flav', default='ele', help="Flavour")
  parser.add_option('--useAllFlavDataset', dest='useAllFlavDataset', default=False, action="store_true", help="When --flav all, read *_all_* datasets directly instead of merging ele,mu")
  parser.add_option("--year", dest='year', default='2016', help="Year")
  parser.add_option("--analysis", dest='analysis', default='fiducialAnalysis', help="Analysis handle: used to specify replacement map and XS*BR normalisations")
  parser.add_option('--massPoints', dest='massPoints', default='120,125,130', help="Mass points to fit")
  parser.add_option('--doEffAccFromJson', dest='doEffAccFromJson', default=False, action="store_true", help="Extract eff x acc from json (produced by getEffAcc). Else, extract from nominal weights in flashgg workspaces")
  parser.add_option('--skipBeamspotReweigh', dest='skipBeamspotReweigh', default=False, action="store_true", help="Compatibility flag for preweighted HZG inputs")
  parser.add_option('--doPlots', dest='doPlots', default=False, action="store_true", help="Produce Signal Fitting plots")
  parser.add_option("--doVoigtian", dest='doVoigtian', default=False, action="store_true", help="Use Voigtians instead of Gaussians for signal models with Higgs width as parameter")
  parser.add_option("--useDCB", dest='useDCB', default=False, action="store_true", help="Use DCB in signal fitting")
  parser.add_option("--useDCBplusGaussian", dest='useDCBplusGaussian', default=False, action="store_true", help="Add a broad Gaussian shoulder to the DCB signal model (requires --useDCB)")
  parser.add_option("--useDiagonalProcForShape", dest='useDiagonalProcForShape', default=False, action="store_true", help="Use shape of diagonal process, keeping normalisation (requires diagonal mapping produced by getDiagProc script)")
  parser.add_option('--skipVertexScenarioSplit', dest='skipVertexScenarioSplit', default=True, action="store_true", help="Compatibility flag for the inclusive HZG vertex model")
  parser.add_option('--skipZeroes', dest='skipZeroes', default=False, action="store_true", help="Skip proc x cat is numEntries = 0., or sumEntries < 0.")
  # For systematics
  parser.add_option('--skipSystematics', dest='skipSystematics', default=False, action="store_true", help="Skip shape systematics in signal model")
  parser.add_option('--systematicsFile', help='Category systematics DataFrame pickle (overrides the campaign path)')
  parser.add_option('--useDiagonalProcForSyst', dest='useDiagonalProcForSyst', default=False, action="store_true", help="Use diagonal process for systematics (requires diagonal mapping produced by getDiagProc script)")
  parser.add_option("--scales", dest='scales', default='CMS_scale_e,CMS_scale_g,CMS_scale_m', help="Photon shape systematics: scales")
  parser.add_option("--scalesCorr", dest='scalesCorr', default='', help='Photon shape systematics: scalesCorr')
  parser.add_option("--scalesGlobal", dest='scalesGlobal', default='', help='Photon shape systematics: scalesGlobal')
  parser.add_option("--smears", dest='smears', default='CMS_res_e,CMS_res_g,CMS_res_m', help='Photon shape systematics: smears')
  # Parameter values
  parser.add_option('--replacementThreshold', dest='replacementThreshold', default=100, type='int', help="Nevent threshold to trigger replacement dataset")
  parser.add_option('--MHPolyOrder', dest='MHPolyOrder', default=1, type='int', help="Order of polynomial for MH dependence")
  parser.add_option("--MHPolyConfig", dest='MHPolyConfig', default=None, help="Read coefficient value from config file")
  parser.add_option("--fitSeedConfig", dest='fitSeedConfig', default=None, help="JSON with floating signal-fit parameter starting values")
  parser.add_option('--nBins', dest='nBins', default=340, type='int', help="Number of bins for fit")
  parser.add_option('--fitMassLow', dest='fitMassLow', default=105., type='float', help="Lower mass bound for signal fit")
  parser.add_option('--fitMassHigh', dest='fitMassHigh', default=150., type='float', help="Upper mass bound for signal fit")
  parser.add_option('--fitObjective', dest='fitObjective', default='chi2', help="Fit objective: chi2 or unbinnedNLL")
  parser.add_option('--unbinnedWeightMode', dest='unbinnedWeightMode', default='signed', help="Unbinned NLL weight treatment: signed, positive, or absolute")
  parser.add_option('--fitPartitionVar', dest='fitPartitionVar', default='fit_partition', help="RooDataSet variable carrying the deterministic tune/validation partition")
  parser.add_option('--fitPartitions', dest='fitPartitions', default='', help="Comma-separated partition indices used for fitting; empty uses all events")
  # Minimizer options
  parser.add_option('--minimizerMethod', dest='minimizerMethod', default='TNC', help="(Scipy) Minimizer method")
  parser.add_option('--minimizerTolerance', dest='minimizerTolerance', default=1e-8, type='float', help="(Scipy) Minimizer toleranve")
  parser.add_option('--runMinuitPolish', dest='runMinuitPolish', default=False, action="store_true", help="Polish the identical objective with ROOT Minuit2 and record EDM/covariance diagnostics")
  parser.add_option('--strictEdmThreshold', dest='strictEdmThreshold', default=1e-6, type='float', help="Strict Minuit2 EDM acceptance threshold")
  parser.add_option('--minuitPolishTolerance', dest='minuitPolishTolerance', default=None, type='float', help="Optional Minuit2 stopping tolerance, independent of the strict EDM acceptance threshold")
  parser.add_option('--requireStrictConvergence', dest='requireStrictConvergence', default=False, action="store_true", help="Fail after writing diagnostics unless the strict ROOT Minuit2 gate passes")
  return parser.parse_args()
(opt,args) = get_options()

if opt.useDCBplusGaussian and not opt.useDCB:
  print(" --> [ERROR] --useDCBplusGaussian requires --useDCB")
  sys.exit(1)
if opt.fitMassLow >= opt.fitMassHigh:
  print(" --> [ERROR] fitMassLow (%.3f) must be below fitMassHigh (%.3f)"%(opt.fitMassLow,opt.fitMassHigh))
  sys.exit(1)
if opt.fitObjective not in ['chi2','unbinnedNLL']:
  print(" --> [ERROR] fitObjective must be chi2 or unbinnedNLL (got %s)"%opt.fitObjective)
  sys.exit(1)
if opt.unbinnedWeightMode not in ['signed','positive','absolute']:
  print(" --> [ERROR] unbinnedWeightMode must be signed, positive, or absolute (got %s)"%opt.unbinnedWeightMode)
  sys.exit(1)
if opt.strictEdmThreshold <= 0.:
  print(" --> [ERROR] strictEdmThreshold must be positive (got %g)"%opt.strictEdmThreshold)
  sys.exit(1)
if opt.requireStrictConvergence and not opt.runMinuitPolish:
  print(" --> [ERROR] --requireStrictConvergence requires --runMinuitPolish")
  sys.exit(1)
if opt.minuitPolishTolerance is not None and (not math.isfinite(opt.minuitPolishTolerance) or opt.minuitPolishTolerance <= 0.):
  raise ValueError("minuitPolishTolerance must be finite and positive")
outputBase = os.path.join(swd__, opt.outputPath)
fitPartitions = []
if opt.fitPartitions:
  try:
    fitPartitions = sorted(set(int(value) for value in opt.fitPartitions.split(",")))
  except ValueError:
    print(" --> [ERROR] fitPartitions must be comma-separated integers")
    sys.exit(1)
  if any(value < 0 for value in fitPartitions):
    print(" --> [ERROR] fitPartitions cannot contain negative values")
    sys.exit(1)

sourceProcsForFit = [p for p in (opt.sourceProcs.split(",") if opt.sourceProcs else opt.proc.split(",")) if p]

def findSignalWSFileNames(input_ws_dir, mass_point, proc):
  pattern = "%s/output*M%s*%s.root"%(input_ws_dir,mass_point,proc)
  files = glob.glob(pattern)
  if not files:
    print(" --> [ERROR] no signal workspace files match: %s"%pattern)
    sys.exit(1)
  return files

def getSignalDatasetName(proc, mass_point, flav, cat):
  return "%s_%s_%s_%s_%s"%(procToData(proc),mass_point,flav,sqrts__,cat)

def getSignalDatasetNameCandidates(proc, mass_point, flav, cat):
  names = [getSignalDatasetName(proc,mass_point,flav,cat)]
  rawName = "%s_%s_%s_%s_%s"%(proc,mass_point,flav,sqrts__,cat)
  if rawName not in names: names.append(rawName)
  return names

def loadReducedDataset(input_ws, dataset_names, reduced_set, source):
  if isinstance(dataset_names,str): dataset_names = [dataset_names]
  dataset = None
  for dataset_name in dataset_names:
    dataset = input_ws.data(dataset_name)
    if dataset: break
  if not dataset:
    print(" --> [ERROR] dataset(s) %s not found in %s"%(",".join(dataset_names),source))
    sys.exit(1)
  return dataset.reduce(reduced_set)

def inputFitArgSet(input_ws):
  mass = input_ws.var(opt.xvar)
  if not mass:
    raise RuntimeError("Input workspace is missing %s"%opt.xvar)
  result = ROOT.RooArgSet(mass)
  if fitPartitions:
    partition = input_ws.var(opt.fitPartitionVar)
    if not partition:
      raise RuntimeError(
        "Input workspace is missing requested fit partition variable %s"%
        opt.fitPartitionVar
      )
    result.add(partition)
  return result

def massOnlyDataset(dataset):
  mass = dataset.get().find(opt.xvar)
  if not mass:
    raise RuntimeError("Dataset %s is missing %s"%(dataset.GetName(),opt.xvar))
  return dataset.reduce(ROOT.RooArgSet(mass))

def selectFitDataset(dataset):
  if not fitPartitions:
    return massOnlyDataset(dataset)
  expression = "||".join(
    "%s==%d"%(opt.fitPartitionVar,value) for value in fitPartitions
  )
  mass = dataset.get().find(opt.xvar)
  selected = dataset.reduce(ROOT.RooArgSet(mass),expression)
  if not selected or selected.numEntries() == 0:
    raise RuntimeError(
      "Fit partition selection %s produced no entries for %s"%
      (expression,dataset.GetName())
    )
  return selected

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETUP: signal fit
print(" --> Running fit for (proc,cat,flav) = (%s,%s,%s)"%(opt.proc,opt.cat,opt.flav))
if( len(opt.massPoints.split(",")) == 1 )&( opt.MHPolyOrder > 0 ):
  print(" --> [WARNING] Attempting to fit polynomials of O(MH^%g) for single mass point. Setting order to 0"%opt.MHPolyOrder)
  opt.MHPolyOrder=0

# Add stopwatch function

# Load replacement map
if opt.analysis not in globalReplacementMap:
  print(" --> [ERROR] replacement map does not exist for analysis (%s). Please add to tools/replacementMap.py"%opt.analysis)
  leave()
else: rMap = globalReplacementMap[opt.analysis]

# Load XSBR map
if opt.analysis not in globalXSBRMap:
  print(" --> [ERROR] XS * BR map does not exist for analysis (%s). Please add to tools/XSBRMap.py"%opt.analysis)
  leave()
else: xsbrMap = globalXSBRMap[opt.analysis]

# Load RooRealVars
nominalWSFileName = findSignalWSFileNames(opt.inputWSDir,MHNominal,sourceProcsForFit[0])[0]
f0 = ROOT.TFile(nominalWSFileName,"read")
inputWS0 = f0.Get(inputWSName__)
if not inputWS0:
  print(" --> [ERROR] workspace %s not found in %s"%(inputWSName__,nominalWSFileName))
  sys.exit(1)
xvar = inputWS0.var(opt.xvar)
if not xvar:
  print(" --> [ERROR] observable %s not found in %s"%(opt.xvar,nominalWSFileName))
  sys.exit(1)
xvarFit = xvar.Clone()
aset = ROOT.RooArgSet(xvar)
f0.Close()

# Create MH var
MH = ROOT.RooRealVar("MH","m_{H}", int(MHLow), int(MHHigh))
MH.setUnit("GeV")
MH.setConstant(True)

if opt.skipZeroes:
  # Extract nominal mass dataset and see if entries == 0
  WSFileName = findSignalWSFileNames(opt.inputWSDir,MHNominal,opt.proc)[0]
  f = ROOT.TFile(WSFileName,"read")
  inputWS = f.Get(inputWSName__)
  if opt.proc.split("_")[-1] in ["i"]:
  # if opt.proc.split("_")[-1] in ["in", "out"]:
    d = inputWS.data("%s_%s_%s_%s_%s_%s"%(procToData(opt.proc.split("_")[0]),procToData(opt.proc.split("_")[-1]),MHNominal,opt.flav,sqrts__,opt.cat)).reduce(aset)
  else:
    d = loadReducedDataset(inputWS,getSignalDatasetName(opt.proc,MHNominal,opt.flav,opt.cat),aset,WSFileName)
  if( d.numEntries() == 0. )|( d.sumEntries <= 0. ):
    print(" --> (%s,%s,%s) has zero events. Will not construct signal model"%(opt.proc,opt.cat,opt.flav))
    sys.exit(1)
  inputWS.Delete()
  f.Close()
 
# Define proc x cat x flav with which to extract shape: if skipVertexScenarioSplit label all events as "RV"
procRVFit, catRVFit, flavRVFit = opt.proc, opt.cat, opt.flav
print(" --> Skipping vertex scenario split")

# Options for using diagonal process from getDiagProc output json
if opt.useDiagonalProcForShape:
  if not os.path.exists("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(outputBase,opt.ext)):
    print(" --> [ERROR] Diagonal process json from getDiagProc does not exist. Using nominal proc x cat x flav for shape")
  else:
    with open("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(outputBase,opt.ext),"r") as jf: dproc = json.load(jf)
    procRVFit = dproc[opt.cat]
    print(" --> Using diagonal proc (%s,%s,%s) for shape"%(procRVFit,opt.cat,opt.flav))

# Process for syst
procSyst = opt.proc
if opt.useDiagonalProcForSyst:
  if not os.path.exists("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(outputBase,opt.ext)):
    print(" --> [ERROR] Diagonal process json from getDiagProc does not exist. Using nominal proc x cat x flav for systematics")
  else:
    with open("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(outputBase,opt.ext),"r") as jf: dproc = json.load(jf)
    procSyst = dproc[opt.cat]
    print(" --> Using diagonal proc (%s,%s,%s) for systematics"%(procSyst,opt.cat,opt.flav))

# Define process with which to extract normalisation: nominal
# procNorm, catNorm, flavNorm = opt.proc, opt.cat, opt.flav

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT DATASETS TO FIT (for each mass point)
nominalDatasets = od()
# For RV (or if skipping vertex scenario split)
datasetRVForFit = od()
print((opt.inputWSDir))
print(" --> Source processes for shape datasets: %s"%(",".join(sourceProcsForFit)))
print(" --> Fit objective: %s"%opt.fitObjective)
if opt.fitObjective == 'unbinnedNLL': print(" --> Unbinned weight mode: %s"%opt.unbinnedWeightMode)
print(" --> Fit partitions: %s"%(",".join(str(value) for value in fitPartitions) if fitPartitions else "all"))
print(" --> Fit mass range: %.3f to %.3f GeV"%(opt.fitMassLow,opt.fitMassHigh))
for mp in opt.massPoints.split(","):
  d_list = []

  for proc in sourceProcsForFit:
    # Read the configured process label directly; e.g. proc=all reads all_* datasets.
    print((" --> Extracting dataset: %s/output*M%s*%s.root" %(opt.inputWSDir,mp,proc)))
    WSFileNames = findSignalWSFileNames(opt.inputWSDir,mp,proc)
    print(WSFileNames)

    flavsToRead = opt.flav.split(",") if opt.useAllFlavDataset else re.sub("all", "ele,mu", opt.flav).split(",")
    for flav in flavsToRead:
    # for flav in opt.flav.split(","):
      for WSFileName in WSFileNames:
        f = ROOT.TFile(WSFileName,"read")
        inputWS = f.Get(inputWSName__)
        # reducedSet = ROOT.RooArgSet(inputWS.var(opt.xvar),inputWS.var("dZ"))
        reducedSet = inputFitArgSet(inputWS)
        datasetName = getSignalDatasetNameCandidates(proc,mp,flav,catRVFit)
        # print(datasetName)
        rdataset = loadReducedDataset(inputWS,datasetName,reducedSet,WSFileName)
        d_list.append(rdataset)
        # inputWS.Delete()
        f.Close()
  d = d_list[0].emptyClone()
  if len(d_list) > 1: 
    for ds in d_list: d.append(ds)
  else: d = d_list[0]
  print(("Number of DataSet: %d"%len(d_list)))
  nominalDatasets[mp] = massOnlyDataset(d)
  fitDataset = selectFitDataset(d)
  datasetRVForFit[mp] = fitDataset
        
# Check if nominal yield > threshold (or if +ve sum of weights). If not then use replacement proc x cat
if( datasetRVForFit[MHNominal].numEntries() < opt.replacementThreshold  )|( datasetRVForFit[MHNominal].sumEntries() < 0. ):
  nominal_numEntries = datasetRVForFit[MHNominal].numEntries()
  procReplacementFit, catReplacementFit = rMap['procRVMap'][opt.cat], opt.cat #opt.cat // rMap['catRVMap'][opt.cat]
  for mp in opt.massPoints.split(","):
    WSFileName = findSignalWSFileNames(opt.inputWSDir,mp,procReplacementFit)[0]
    f = ROOT.TFile(WSFileName,"read")
    inputWS = f.Get(inputWSName__)
    if procReplacementFit.split("_")[-1] in ["i"]:
    # if procReplacementFit.split("_")[-1] in ["in", "out"]:
      d = inputWS.data("%s_%s_%s_%s_%s_%s"%(procToData(procReplacementFit.split("_")[0]),procToData(procReplacementFit.split("_")[-1]),mp,flavRVFit,sqrts__,catReplacementFit)).reduce(aset)
    else:
      d = loadReducedDataset(
        inputWS,getSignalDatasetName(procReplacementFit,mp,flavRVFit,catReplacementFit),
        inputFitArgSet(inputWS),WSFileName
      )
    fitDataset = selectFitDataset(d)
    datasetRVForFit[mp] = fitDataset
    inputWS.Delete()
    f.Close() 

  # Check if replacement dataset has too few entries: if so throw error
  if( datasetRVForFit[MHNominal].numEntries() < opt.replacementThreshold )|( datasetRVForFit[MHNominal].sumEntries() < 0. ):
    print(" --> [ERROR] replacement dataset (%s,%s) has too few entries (%g < %g)"%(procReplacementFit,catReplacementFit,datasetRVForFit[MHNominal].numEntries(),opt.replacementThreshold))
    sys.exit(1)

  else:
    procRVFit, catRVFit = procReplacementFit, catReplacementFit
    print(" --> Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat,flav) = (%s,%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procRVFit,catRVFit,flavRVFit))
    for mp in opt.massPoints.split(","):
      print("     * MH = %s GeV: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit[mp].numEntries(),datasetRVForFit[mp].sumEntries()))

else:
  print(" --> Using (proc,cat,flav) = (%s,%s,%s) for extracting shape"%(procRVFit,catRVFit,flavRVFit))
  for mp in opt.massPoints.split(","):
    print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit[mp].numEntries(),datasetRVForFit[mp].sumEntries()))

# Repeat for WV scenario

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If using nGaussian fit then extract nGaussians from fTest json file
if not opt.useDCB:
  with open("%s/outdir_%s/fTest/json/nGauss_%s.json"%(outputBase,opt.ext,catRVFit)) as jf: ngauss = json.load(jf)
  if( datasetRVForFit[MHNominal].numEntries() < opt.replacementThreshold  )|( datasetRVForFit[MHNominal].sumEntries() < 0. ):
    nRV = int(ngauss["%s_%s_%s"%(rMap['procRVMap'][catRVFit],catRVFit,flavRVFit)]['nRV'])
  else:
    nRV = int(ngauss["%s_%s_%s"%(procRVFit,catRVFit,flavRVFit)]['nRV'])
  print(" --> Fitting function: convolution of nGaussians (%g)"%nRV)
else:
  print(" --> Fitting function: %s"%("DCB + Gaussian" if opt.useDCBplusGaussian else "DCB"))

if opt.doVoigtian:
  print(" --> Will add natural Higgs width as parameter in Pdf (Gaussians -> Voigtians)")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT: simultaneous signal fit (ssf)
ssfMap = od()
name = "Total"
# FIXME: add flav to SimultaneousFit function
ssfRV = SimultaneousFit(name,opt.proc,opt.cat,opt.flav,datasetRVForFit,xvar.Clone(),MH,MHLow,MHHigh,opt.massPoints,opt.nBins,opt.MHPolyOrder,opt.minimizerMethod,opt.minimizerTolerance,opt.MHPolyConfig,verbose=False,_fitRange=[opt.fitMassLow,opt.fitMassHigh],_fitObjective=opt.fitObjective,_unbinnedWeightMode=opt.unbinnedWeightMode,_runMinuitPolish=opt.runMinuitPolish,_strictEdmThreshold=opt.strictEdmThreshold,_minuitPolishTolerance=opt.minuitPolishTolerance)
if opt.useDCB: 
  if opt.MHPolyConfig: 
    ssfRV.buildFromConfig(opt.MHPolyConfig)
  else: 
    ssfRV.buildDCBplusGaussian(addGaussian=opt.useDCBplusGaussian)
    if opt.fitSeedConfig: ssfRV.setFitParameterValues(opt.fitSeedConfig)
    ssfRV.runFit()
else: 
  ssfRV.buildNGaussians(nRV)
  ssfRV.runFit()
ssfRV.buildSplines()
ssfMap[name] = ssfRV


# Write optimizer evidence before final-model construction so a failed strict
# gate still leaves a machine-readable diagnostic for the rejected attempt.
diagnosticsDir = "%s/outdir_%s/signalFit/diagnostics"%(outputBase,opt.ext)
if not os.path.isdir(diagnosticsDir): os.makedirs(diagnosticsDir)
diagnosticsPath = "%s/optimizer_%s_%s_%s.json"%(diagnosticsDir,opt.proc,opt.cat,opt.flav)
optimizerPayload = {
  "schema_version": 1,
  "process": opt.proc,
  "source_processes": sourceProcsForFit,
  "category": opt.cat,
  "flavour": opt.flav,
  "objective": opt.fitObjective,
  "unbinned_weight_mode": opt.unbinnedWeightMode,
  "fit_partition_variable": opt.fitPartitionVar,
  "fit_partitions": fitPartitions if fitPartitions else "all",
  "strict_edm_threshold": opt.strictEdmThreshold,
  "require_strict_convergence": bool(opt.requireStrictConvergence),
  "fits": {fitName: fit.getOptimizerDiagnostics() for fitName,fit in ssfMap.items()},
  "strict_convergence_pass": all(
    fit.passesStrictConvergence() for fit in ssfMap.values()
  ) if opt.runMinuitPolish else None,
}
with open(diagnosticsPath,"w") as diagnosticsHandle:
  json.dump(optimizerPayload,diagnosticsHandle,indent=2,sort_keys=True)
  diagnosticsHandle.write("\n")
print(" --> Optimizer diagnostics: %s"%diagnosticsPath)
if opt.requireStrictConvergence and not optimizerPayload["strict_convergence_pass"]:
  print(" --> [ERROR] strict optimizer validation failed; see %s"%diagnosticsPath)
  sys.exit(2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FINAL MODEL: construction
print("\n --> Constructing final model")
# FIXME: add flav to FinalModel function
systematicsFile = opt.systematicsFile
if not opt.skipSystematics and systematicsFile is None:
  relative = os.path.join("calcPhotonSyst", "pkl", opt.cat + ".pkl")
  paths = [os.path.join(outputBase, "outdir_" + opt.ext, relative),
           os.path.join(swd__, "outdir_" + opt.ext, relative),
           os.path.join(swd__, "output_" + opt.ext, relative)]
  systematicsFile = next((path for path in paths if os.path.isfile(path)), paths[0])
fm = FinalModel(ssfMap,opt.proc,opt.cat,opt.flav,opt.ext,opt.year,sqrts__,nominalDatasets,xvar,MH,MHLow,MHHigh,opt.massPoints,xsbrMap,procSyst,opt.scales,opt.scalesCorr,opt.scalesGlobal,opt.smears,opt.doVoigtian,opt.useDCB,opt.skipSystematics,opt.doEffAccFromJson,systematicsFile)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE: to output workspace
foutDir = "%s/outdir_%s/signalFit/output"%(outputBase,opt.ext)
foutName = "%s/CMS-HZG_sigfit_%s_%s_%s_%s_%s.root"%(foutDir,opt.ext,opt.proc,opt.year,opt.flav,opt.cat)
print("\n --> Saving output workspace to file: %s"%foutName)
os.makedirs(foutDir, exist_ok=True)
fout = ROOT.TFile(foutName,"RECREATE")
outWS = ROOT.RooWorkspace("%s_%s"%(outputWSName__,sqrts__),"%s_%s"%(outputWSName__,sqrts__))
# The interpolation scan leaves MH at its last grid point.
MH.setVal(float(MHNominal))
fm.save(outWS)
outWS.Write()
fout.Close()
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTTING
if opt.doPlots:
  print("\n --> Making plots...")
  os.makedirs("%s/outdir_%s/signalFit/Plots"%(outputBase,opt.ext), exist_ok=True)
  # FIXME: add flav to plotPdfComponents function
  plotPdfComponents(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots"%(outputBase,opt.ext),_extension="total_",_proc=procRVFit,_cat=catRVFit,_flav=opt.flav) 
  # Plot interpolation
  # FIXME: add flav to plotInterpolation and plotSplines functions
  try:
    plotInterpolation(fm,_outdir="%s/outdir_%s/signalFit/Plots"%(outputBase,opt.ext)) 
    plotSplines(fm,_outdir="%s/outdir_%s/signalFit/Plots"%(outputBase,opt.ext),_nominalMass=MHNominal) 
  except:
    pass

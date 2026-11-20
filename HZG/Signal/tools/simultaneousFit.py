# Class for performing simultaneous signal fit
import ROOT
import json
import numpy as np
from scipy.optimize import minimize
import scipy.stats
import ctypes
from collections import OrderedDict as od
from array import array
import os
import math

# Parameter lookup table for initialisation
# So far defined up to MHPolyOrder=2
pLUT = od()
pLUT['DCB'] = od()
pLUT['DCB']['dm_p0'] = [0.1,-1.,1.]
pLUT['DCB']['dm_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['dm_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['sigma_p0'] = [1.,0.,5.]
pLUT['DCB']['sigma_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['sigma_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['n1_p0'] = [20.,1.00001,500]
pLUT['DCB']['n1_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['n1_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['n2_p0'] = [20.,1.00001,500]
pLUT['DCB']['n2_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['n2_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['a1_p0'] = [1.,1.,10.]
pLUT['DCB']['a1_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['a1_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['a2_p0'] = [1.,1.,20.]
pLUT['DCB']['a2_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['a2_p2'] = [0.0,-0.001,0.001]
pLUT['Gaussian_wdcb'] = od()
pLUT['Gaussian_wdcb']['dm_p0'] = [0.1,-1.5,1.5]
pLUT['Gaussian_wdcb']['dm_p1'] = [0.01,-0.01,0.01]
pLUT['Gaussian_wdcb']['dm_p2'] = [0.01,-0.01,0.01]
pLUT['Gaussian_wdcb']['sigma_p0'] = [1.5,1.0,4.]
pLUT['Gaussian_wdcb']['sigma_p1'] = [0.0,-0.1,0.1]
pLUT['Gaussian_wdcb']['sigma_p2'] = [0.0,-0.001,0.001]
pLUT['Frac'] = od()
pLUT['Frac']['p0'] = [0.25,0.01,0.99]
pLUT['Frac']['p1'] = [0.,-0.05,0.05]
pLUT['Frac']['p2'] = [0.,-0.0001,0.0001]
pLUT['Gaussian'] = od()
pLUT['Gaussian']['dm_p0'] = [0.1,-5.,5.]
pLUT['Gaussian']['dm_p1'] = [0.0,-0.01,0.01]
pLUT['Gaussian']['dm_p2'] = [0.0,-0.01,0.01]
pLUT['Gaussian']['sigma_p0'] = ['func',0.5,10.0]
pLUT['Gaussian']['sigma_p1'] = [0.0,-0.01,0.01]
pLUT['Gaussian']['sigma_p2'] = [0.0,-0.01,0.01]
pLUT['FracGaussian'] = od()
pLUT['FracGaussian']['p0'] = ['func',0.01,0.99]
pLUT['FracGaussian']['p1'] = [0.01,-0.005,0.005]
pLUT['FracGaussian']['p2'] = [0.00001,-0.00001,0.00001]

if os.environ.get("HZG_DSCB_RELAX_TAIL_BOUNDS", ""):
  a_min = float(os.environ.get("HZG_DSCB_A_MIN", "0.5"))
  n_max = float(os.environ.get("HZG_DSCB_N_MAX", "2000"))
  for key in ["a1_p0", "a2_p0"]:
    pLUT["DCB"][key][1] = a_min
  for key in ["n1_p0", "n2_p0"]:
    pLUT["DCB"][key][2] = n_max

# Function to convert sumw2 variance to poisson interval
def poisson_interval(x,eSumW2,level=0.68):
  neff = x**2/(eSumW2**2)
  scale = abs(x)/neff
  l = scipy.stats.gamma.interval(level, neff, scale=scale,)[0]
  u = scipy.stats.gamma.interval(level, neff+1, scale=scale,)[1]
  # protect against no effective entries
  l[neff==0] = 0.
  # protect against no variance
  l[eSumW2==0.] = 0.
  u[eSumW2==0.] = np.inf
  # convert to upper and lower errors
  eLo, eHi = abs(l-x),abs(u-x)
  return eLo, eHi

# Function to calc chi2 for binned fit given pdf, RooDataHist and xvar as inputs
#def calcChi2(x,pdf,d,errorType="Sumw2",_verbose=False,fitRange=[100,180]):
#def calcChi2(x,pdf,d,errorType="Poisson",_verbose=False,fitRange=[110,140]):
def calcChi2(x,pdf,d,errorType="Poisson",_verbose=False,fitRange=[105,150]):

  k = 0. # number of non empty bins (for calc degrees of freedom)
  normFactor = d.sumEntries()
  normalization = ROOT.RooArgSet(x)
  
  # Using numpy and poisson error
  bins, nPdf, nData, eDataSumW2 = [], [],[],[]
  for i in range(d.numEntries()):
    p = d.get(i)
    x.setVal(p.getRealValue(x.GetName()))
    if( x.getVal() < fitRange[0] )|( x.getVal() > fitRange[1] ): continue
    ndata = d.weight()
    if ndata*ndata == 0: continue
    npdf = pdf.getVal(normalization)*normFactor*d.binVolume()
    eLo, eHi = ctypes.c_double(0.), ctypes.c_double(0.)
    d.weightError(eLo,eHi,ROOT.RooAbsData.SumW2)
    bins.append(i)
    nPdf.append(npdf)
    nData.append(ndata)
    eDataSumW2.append(eHi.value) if npdf>ndata else eDataSumW2.append(eLo.value)
    k += 1

  # Convert to numpy array
  nPdf = np.asarray(nPdf)
  nData = np.asarray(nData)
  eDataSumW2 = np.asarray(eDataSumW2)
  if (len(nPdf) == 0) or (not np.all(np.isfinite(nPdf))):
    return 1.e30,k

  if errorType == 'Poisson':
    # Change error to poisson intervals: take max interval as error
    eLo,eHi = poisson_interval(nData,eDataSumW2,level=0.68)
    #eDataPoisson = 0.5*(eHi+eLo)
    eDataPoisson = np.maximum(eHi,eLo) 
    #eDataPoisson = (nPdf>nData)*eHi + (nPdf<=nData)*eLo 
    e = eDataPoisson
    # Calculate chi2 terms
    terms = (nPdf-nData)**2/(eDataPoisson**2)
  elif errorType == "Expected":
    # Change error to sqrt pdf entries
    eExpected = np.sqrt(nPdf)
    e = eExpected
    # Calculate chi2 terms
    terms = (nPdf-nData)**2/(eExpected**2)
  else:
    # Use SumW2 terms to calculate chi2
    e = eDataSumW2
    terms = (nPdf-nData)**2/(eDataSumW2**2)
   
  if not np.all(np.isfinite(terms)):
    return 1.e30,k

  # If verbose: print to screen
  if _verbose:
    for i in range(len(terms)):
      print(" --> [DEBUG] Bin %g : nPdf = %.6f, nData = %.6f, e(%s) = %.6f --> chi2 term = %.6f"%(bins[i],nPdf[i],nData[i],errorType,e[i],terms[i]))

  # Sum terms
  result = terms.sum()

  return result,k

# Function to add chi2 for multiple mass points
def nChi2Addition(X,ssf,verbose=False):
  # X: vector of param values (updated with minimise function)
  # Loop over parameters and set RooVars
  for i in range(len(X)): ssf.FitParameters[i].setVal(X[i])
  # Loop over datasets: adding chi2 for each mass point
  chi2sum = 0
  K = 0 # number of non empty bins
  C = sum(1 for p in ssf.FitParameters if p.GetName() != ssf.MH.GetName())
  for mp,d in ssf.DataHists.items():
    ssf.MH.setVal(int(mp))
    chi2, k  = calcChi2(ssf.xvar,ssf.Pdfs['final'],d,_verbose=verbose,fitRange=ssf.fitRange)
    chi2sum += chi2
    K += k
  # N degrees of freedom
  ndof = K-C
  ssf.setNdof(ndof)
  return chi2sum

# Function to add weighted unbinned NLL values for multiple mass points.
# The input datasets are normalized per mass point in prepareDataHists, matching
# the equal-mass-point treatment used by the existing binned chi2 objective.
def nNLLAddition(X,ssf):
  for i in range(len(X)): ssf.FitParameters[i].setVal(X[i])
  nllsum = 0.
  for mp,nll in ssf.NLLs.items():
    ssf.MH.setVal(int(mp))
    # createNLL owns a cloned computation graph. Keep its coefficient
    # variables synchronized with the original variables controlled by scipy.
    for i in range(len(X)):
      target = ssf.NLLParameterMaps[mp].get(ssf.FitParameters[i].GetName())
      if target: target.setVal(X[i])
    nllMH = ssf.NLLParameterMaps[mp].get(ssf.MH.GetName())
    if nllMH: nllMH.setVal(int(mp))
    # Polynomial coefficients can each be in range while their value at a
    # non-central mass point leaves the physical RooDoubleCBFast domain.
    # Reject those trial points before RooFit evaluates every event.
    dcbPhysicalMins = {
      'sigma_dcb': 0.,
      'a1_dcb': 0.,
      'a2_dcb': 0.,
      'n1_dcb': 1.,
      'n2_dcb': 1.
    }
    for name,minimum in dcbPhysicalMins.items():
      if name in ssf.Polynomials:
        value = ssf.Polynomials[name].getVal()
        if (not np.isfinite(value)) or (value <= minimum): return 1.e30
    value = nll.getVal()
    if not np.isfinite(value): return 1.e30
    nllsum += ssf.NLLScales[mp]*value
  return nllsum

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
class SimultaneousFit:
  # Constructor
  def __init__(self,_name,_proc,_cat,_flav,_datasetForFit,_xvar,_MH,_MHLow,_MHHigh,_massPoints,_nBins,_MHPolyOrder,_minimizerMethod,_minimizerTolerance,_MHPolyConfig=None,verbose=False,_fitRange=None,_fitObjective="chi2",_unbinnedWeightMode="signed",_runMinuitPolish=False,_strictEdmThreshold=1.e-6,_minuitPolishTolerance=None):
    self.name = _name
    self.proc = _proc
    self.cat = _cat
    self.flav = _flav
    self.datasetForFit = _datasetForFit
    self.xvar = _xvar
    self.MH = _MH
    self.MHLow = _MHLow
    self.MHHigh = _MHHigh
    self.massPoints = _massPoints
    self.nBins = _nBins
    self.MHPolyOrder = _MHPolyOrder
    self.minimizerMethod = _minimizerMethod
    self.minimizerTolerance = _minimizerTolerance
    self.MHPolyConfig = _MHPolyConfig
    self.verbose = verbose
    self.fitRange = list(_fitRange) if _fitRange else [105.,150.]
    self.fitObjective = _fitObjective
    self.unbinnedWeightMode = _unbinnedWeightMode
    self.runMinuitPolish = bool(_runMinuitPolish)
    self.strictEdmThreshold = float(_strictEdmThreshold)
    self.minuitPolishTolerance = (max(1.e-6,self.strictEdmThreshold*1000.)
                                  if _minuitPolishTolerance is None else float(_minuitPolishTolerance))
    if not math.isfinite(self.minuitPolishTolerance) or self.minuitPolishTolerance <= 0.:
      raise ValueError("Minuit polish tolerance must be finite and positive")
    # Prepare vars
    self.MH.setConstant(False)
    self.MH.setVal(125)
    self.MH.setBins(10)
    self.dMH = ROOT.RooFormulaVar("dMH","dMH","@0-125.0",ROOT.RooArgList(self.MH)) 
    self.xvar.setVal(125)
    self.xvar.setBins(self.nBins)
    # Dicts to store all fit vars, polynomials, pdfs and splines
    self.nGaussians = 1
    self.Vars = od()
    self.Varlists = od()
    self.Polynomials = od()
    self.Pdfs = od()
    self.Coeffs = od()
    self.Splines = od()
    # Prepare normalized RooDataSets and RooDataHists for fit
    self.DataSets = od()
    self.DataHists = od()
    self.NLLScales = od()
    self.prepareDataHists()
    # Fit containers
    self.FitParameters = None
    self.Ndof = None
    self.Chi2 = None
    self.NLL = None
    self.NLLs = od()
    self.NLLParameterMaps = od()
    self.FitResult = None
    self.OptimizerDiagnostics = od()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function for setting N degrees of freedom
  def setNdof(self,_ndof): self.Ndof = _ndof  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to extract param bounds
  def extractXBounds(self):
    XBounds = []
    for i in range(len(self.FitParameters)): XBounds.append((self.FitParameters[i].getMin(),self.FitParameters[i].getMax()))
    return np.asarray(XBounds)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to extract initial param value vector
  def extractX0(self):
    X0 = []
    for i in range(len(self.FitParameters)): X0.append(self.FitParameters[i].getVal())
    return np.asarray(X0)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to normalise datasets and convert to RooDataHists for calc chi2
  def prepareDataHists(self):
    # Loop over datasets and normalise to 1
    for k,d in self.datasetForFit.items():
      sumw = d.sumEntries()
      if self.fitObjective == "unbinnedNLL":
        fitSumW = 0.
        fitSumW2 = 0.
        for i in range(0,d.numEntries()):
          d.get(i)
          fitWeight = d.weight()
          if self.unbinnedWeightMode == "positive": fitWeight = max(fitWeight,0.)
          elif self.unbinnedWeightMode == "absolute": fitWeight = abs(fitWeight)
          fitSumW += fitWeight
          fitSumW2 += fitWeight*fitWeight
        if fitSumW <= 0. or fitSumW2 <= 0.:
          raise ValueError("Unbinned fit dataset %s has non-positive usable weight"%str(k))
        self.NLLScales[k] = (fitSumW*fitSumW/fitSumW2) if fitSumW2 > 0. else 1.
      drw = d.emptyClone()
      if self.fitObjective == "unbinnedNLL": dnll = d.emptyClone()
      self.Vars['weight'] = ROOT.RooRealVar("weight","weight",-10000,10000)
      for i in range(0,d.numEntries()):
        self.xvar.setVal(d.get(i).getRealValue(self.xvar.GetName()))
        self.Vars['weight'].setVal((1/sumw)*d.weight())
        drw.add(ROOT.RooArgSet(self.xvar,self.Vars['weight']),self.Vars['weight'].getVal())
        if self.fitObjective == "unbinnedNLL":
          fitWeight = d.weight()
          if self.unbinnedWeightMode == "positive": fitWeight = max(fitWeight,0.)
          elif self.unbinnedWeightMode == "absolute": fitWeight = abs(fitWeight)
          if fitWeight != 0.:
            self.Vars['weight'].setVal((1/fitSumW)*fitWeight)
            dnll.add(ROOT.RooArgSet(self.xvar,self.Vars['weight']),self.Vars['weight'].getVal())
      if self.fitObjective == "unbinnedNLL": self.DataSets[k] = dnll
      # Convert to RooDataHist
      self.DataHists[k] = ROOT.RooDataHist("%s_hist"%d.GetName(),"%s_hist"%d.GetName(),ROOT.RooArgSet(self.xvar),drw)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Build native RooFit NLL objects once so repeated scipy evaluations remain
  # in C++. MH is deliberately non-constant and is set for each mass point.
  def prepareNLLs(self):
    self.xvar.setRange("signalFitRange",self.fitRange[0],self.fitRange[1])
    self.NLLs = od()
    for mp,d in self.DataSets.items():
      self.MH.setVal(int(mp))
      self.NLLs[mp] = self.Pdfs['final'].createNLL(
        d,
        ROOT.RooFit.Range("signalFitRange")
      )
      nllVars = self.NLLs[mp].getVariables()
      self.NLLParameterMaps[mp] = {
        p.GetName(): nllVars.find(p.GetName()) for p in self.FitParameters
      }
      if os.environ.get("HZG_SIGNAL_FIT_DIAGNOSTICS", ""):
        missing = [name for name,target in self.NLLParameterMaps[mp].items() if not target]
        print(" --> NLL_PARAMETER_MAP mass=%s mapped=%d missing=%s"%(
          str(mp), len(self.NLLParameterMaps[mp])-len(missing), ",".join(missing) if missing else "none"
        ))
  
  # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # def buildDCBplusGaussian(self,_recursive=True):

  #   # DCB
  #   # Define polynominal functions (in dMH)
  #   for f in ['dm','sigma','n1','n2','a1','a2']: 
  #     k = "%s_dcb"%f
  #     self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
  #     # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
  #     for po in range(0,self.MHPolyOrder+1):
  #       self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['DCB']["%s_p%s"%(f,po)][0],pLUT['DCB']["%s_p%s"%(f,po)][1],pLUT['DCB']["%s_p%s"%(f,po)][2])
	# self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] ) 
  #     # Define polynominal
  #     self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
  #   # Mean function
  #   self.Polynomials['mean_dcb'] = ROOT.RooFormulaVar("mean_dcb","mean_dcb","(@0+@1)",ROOT.RooArgList(self.MH,self.Polynomials['dm_dcb']))
  #   # Build DCB
  #   self.Pdfs['dcb'] = ROOT.RooDoubleCBFast("dcb","dcb",self.xvar,self.Polynomials['mean_dcb'],self.Polynomials['sigma_dcb'],self.Polynomials['a1_dcb'],self.Polynomials['n1_dcb'],self.Polynomials['a2_dcb'],self.Polynomials['n2_dcb'])

  #   # Gaussian
  #   # Define polynomial function for sigma (in dMH). Gaussian defined to have same mean as DCB
  #   f = "sigma"
  #   k = "%s_gaus"%f
  #   self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)    
  #   # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
  #   for po in range(0,self.MHPolyOrder+1):
  #     self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['Gaussian_wdcb']["%s_p%s"%(f,po)][0],pLUT['Gaussian_wdcb']["%s_p%s"%(f,po)][1],pLUT['Gaussian_wdcb']["%s_p%s"%(f,po)][2])
  #     self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] )
  #   # Define polynomial
  #   self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
  #   # Build Gaussian
  #   self.Pdfs['gaus'] = ROOT.RooGaussian("gaus","gaus",self.xvar,self.Polynomials['mean_dcb'],self.Polynomials['sigma_gaus'])
    
        
  #   # Relative fraction: also polynomial of order MHPolyOrder
  #   self.Varlists['frac'] = ROOT.RooArgList("frac_coeffs")
  #   for po in range(0,self.MHPolyOrder+1):
  #     self.Vars['frac_p%g'%po] = ROOT.RooRealVar("frac_p%g"%po,"frac_p%g"%po,pLUT['Frac']['p%g'%po][0],pLUT['Frac']['p%g'%po][1],pLUT['Frac']['p%g'%po][2])
  #     self.Varlists['frac'].add( self.Vars['frac_p%g'%po] )
  #   # Define Polynomial
  #   self.Polynomials['frac'] = ROOT.RooPolyVar('frac','frac',self.dMH,self.Varlists['frac'])
  #   # Constrain fraction to not be above 1 or below 0
  #   self.Polynomials['frac_constrained'] = ROOT.RooFormulaVar("frac_constrained","frac_constrained","(@0>0)*(@0<1)*@0+(@0>1.0)*0.9999",ROOT.RooArgList(self.Polynomials['frac']))
  #   self.Coeffs['frac_constrained'] = self.Polynomials['frac_constrained' ]

  #   # Define total PDF
  #   _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
  #   for pdf in ['dcb','gaus']: _pdfs.add(self.Pdfs[pdf])
  #   _coeffs.add(self.Coeffs['frac_constrained'])
  #   self.Pdfs['final'] = ROOT.RooAddPdf("%s_%s"%(self.proc,self.cat),"%s_%s"%(self.proc,self.cat),_pdfs,_coeffs,_recursive)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  def buildFromConfig(self, config_file):
    if self.verbose: 
        print(" --> (%s) Building DCB from config: %s" % (self.name, config_file))
    
    # Read from json file
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    # DCB
    for f in ['dm','sigma','n1','n2','a1','a2']: 
            
        k = "%s_dcb" % f
        self.Varlists[k] = ROOT.RooArgList("%s_coeffs" % k)
        
        # Read parameter values from config
        for po in range(0, self.MHPolyOrder+1):
            param_name = '%s_p%s' % (f, po)
            if param_name in config:
                param_value = config[param_name]
                # Set to constant(min=max=fix value)
                self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar(
                    "%s_p%g"%(k,po), "%s_p%g"%(k,po),
                    param_value, param_value, param_value  # Constant
                )
            else:
                # Not found parameter in config file
                default_val = pLUT['DCB']["%s_p%s"%(f,po)][0]
                self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar(
                    "%s_p%g"%(k,po), "%s_p%g"%(k,po),
                    default_val, default_val, default_val  # Constant
                )
            
            self.Varlists[k].add(self.Vars['%s_p%g'%(k,po)])
        
        # Polynomial
        self.Polynomials[k] = ROOT.RooPolyVar(k, k, self.dMH, self.Varlists[k])
    
    # Mean function
    self.Polynomials['mean_dcb'] = ROOT.RooFormulaVar(
        "mean_dcb", "mean_dcb", "(@0+@1)",
        ROOT.RooArgList(self.MH, self.Polynomials['dm_dcb'])
    )
    
    # DCB PDF
    self.Pdfs['dcb'] = ROOT.RooDoubleCBFast(
        "dcb", "dcb", self.xvar, self.Polynomials['mean_dcb'],
        self.Polynomials['sigma_dcb'], self.Polynomials['a1_dcb'],
        self.Polynomials['n1_dcb'], self.Polynomials['a2_dcb'], 
        self.Polynomials['n2_dcb']
    )
    
    self.Pdfs['final'] = self.Pdfs['dcb']
    fv = self.Pdfs['final'].getVariables().Clone()
    fv.remove(self.xvar)
    self.FitParameters = ROOT.RooArgList(fv)
    self.Chi2 = self.getChi2()
    
    if self.verbose:
        print(" --> (%s) DCB built from config, all parameters fixed" % self.name)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  def buildDCBplusGaussian(self,_recursive=True,addGaussian=False):

    # DCB
    # Define polynominal functions (in dMH)
    for f in ['dm','sigma','n1','n2','a1','a2']: 
      k = "%s_dcb"%f
      self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
      # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
      for po in range(0,self.MHPolyOrder+1):
        self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['DCB']["%s_p%s"%(f,po)][0],pLUT['DCB']["%s_p%s"%(f,po)][1],pLUT['DCB']["%s_p%s"%(f,po)][2])
        self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] ) 
      # Define polynominal
      self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
    # Mean function
    self.Polynomials['mean_dcb'] = ROOT.RooFormulaVar("mean_dcb","mean_dcb","(@0+@1)",ROOT.RooArgList(self.MH,self.Polynomials['dm_dcb']))
    # Build DCB
    self.Pdfs['dcb'] = ROOT.RooDoubleCBFast("dcb","dcb",self.xvar,self.Polynomials['mean_dcb'],self.Polynomials['sigma_dcb'],self.Polynomials['a1_dcb'],self.Polynomials['n1_dcb'],self.Polynomials['a2_dcb'],self.Polynomials['n2_dcb'])
    
    if not addGaussian:
      self.Pdfs['final'] = self.Pdfs['dcb']
      return

    # Optional broad Gaussian shoulder sharing the DSCB mean.  This keeps the
    # pure-DSCB path unchanged while allowing categories whose core and
    # shoulders cannot be described by a single DSCB to be tuned explicitly.
    k = "sigma_gaus"
    self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
    for po in range(0,self.MHPolyOrder+1):
      self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar(
        "%s_p%g"%(k,po), "%s_p%g"%(k,po),
        pLUT['Gaussian_wdcb']["sigma_p%s"%po][0],
        pLUT['Gaussian_wdcb']["sigma_p%s"%po][1],
        pLUT['Gaussian_wdcb']["sigma_p%s"%po][2]
      )
      self.Varlists[k].add(self.Vars['%s_p%g'%(k,po)])
    self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
    self.Pdfs['gaus'] = ROOT.RooGaussian(
      "gaus","gaus",self.xvar,self.Polynomials['mean_dcb'],self.Polynomials[k]
    )

    self.Varlists['frac'] = ROOT.RooArgList("frac_coeffs")
    for po in range(0,self.MHPolyOrder+1):
      self.Vars['frac_p%g'%po] = ROOT.RooRealVar(
        "frac_p%g"%po,"frac_p%g"%po,
        pLUT['Frac']['p%g'%po][0],pLUT['Frac']['p%g'%po][1],pLUT['Frac']['p%g'%po][2]
      )
      self.Varlists['frac'].add(self.Vars['frac_p%g'%po])
    self.Polynomials['frac'] = ROOT.RooPolyVar('frac','frac',self.dMH,self.Varlists['frac'])
    self.Polynomials['frac_constrained'] = ROOT.RooFormulaVar(
      "frac_constrained","frac_constrained",
      "(@0>0)*(@0<1)*@0+(@0>=1)*0.9999",
      ROOT.RooArgList(self.Polynomials['frac'])
    )
    self.Coeffs['frac_constrained'] = self.Polynomials['frac_constrained']

    _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
    _pdfs.add(self.Pdfs['dcb'])
    _pdfs.add(self.Pdfs['gaus'])
    _coeffs.add(self.Coeffs['frac_constrained'])
    self.Pdfs['final'] = ROOT.RooAddPdf(
      "%s_%s_%s"%(self.proc,self.cat,self.flav),
      "%s_%s_%s"%(self.proc,self.cat,self.flav),
      _pdfs,_coeffs,_recursive
    )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  def buildNGaussians(self,nGaussians,_recursive=True):

    # Set number of gaussians
    self.nGaussians = nGaussians

    # Loop over NGaussians
    for g in range(0,nGaussians):
      # Define polynominal functions for mean and sigma (in MH)
      for f in ['dm','sigma']: 
        k = "%s_g%g"%(f,g)
        self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
        # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
        for po in range(0,self.MHPolyOrder+1):
          # p0 value of sigma is function of g (creates gaussians of increasing width)
          if(f == "sigma")&(po==0): 
            self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),(g+1)*1.0,pLUT['Gaussian']["%s_p%s"%(f,po)][1],pLUT['Gaussian']["%s_p%s"%(f,po)][2])
          else:
            self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['Gaussian']["%s_p%s"%(f,po)][0],pLUT['Gaussian']["%s_p%s"%(f,po)][1],pLUT['Gaussian']["%s_p%s"%(f,po)][2])
          self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] ) 
        # Define polynominal
        self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
      # Mean function
      self.Polynomials['mean_g%g'%g] = ROOT.RooFormulaVar("mean_g%g"%g,"mean_g%g"%g,"(@0+@1)",ROOT.RooArgList(self.MH,self.Polynomials['dm_g%g'%g]))
      # Build Gaussian
      self.Pdfs['gaus_g%g'%g] = ROOT.RooGaussian("gaus_g%g"%g,"gaus_g%g"%g,self.xvar,self.Polynomials['mean_g%g'%g],self.Polynomials['sigma_g%g'%g])

      # Relative fractions: also polynomials of order MHPolyOrder (define up to n=nGaussians-1)
      if g < nGaussians-1:
        self.Varlists['frac_g%g'%g] = ROOT.RooArgList("frac_g%g_coeffs"%g)
        for po in range(0,self.MHPolyOrder+1):
          if po == 0:
            self.Vars['frac_g%g_p%g'%(g,po)] = ROOT.RooRealVar("frac_g%g_p%g"%(g,po),"frac_g%g_p%g"%(g,po),0.5-0.05*g,pLUT['FracGaussian']['p%g'%po][1],pLUT['FracGaussian']['p%g'%po][2])
          else:
            self.Vars['frac_g%g_p%g'%(g,po)] = ROOT.RooRealVar("frac_g%g_p%g"%(g,po),"frac_g%g_p%g"%(g,po),pLUT['FracGaussian']['p%g'%po][0],pLUT['FracGaussian']['p%g'%po][1],pLUT['FracGaussian']['p%g'%po][2])
          self.Varlists['frac_g%g'%g].add( self.Vars['frac_g%g_p%g'%(g,po)] )
        # Define Polynomial
        self.Polynomials['frac_g%g'%g] = ROOT.RooPolyVar("frac_g%g"%g,"frac_g%g"%g,self.dMH,self.Varlists['frac_g%g'%g])
        # Constrain fraction to not be above 1 or below 0
        self.Polynomials['frac_g%g_constrained'%g] = ROOT.RooFormulaVar('frac_g%g_constrained'%g,'frac_g%g_constrained'%g,"(@0>0)*(@0<1)*@0+ (@0>1.0)*0.9999",ROOT.RooArgList(self.Polynomials['frac_g%g'%g]))
        self.Coeffs['frac_g%g_constrained'%g] = self.Polynomials['frac_g%g_constrained'%g]
    # End of loop over n Gaussians
    
    # Define total PDF
    _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
    for g in range(0,nGaussians): 
      _pdfs.add(self.Pdfs['gaus_g%g'%g])
      if g < nGaussians-1: _coeffs.add(self.Coeffs['frac_g%g_constrained'%g])
    # Use the standard RooAddPdf coefficient convention for the Gaussian-like
    # model. In current PyROOT/ROOT this avoids the recursive-fraction overload
    # ambiguity observed for fixed nGaussian fits.
    self.Pdfs['final'] = ROOT.RooAddPdf("%s_%s_%s"%(self.proc,self.cat,self.flav),"%s_%s_%s"%(self.proc,self.cat,self.flav),_pdfs,_coeffs)
    
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  def runFit(self):
    # Extract fit variables: remove xvar from fit parameters
    fv = self.Pdfs['final'].getVariables().Clone()
    fv.remove(self.xvar)
    # MH is set explicitly for each mass-point contribution. Leaving it in the
    # vector gives the optimizer an exactly flat direction and prevents a
    # meaningful covariance/EDM check. Preserve legacy behavior unless the
    # scaled optimizer or the opt-in strict validation path is requested.
    if self.minimizerMethod == "L-BFGS-B" or self.runMinuitPolish:
      fv.remove(self.MH)
    # A seed configuration may deliberately fix a weakly identified DSCB
    # coefficient.  Keep that coefficient in the PDF, but do not give scipy or
    # Minuit2 a zero-dimensional/constant direction.
    constantVars = [var for var in fv if var.isConstant()]
    for var in constantVars:
      fv.remove(var)
    self.FitParameters = ROOT.RooArgList(fv)
    
    # Create initial vector of parameters and calculate initial Chi2
    if self.verbose: print("\n --> (%s) Initialising fit parameters"%self.name)
    x0 = self.extractX0()
    xbounds = self.extractXBounds()
    self.Chi2 = self.getChi2()
    if self.fitObjective == "unbinnedNLL":
      self.prepareNLLs()
      self.NLL = self.getNLL()
      probe = x0.copy()
      probeIndex = next((i for i,p in enumerate(self.FitParameters) if p.GetName() != self.MH.GetName()),None)
      if probeIndex is not None:
        span = xbounds[probeIndex][1]-xbounds[probeIndex][0]
        probe[probeIndex] = min(xbounds[probeIndex][1],probe[probeIndex]+max(1.e-5,1.e-4*span))
        probeNLL = nNLLAddition(probe,self)
        nNLLAddition(x0,self)
        if os.environ.get("HZG_SIGNAL_FIT_DIAGNOSTICS", ""):
          print(" --> NLL_SENSITIVITY parameter=%s initial=%.9g probe=%.9g delta=%.9g"%(
            self.FitParameters[probeIndex].GetName(),self.NLL,probeNLL,probeNLL-self.NLL
          ))
        if abs(probeNLL-self.NLL) <= 1.e-10*max(1.,abs(self.NLL)):
          raise RuntimeError("Unbinned NLL does not respond to parameter %s"%self.FitParameters[probeIndex].GetName())
    # Print parameter pre-fit values
    if self.verbose: self.printFitParameters(title="Pre-fit")

    # Run fit
    if self.verbose: print(" --> (%s) Running fit"%self.name)
    objective = nNLLAddition if self.fitObjective == "unbinnedNLL" else nChi2Addition
    if self.minimizerMethod == "L-BFGS-B":
      # Optimize dimensionless coordinates in [0,1]. The DSCB coefficients
      # span very different physical scales, so raw finite differences can be
      # numerically zero for some parameters while overshooting others.
      lower = xbounds[:,0]
      upper = xbounds[:,1]
      spans = upper-lower
      if np.any(spans <= 0.):
        raise RuntimeError("Unbinned fit has a non-positive parameter span")
      u0 = np.clip((x0-lower)/spans,0.,1.)
      ubounds = [(0.,1.) for _ in range(len(u0))]
      def scaledObjective(U):
        return objective(lower+np.asarray(U)*spans,self)
      result = minimize(
        scaledObjective,u0,bounds=ubounds,method=self.minimizerMethod,
        tol=self.minimizerTolerance,
        options={
          'ftol': max(self.minimizerTolerance,1.e-12),
          'gtol': 1.e-6,
          'eps': 1.e-5,
          'maxiter': 2000,
          'maxls': 50
        }
      )
      result['x_scaled'] = np.asarray(result['x']).copy()
      result['x'] = lower+np.asarray(result['x'])*spans
      self.FitResult = result
    else:
      self.FitResult = minimize(
        objective,x0,args=self,bounds=xbounds,method=self.minimizerMethod,
        tol=self.minimizerTolerance
      )
    if np.all(np.isfinite(self.FitResult['x'])):
      for i in range(len(self.FitResult['x'])):
        self.FitParameters[i].setVal(self.FitResult['x'][i])
    self._recordScipyDiagnostics()
    if self.runMinuitPolish:
      self._runMinuitPolish(objective,xbounds)
    self.Chi2 = self.getChi2()
    if self.fitObjective == "unbinnedNLL": self.NLL = self.getNLL()
    #self.Chi2 = nChi2Addition(self.FitResult['x'],self)
    # Print parameter post-fit values
    if self.verbose: self.printFitParameters(title="Post-fit")
    if os.environ.get("HZG_SIGNAL_FIT_DIAGNOSTICS", ""):
      status = "NA"
      message = "NA"
      if self.FitResult is not None:
        status = str(self.FitResult.get("status", "NA"))
        message = str(self.FitResult.get("message", "NA"))
      reduced = self.Chi2/int(self.Ndof) if self.Ndof else float("nan")
      objectiveValue = self.NLL if self.fitObjective == "unbinnedNLL" else self.Chi2
      print(" --> FIT_DIAGNOSTIC name=%s objective=%s objective_value=%.6f chi2=%.6f ndof=%s reduced_chi2=%.6f status=%s message=%s"%(
        self.name, self.fitObjective, objectiveValue, self.Chi2, str(self.Ndof), reduced, status, message
      ))
      if self.FitResult is not None:
        print(" --> FIT_OPTIMIZER nit=%s nfev=%s success=%s"%(
          str(self.FitResult.get('nit','NA')),
          str(self.FitResult.get('nfev','NA')),
          str(self.FitResult.get('success','NA'))
        ))
      print(" --> FIT_PARAMETERS %s"%" ".join(
        "%s=%.9g"%(p.GetName(),p.getVal()) for p in self.FitParameters
        if p.GetName() != self.MH.GetName()
      ))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Persist optimizer evidence in a framework-neutral, JSON-serializable form.
  def _recordScipyDiagnostics(self):
    result = self.FitResult
    jac = np.asarray(result.get("jac",[]),dtype=float)
    finite_jac = jac.size > 0 and np.all(np.isfinite(jac))
    self.OptimizerDiagnostics["scipy"] = {
      "method": self.minimizerMethod,
      "tolerance": float(self.minimizerTolerance),
      "success": bool(result.get("success",False)),
      "status": int(result.get("status",-1)),
      "message": str(result.get("message","")),
      "objective": float(result.get("fun",float("nan"))),
      "iterations": int(result.get("nit",-1)),
      "function_evaluations": int(result.get("nfev",-1)),
      "gradient_norm": float(np.linalg.norm(jac)) if finite_jac else None,
      "gradient_max_abs": float(np.max(np.abs(jac))) if finite_jac else None,
    }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Polish the identical objective with ROOT Minuit2 so the signal fit has a
  # real EDM and covariance diagnostic. This is opt-in and leaves legacy fits
  # unchanged unless --runMinuitPolish is requested.
  def _runMinuitPolish(self,objective,xbounds):
    names = [self.FitParameters[i].GetName() for i in range(len(self.FitParameters))]
    start = self.extractX0()

    def minuitObjective(values):
      point = np.asarray([values[i] for i in range(len(start))],dtype=float)
      return float(objective(point,self))

    minuit = ROOT.Math.Factory.CreateMinimizer("Minuit2","Migrad")
    if not minuit:
      raise RuntimeError("Cannot create ROOT Minuit2 minimizer")
    error_def = 0.5 if self.fitObjective == "unbinnedNLL" else 1.0
    minuit.SetErrorDef(error_def)
    minuit.SetStrategy(1)
    minuit.SetMaxFunctionCalls(200000)
    minuit.SetMaxIterations(200000)
    minuit.SetPrintLevel(0)
    # Minuit2 uses tolerance in its EDM stopping rule. The explicit EDM gate
    # below remains authoritative.
    minuit.SetTolerance(self.minuitPolishTolerance)
    functor = ROOT.Math.Functor(ROOT.std.function["double(const double*)"](minuitObjective),len(start))
    minuit.SetFunction(functor)
    for i,name in enumerate(names):
      low,high = float(xbounds[i][0]),float(xbounds[i][1])
      minuit.SetLimitedVariable(
        i,name,float(start[i]),max((high-low)*0.01,1.e-6),low,high
      )
    minimum_valid = bool(minuit.Minimize())
    hesse_valid = bool(minuit.Hesse())
    values = minuit.X()
    errors = minuit.Errors()

    for i,name in enumerate(names):
      self.FitParameters[i].setVal(float(values[i]))

    covariance_status = int(minuit.CovMatrixStatus())
    parameters = []
    parameters_at_limit = False
    for i,name in enumerate(names):
      value = float(values[i])
      error = float(errors[i])
      low,high = float(xbounds[i][0]),float(xbounds[i][1])
      span = high-low
      absolute_margin = min(value-low,high-value)/span if span > 0. else None
      if absolute_margin is not None and absolute_margin <= 1.e-5:
        parameters_at_limit = True
      parameters.append({
        "name": name,
        "value": value,
        "error": error,
        "min": low,
        "max": high,
        "absolute_fractional_margin": absolute_margin,
        "uncertainty_margin": (
          min(value-low,high-value)/error if error > 0. else None
        ),
      })
    edm = float(minuit.Edm())
    self.OptimizerDiagnostics["minuit"] = {
      "backend": "ROOT::Math::Minuit2",
      "valid": bool(minimum_valid and int(minuit.Status()) == 0),
      "status": int(minuit.Status()),
      "edm": edm,
      "edm_goal": float(self.strictEdmThreshold),
      "minimizer_tolerance": self.minuitPolishTolerance,
      "has_covariance": covariance_status >= 1,
      "has_accurate_covariance": covariance_status == 3,
      "has_positive_definite_covariance": covariance_status >= 2,
      "has_made_positive_definite_covariance": covariance_status == 2,
      "covariance_status": covariance_status,
      "hesse_failed": not hesse_valid,
      "has_parameters_at_limit": parameters_at_limit,
      "function_evaluations": int(minuit.NCalls()),
      "objective": float(minuit.MinValue()),
      "parameters": parameters,
    }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getOptimizerDiagnostics(self):
    diagnostics = od(self.OptimizerDiagnostics)
    diagnostics["fixed_parameters"] = [
      {
        "name": name,
        "value": float(var.getVal()),
        "min": float(var.getMin()),
        "max": float(var.getMax()),
      }
      for name,var in self.Vars.items()
      if hasattr(var,"isConstant") and var.isConstant()
    ]
    diagnostics["fit"] = {
      "name": self.name,
      "process": self.proc,
      "category": self.cat,
      "flavour": self.flav,
      "objective": self.fitObjective,
      "unbinned_weight_mode": self.unbinnedWeightMode,
      "bins": int(self.nBins),
      "fit_range": [float(self.fitRange[0]),float(self.fitRange[1])],
      "chi2": float(self.Chi2) if self.Chi2 is not None else None,
      "ndof": int(self.Ndof) if self.Ndof is not None else None,
      "reduced_chi2": (
        float(self.Chi2)/int(self.Ndof)
        if self.Chi2 is not None and self.Ndof else None
      ),
      "nll": float(self.NLL) if self.NLL is not None else None,
    }
    return diagnostics

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def passesStrictConvergence(self):
    minuit = self.OptimizerDiagnostics.get("minuit")
    if not minuit:
      return False
    return (
      minuit["valid"]
      and minuit["has_accurate_covariance"]
      and minuit["has_positive_definite_covariance"]
      and not minuit["hesse_failed"]
      and not minuit["has_parameters_at_limit"]
      and math.isfinite(minuit["edm"])
      and minuit["edm"] < self.strictEdmThreshold
    )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function build RooSpline1D from Polynomials to model MH dependence on fit params
  def buildSplines(self):
    # Loop over polynomials
    for k, poly in self.Polynomials.items():
      _x, _y = [], []
      _mh = 100.
      while(_mh<180.1):
        self.MH.setVal(_mh)
        _x.append(_mh)
        _y.append(poly.getVal())
        _mh += 0.1
      # Convert to arrays
      arr_x, arr_y = array('f',_x), array('f',_y)
      # Create spline and save to dict
      self.Splines[k] = ROOT.RooSpline1D(poly.GetName(),poly.GetName(),self.MH,len(_x),arr_x,arr_y)
           
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to print fit values
  def printFitParameters(self,title="Fit"):
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(" --> (%s) %s parameter values:"%(self.name,title))
    for i in range(len(self.FitParameters)):
      if self.FitParameters[i].GetName() == self.MH.GetName(): continue
      print("    * %-20s = %.6f"%(self.FitParameters[i].GetName(),self.FitParameters[i].getVal()))
    print("    ~~~~~~~~~~~~~~~~")
    print("    * chi2 = %.6f, n(dof) = %g --> chi2/n(dof) = %.3f"%(self.getChi2(),int(self.Ndof),self.getChi2()/int(self.Ndof)))
    print("    ~~~~~~~~~~~~~~~~")
    print("    * [VERBOSE] chi2 = %.6f"%(self.getChi2(verbose=True)))
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to set vars from json file
  def setVars(self,_json):
    with open(_json,"r") as jf: _vals = json.load(jf)
    for k,v in _vals.items(): self.Vars[k].setVal(v)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set fit-parameter starting values and optionally fix explicitly documented
  # weak parameters.  Parameters remain floating unless "constant": true is
  # present in the seed JSON.
  def setFitParameterValues(self,_json):
    with open(_json,"r") as jf: values = json.load(jf)
    for name,spec in values.items():
      if name not in self.Vars:
        raise KeyError("Unknown signal fit seed parameter: %s"%name)
      var = self.Vars[name]
      if isinstance(spec,dict):
        value = float(spec.get('value',var.getVal()))
        lower = float(spec.get('min',var.getMin()))
        upper = float(spec.get('max',var.getMax()))
        if lower >= upper:
          raise ValueError("Signal fit range %s has min >= max"%name)
        var.setMin(lower)
        var.setMax(upper)
      else:
        value = float(spec)
      if value < var.getMin() or value > var.getMax():
        raise ValueError("Signal fit seed %s=%g is outside [%g,%g]"%(
          name,value,var.getMin(),var.getMax()
        ))
      var.setVal(value)
      var.setConstant(bool(spec.get("constant",False)) if isinstance(spec,dict) else False)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to re-calculate chi2 after setting vars
  def getChi2(self,verbose=False):
    x = self.extractX0()
    self.Chi2 = nChi2Addition(x,self,verbose=verbose)
    return self.Chi2

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to re-calculate the summed unbinned NLL after setting vars
  def getNLL(self):
    x = self.extractX0()
    self.NLL = nNLLAddition(x,self)
    return self.NLL

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function to re-calculate chi2/ndof after setting vars
  def getReducedChi2(self):
    x = self.extractX0()
    self.Chi2 = nChi2Addition(x,self)
    return self.Chi2/int(self.Ndof)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

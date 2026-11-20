import ROOT
import os
import sys
import json
import re
import numpy as np
import pandas
import pickle
from collections import OrderedDict as od
from commonObjects import *
from commonTools import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
class FinalModel:
  # Constructor
  def __init__(self,_ssfMap,_proc,_cat,_flav,_ext,_year,_sqrts,_datasets,_xvar,_MH,_MHLow,_MHHigh,_massPoints,_xsbrMap,_procSyst,_scales,_scalesCorr,_scalesGlobal,_smears,_doVoigtian,_useDCB,_skipSystematics,_doEffAccFromJson,_systematicsFile=None):
    self.ssfMap = _ssfMap
    self.proc = _proc
    self.procSyst = _procSyst # Signal process used for systematics (useful for low stat cases)
    self.cat = _cat
    self.flav = _flav
    self.ext = _ext
    self.year = _year
    self.sqrts = _sqrts
    self.name = "%s_%s_%s_%s_%s"%(self.proc,self.year,self.cat,self.flav,self.sqrts)
    self.datasets = _datasets
    self.xvar = _xvar
    self.aset = ROOT.RooArgSet(self.xvar)
    self.MH = _MH
    self.MHLow = _MHLow
    self.MHHigh = _MHHigh
    self.massPoints = _massPoints
    self.intLumi = ROOT.RooRealVar("IntLumi","IntLumi",1.,0.,999999999.) # in pb^-1
    self.xsbrMap = _xsbrMap
    # Systematics
    self.skipSystematics = _skipSystematics
    self.scales = _scales
    self.scalesCorr = _scalesCorr
    self.scalesGlobal = _scalesGlobal
    self.smears = _smears
    # Options:
    self.useDCB = _useDCB
    self.doVoigtian = _doVoigtian
    if self.doVoigtian: self.GammaH = ROOT.RooRealVar("GammaH","GammaH",0.004,0.,5.)
    self.doEffAccFromJson = _doEffAccFromJson
    self.verbose = True
    # Dict to store objects
    self.Splines = od()
    self.Functions = od()
    self.NuisanceMap = od()
    self.Pdfs = od()
    self.Datasets = od()
    # Build XS/BR/EA splines
    self.buildXSBRSplines()
    self.buildEffAccSpline()
    # If not skip systematics: add nuisance params to dict
    if not self.skipSystematics: self.buildNuisanceMap(_systematicsFile)
    # Build final pdfs
    self.buildPdf(self.ssfMap['Total'],ext='total',useDCB=self.useDCB)
    self.Pdfs['final'] = self.Pdfs['total']
    # Build final normalisation, datasets and extended Pdfs
    self.buildNorm()
    self.buildDatasets()
    self.buildExtended() 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Functions to get XS, BR and EA splines for given proc/decay from map
  def buildXSBRSplines(self):
    mh = np.linspace(120.,130.,101)
    # These inputs carry their physical weights. Keep the spline objects and
    # the factor of 1000 used by the existing HZG workspace/datacard contract.
    proc = self.xsbrMap[self.proc]
    if proc["mode"] != "constant":
      raise ValueError("HZG preweighted inputs require a constant process factor")
    xs = np.full_like(mh, proc.get("factor", 1.0))
    self.Splines['xs'] = ROOT.RooSpline1D("fxs_%s_%s"%(self.proc,self.sqrts),"fxs_%s_%s"%(self.proc,self.sqrts),self.MH,len(mh),mh,xs)
    br = np.full_like(mh, self.xsbrMap["decay"].get("factor", 1.0) / 1000.)
    self.Splines['br'] = ROOT.RooSpline1D("fbr_%s"%self.sqrts,"fbr_%s"%self.sqrts,self.MH,len(mh),mh,br)

  def buildEffAccSpline(self):
    # Two treatments: load from json created with getEffAcc.py script or calc from sum of weights
    # Loop over mass points  
    ea, mh = [], []
    for mp in self.massPoints.split(","):
      mh.append(float(mp))
      if self.doEffAccFromJson:
        jfname = "%s/outdir_%s/getEffAcc/json/effAcc_M%s_%s.json"%(swd__,self.ext,mp,self.ext)
        if not os.path.exists(jfname):
          print(" --> [ERROR] effAcc json file (%s) does not exist for mass point = %s. Run getEffAcc first."%(jfname,mp))
          sys.exit(1)
        with open(jfname,'r') as jf: ea_data = json.load(jf)
        ea.append(float(ea_data['%s_%s_%s'%(self.proc,self.cat,self.flav)]))
      else:
        sumw = self.datasets[mp].sumEntries()
        ea.append(sumw) 
    # If single mass point then add MHLow and MHHigh dummy points for constant ea
    if len(ea) == 1: ea, mh = [ea[0],ea[0],ea[0]], [float(self.MHLow),mh[0],float(self.MHHigh)]
    # Convert to numpy arrays and make spline
    ea, mh = np.asarray(ea), np.asarray(mh)
    self.Splines['ea'] = ROOT.RooSpline1D("fea_%s"%(self.name),"fea_%s"%(self.name),self.MH,len(mh),mh,ea)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to build final normalisation: XS * BR * eff * acc * rate
  def buildNorm(self):
    # Build rate function: encode affect of nuisances on signal rate
    self.buildRate("rate_%s"%self.name,skipSystematics=self.skipSystematics)
    finalPdfName = self.Pdfs['final'].GetName()
    self.Functions['final_norm'] = ROOT.RooFormulaVar("%s_norm"%finalPdfName,"%s_norm"%finalPdfName,"@0*@1*@2*@3",ROOT.RooArgList(self.Splines['xs'],self.Splines['br'],self.Splines['ea'],self.Functions['rate_%s'%self.name]))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function for making nuisance param w/ info to add to Nuisance dict
  def makeNuisance(self,nuisanceName,nuisanceMeanConst,nuisanceSigmaConst,nuisanceRateConst,nuisanceType,nuisanceOpts=[]):
    self.NuisanceMap[nuisanceType][nuisanceName] = {
      'name':nuisanceName,
      # 'param':ROOT.RooRealVar("%s_%s"%(outputWSNuisanceTitle__,nuisanceName),"%s_%s"%(outputWSNuisanceTitle__,nuisanceName),0.,-5.,5.),
      'param':ROOT.RooRealVar(nuisanceName,nuisanceName,0.,-5.,5.),
      'meanConst':ROOT.RooConstVar("const_%s_mean_%s"%(self.name,nuisanceName),"const_%s_mean_%s"%(self.name,nuisanceName),nuisanceMeanConst),
      'sigmaConst':ROOT.RooConstVar("const_%s_sigma_%s"%(self.name,nuisanceName),"const_%s_sigma_%s"%(self.name,nuisanceName),nuisanceSigmaConst),
      'rateConst':ROOT.RooConstVar("const_%s_rate_%s"%(self.name,nuisanceName),"const_%s_rate_%s"%(self.name,nuisanceName),nuisanceRateConst),
      'opts':nuisanceOpts}
    self.NuisanceMap[nuisanceType][nuisanceName]['param'].setConstant(True)
  
  # Function for building Nuisance param map:
  def buildNuisanceMap(self,psname):
    # Dict to store nuisances of different type in map
    for sType in ['scales','scalesCorr','scalesGlobal','smears']:
      if getattr(self,sType) != '': self.NuisanceMap[sType] = od()

    # Extract calcPhotonSyst output
    if psname is None:
      psname = os.path.join(swd__, "output_" + self.ext, "calcPhotonSyst", "pkl", self.cat + ".pkl")
    with open(psname, "rb") as stream:
      table = pickle.load(stream)
    rows = table.loc[table["proc"] == self.procSyst]
    for column, value in (("cat", self.cat), ("flav", self.flav)):
      if column in table:
        rows = rows.loc[rows[column] == value]
    if len(rows) != 1:
      raise ValueError(f"{psname}: expected one systematics row for process {self.procSyst}, "
                       f"category {self.cat}, flavor {self.flav}; found {len(rows)}")
    row = rows.iloc[0]

    for systType in self.NuisanceMap:
      for spec in getattr(self, systType).split(","):
        if not spec: continue
        name, *opts = spec.split(":")
        if systType == "scalesGlobal":
          mean, sigma, rate = 0., 0., 0.
        else:
          # Older calcPhotonSyst tables qualify columns with energy and type.
          suffix = {"scales": "scale", "scalesCorr": "scaleCorr", "smears": "smear"}[systType]
          prefix = name if name + "_mean" in row else f"{name}_{self.sqrts}{suffix}"
          mean, sigma, rate = (float(row[f"{prefix}_{field}"]) for field in ("mean", "sigma", "rate"))
          if not np.isfinite([mean, sigma, rate]).all():
            raise ValueError(f"{psname}: non-finite {name} coefficient for {self.procSyst}/{self.cat}/{self.flav}")
        self.makeNuisance(name, mean, sigma, rate, systType, opts)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to get RV fraction func
  def rebindSplineToMH(self,ssf,key,name):
    oldMH = ssf.MH.getVal()
    xs = np.linspace(100.,180.,801)
    ys = []
    for mh in xs:
      ssf.MH.setVal(float(mh))
      ys.append(ssf.Splines[key].getVal())
    ssf.MH.setVal(oldMH)
    return ROOT.RooSpline1D(name,name,self.MH,len(xs),xs,np.asarray(ys))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to build final PDFs from input SimultaneousFit object splines
  def buildPdf(self,ssf,ext='',useDCB=False,_recursive=True):
    extStr = "%s_%s"%(self.name,ext) if ext!='total' else '%s'%self.name
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # For double CB + Gaussian
    if useDCB:
      hasGaussian = 'sigma_gaus' in ssf.Splines and 'frac_constrained' in ssf.Splines
      # Extract splines
      for f in ['dm','sigma','n1','n2','a1','a2']:
        k = "%s_dcb"%f
        outName = "%s_%s"%(re.sub("sigma","sigma_fit",k),extStr)
        self.Splines["%s_%s"%(k,extStr)] = self.rebindSplineToMH(ssf,k,outName)
      # Build mean and sigma functions: including systematics
      self.buildMean('dm_dcb_%s'%extStr,skipSystematics=self.skipSystematics)
      self.buildSigma('sigma_dcb_%s'%extStr,skipSystematics=self.skipSystematics)
      # Build DCB pdf
      # In the mixture, keep the DSCB component name distinct from the final
      # RooAddPdf name.  Pure DSCB-SC retains the standard dcb_<extStr> name.
      dcbName = ("dcbcore_%s" if hasGaussian else "dcb_%s")%extStr
      self.Pdfs['dcb_%s'%extStr] = ROOT.RooDoubleCBFast(dcbName,dcbName,self.xvar,self.Functions["mean_dcb_%s"%extStr],self.Functions["sigma_dcb_%s"%extStr],self.Splines['a1_dcb_%s'%extStr],self.Splines['n1_dcb_%s'%extStr],self.Splines['a2_dcb_%s'%extStr],self.Splines['n2_dcb_%s'%extStr])
      
      # + Gaussian: shares mean with DCB
      if hasGaussian:
        # Original DCB+Gaussian model
        self.Splines['sigma_gaus_%s'%extStr] = self.rebindSplineToMH(
          ssf,'sigma_gaus',"sigma_fit_gaus_%s"%extStr
        )
        self.buildSigma('sigma_gaus_%s'%extStr,skipSystematics=self.skipSystematics)
        if self.doVoigtian:
          self.Pdfs['gaus_%s'%extStr] = ROOT.RooVoigtian("gaus_%s"%extStr,"gaus_%s"%extStr,self.xvar,self.Functions["mean_dcb_%s"%extStr],self.GammaH,self.Functions["sigma_gaus_%s"%extStr])
        else:
          self.Pdfs['gaus_%s'%extStr] = ROOT.RooGaussian("gaus_%s"%extStr,"gaus_%s"%extStr,self.xvar,self.Functions["mean_dcb_%s"%extStr],self.Functions["sigma_gaus_%s"%extStr])

        # Fraction
        self.Splines['frac_%s'%extStr] = self.rebindSplineToMH(
          ssf,'frac_constrained',"frac_%s"%extStr
        )

        # Define total pdf - DCB + Gaussian
        _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
        for pdf in ['dcb','gaus']: _pdfs.add(self.Pdfs['%s_%s'%(pdf,extStr)])
        _coeffs.add(self.Splines['frac_%s'%extStr])
        self.Pdfs[ext] = ROOT.RooAddPdf("%s_%s"%(outputWSObjectTitle__,extStr),"%s_%s"%(outputWSObjectTitle__,extStr),_pdfs,_coeffs,_recursive)
      else:
        # Use only DCB model, similar to the simplified version in simultaneousFit.py
        self.Pdfs[ext] = self.Pdfs['dcb_%s'%extStr]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # For nGaussians:
    else:
      # For total pdf
      _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()

      for g in range(0,ssf.nGaussians):
        # Extract splines
        for f in ['dm','sigma']:
          k = "%s_g%g"%(f,g)
          self.Splines["%s_%s"%(k,extStr)] = ssf.Splines[k].Clone()
          self.Splines["%s_%s"%(k,extStr)].SetName("%s_%s"%(re.sub("sigma","sigma_fit",k),extStr))
        # Build mean and sigma functions for gaussian, g: including systematics
        self.buildMean('dm_g%g_%s'%(g,extStr),skipSystematics=self.skipSystematics)
        self.buildSigma('sigma_g%g_%s'%(g,extStr),skipSystematics=self.skipSystematics)
        # Build Gaussian
        if self.doVoigtian: 
          self.Pdfs['gaus_g%g_%s'%(g,extStr)] = ROOT.RooVoigtian("gaus_g%g_%s"%(g,extStr),"gaus_g%g_%s"%(g,extStr),self.xvar,self.Functions["mean_g%g_%s"%(g,extStr)],self.GammaH,self.Functions["sigma_g%g_%s"%(g,extStr)])
        else: 
          self.Pdfs['gaus_g%g_%s'%(g,extStr)] = ROOT.RooGaussian("gaus_g%g_%s"%(g,extStr),"gaus_g%g_%s"%(g,extStr),self.xvar,self.Functions["mean_g%g_%s"%(g,extStr)],self.Functions["sigma_g%g_%s"%(g,extStr)])
        _pdfs.add(self.Pdfs['gaus_g%g_%s'%(g,extStr)])
        
        # Fractions
        if g < ssf.nGaussians-1:
          self.Splines['frac_g%g_%s'%(g,extStr)] = ssf.Splines['frac_g%g_constrained'%g]
          self.Splines['frac_g%g_%s'%(g,extStr)].SetName("frac_g%g_%s"%(g,extStr))
          _coeffs.add(self.Splines['frac_g%g_%s'%(g,extStr)])

      # Define total pdf. Keep the nGaussian branch on the standard
      # non-recursive coefficient convention used by the fixed-Gaussian fits.
      self.Pdfs[ext] = ROOT.RooAddPdf("%s_%s"%(outputWSObjectTitle__,extStr),"%s_%s"%(outputWSObjectTitle__,extStr),_pdfs,_coeffs)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Functions to build mean, sigma and rate functions with systematics
  def buildMean(self,dmSplineName="",skipSystematics=False):
    meanName = re.sub("dm","mean",dmSplineName) 
    # Build formula string and dependents list
    dependents = ROOT.RooArgList()
    formula = "(@0+@1)"
    dependents.add(self.MH)
    dependents.add(self.Splines[dmSplineName])
    if not skipSystematics:
      # Add systematics
      formula += "*(1."
      # Global
      if 'scalesGlobal' in self.NuisanceMap:
        for sName, sInfo in self.NuisanceMap['scalesGlobal'].items():
          formula += "+@%g"%dependents.getSize()
          # For adding additional factor
          for so in sInfo['opts']: 
            if "factor_%s"%self.cat in so:
              additionalFactor = float(so.split("=")[-1])
              formula += "*%3.1f"%additionalFactor
          dependents.add(sInfo['param'])
      # Other systs: scales, scalesCorr, smears
      for sType in ['scales','scalesCorr','smears']:
        if sType in self.NuisanceMap:
          for sName, sInfo in self.NuisanceMap[sType].items():
            c = sInfo['meanConst'].getVal()
            if abs(c)>=5.e-5:
              formula += "+@%g*@%g"%(dependents.getSize(),dependents.getSize()+1)
              dependents.add(sInfo['meanConst'])
              dependents.add(sInfo['param'])
      formula += ")"
    self.Functions[meanName] = ROOT.RooFormulaVar(meanName,meanName,formula,dependents)

  def buildSigma(self,sigmaSplineName="",skipSystematics=False):
    sigmaName = sigmaSplineName
    # Build formula string and dependents list
    dependents = ROOT.RooArgList()
    formula = "@0"
    dependents.add(self.Splines[sigmaSplineName])
    if not skipSystematics:
      # Add systematics
      formula += "*TMath::Max(1.e-2,(1."
      for sType in ['scales','scalesCorr','smears']:
        if sType in self.NuisanceMap:
          for sName, sInfo in self.NuisanceMap[sType].items():
            c = sInfo['sigmaConst'].getVal()
            if c>=1e-4:
              formula += "+@%g*@%g"%(dependents.getSize(),dependents.getSize()+1)
              dependents.add(sInfo['sigmaConst'])
              dependents.add(sInfo['param'])
      formula += "))"
    self.Functions[sigmaName] = ROOT.RooFormulaVar(sigmaName,sigmaName,formula,dependents)

  def buildRate(self,rateName,skipSystematics=False):
    # Build formula string and dependents list
    dependents = ROOT.RooArgList()
    formula = "(1."
    if not skipSystematics:
      # Add systematics
      for sType in ['scales','scalesCorr','smears']:
        if sType in self.NuisanceMap:
          for sName, sInfo in self.NuisanceMap[sType].items():
            c = sInfo['rateConst'].getVal()
            if c>=5.e-4:
              formula += "+@%g*@%g"%(dependents.getSize(),dependents.getSize()+1)
              dependents.add(sInfo['rateConst'])
              dependents.add(sInfo['param'])
    formula += ")"
    self.Functions[rateName] = ROOT.RooFormulaVar(rateName,rateName,formula,dependents)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to build datasets to add to workspace
  def buildDatasets(self):
    for mp, d in self.datasets.items(): 
      self.Datasets[mp] = d.Clone("sig_mass_m%s_%s"%(mp,self.name))
      self.Datasets['%s_copy'%mp] = d.Clone("sig_mass_m%s_%s"%(mp,self.cat))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to build extended Pdfs and normalisation with luminosity
  def buildExtended(self):
    finalPdfName = self.Pdfs['final'].GetName()
    self.Functions['final_normThisLumi'] = ROOT.RooFormulaVar("%s_normThisLumi"%finalPdfName,"%s_normThisLumi"%finalPdfName,"@0*@1*@2*@3*@4",ROOT.RooArgList(self.Splines['xs'],self.Splines['br'],self.Splines['ea'],self.Functions['rate_%s'%self.name],self.intLumi))
    self.Pdfs['final_extend'] = ROOT.RooExtendPdf("extend%s"%finalPdfName,"extend%s"%finalPdfName,self.Pdfs['final'],self.Functions['final_norm'])
    self.Pdfs['final_extendThisLumi'] = ROOT.RooExtendPdf("extend%sThisLumi"%finalPdfName,"extend%sThisLumi"%finalPdfName,self.Pdfs['final'],self.Functions['final_normThisLumi'])
 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function for saving to output workspace
  def save(self,wsout):
    wsout.imp = getattr(wsout,"import")
    wsout.imp(self.Pdfs['final'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['final_norm'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['final_normThisLumi'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['final_extend'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['final_extendThisLumi'],ROOT.RooFit.RecycleConflictNodes())
    for d in self.Datasets.values(): wsout.imp(d)

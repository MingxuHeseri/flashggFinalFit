from collections import OrderedDict as od

globalXSBRMap = od()
globalXSBRMap['fiducialAnalysis'] = od()
globalXSBRMap['fiducialAnalysis']['decay'] = {'mode': "hzg"}
globalXSBRMap['fiducialAnalysis']['GG2H_in'] = {'mode':'constant','factor':51.96}
globalXSBRMap['fiducialAnalysis']['VBF_in'] = {'mode':'constant','factor':4.067}
globalXSBRMap['fiducialAnalysis']['VH_in'] = {'mode':'constant','factor':2.3781}
globalXSBRMap['fiducialAnalysis']['TTH_in'] = {'mode':'constant','factor':0.5638}
globalXSBRMap['fiducialAnalysis']['GG2H_out'] = {'mode':'constant','factor':51.96}
globalXSBRMap['fiducialAnalysis']['VBF_out'] = {'mode':'constant','factor':4.067}
globalXSBRMap['fiducialAnalysis']['VH_out'] = {'mode':'constant','factor':2.3781}
globalXSBRMap['fiducialAnalysis']['TTH_out'] = {'mode':'constant','factor':0.5638}
globalXSBRMap['fiducialAnalysis']['ggH'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['VBF'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['WH'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['WplusH'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['WminusH'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['ZH'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['ttH'] = {'mode':'constant','factor': 1}
globalXSBRMap['fiducialAnalysis']['all'] = {'mode':'constant','factor': 1}

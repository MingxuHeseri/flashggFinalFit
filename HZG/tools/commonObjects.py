"""Workspace names and luminosity convention for preweighted HZG inputs."""
import os

cwd__ = os.environ.get('FLASHGGFINALFIT_HZG_BASE', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
swd__ = "%s/Signal"%cwd__
sqrts__ = "13TeV"
lumiMap = {'2016preVFP':1 , '2016postVFP': 1, '2016': 1, '2017':1, '2018':1, 'combined':1, 'merged':1, '2022':1, '2022preEE':1, '2022postEE':1,'2023preBPix':1, '2023postBPix':1, 'all':1}
lumiScaleFactor = 1000.
inputWSName__ = "CMS_hzg_workspace"
outputWSName__ = "wsig"
outputWSObjectTitle__ = "dcb"

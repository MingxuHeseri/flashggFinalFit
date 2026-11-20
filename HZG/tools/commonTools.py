"""Process names used in HZG input datasets."""

procToDataMap = {
    "GG2H": "ggh",
    "VH": "vh",
    "GG2H_in": "ggh_in",
    "VBF_in": "vbf_in",
    "VH_in": "vh_in",
    "GG2H_out": "ggh_out",
    "VBF_out": "vbf_out",
    "VH_out": "vh_out",
    "WH2HQQ": "wh",
    "ZH2HQQ": "zh",
    "QQ2HLNU": "wh",
    "QQ2HLL": "zh",
    "TTH": "tth",
    "TTH_in": "tth_in",
    "TTH_out": "tth_out",
    "BBH": "bbh",
    "THQ": "thq",
    "THW": "thw",
    "GG2HQQ": "ggzh",
    "GG2HLL": "ggzh",
    "GG2HNUNU": "ggzh",
}


def procToData(proc):
    return procToDataMap.get(proc, proc)

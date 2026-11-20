"""Rebuild a frozen FinalFit envelope with the CombinedLimit native erf factory."""
import re

# Displayed orders for the two mixtures count components; the registry counts
# independent shape coefficients. Bernstein includes the fixed leading term.
FAMILIES = {
    "bernerf": ("BernsteinErfTurnOn", "RooErfProductPdf"),
    "experf": ("ExponentialErfTurnOn", "RooErfProductPdf"),
    "powerf": ("PowerLawErfTurnOn", "RooPowerErfPdf"),
    "lauerf": ("LaurentErfTurnOn", "RooPowerErfPdf"),
    "expbernerf": ("ExpBernsteinErfTurnOn", "RooExpPolyErfPdf"),
    "expcheberf": ("ExpChebychevErfTurnOn", "RooExpPolyErfPdf"),
}


def copy_parameter(src, dst):
    if src.hasMin():
        dst.setMin(src.getMin())
    else:
        dst.removeMin()
    if src.hasMax():
        dst.setMax(src.getMax())
    else:
        dst.removeMax()
    dst.setVal(src.getVal())
    dst.setConstant(src.isConstant())
    if src.hasError():
        dst.setError(src.getError())
    else:
        dst.removeError()
    if src.hasAsymError():
        dst.setAsymError(src.getAsymErrorLo(), src.getAsymErrorHi())
    else:
        dst.removeAsymError()


def native_envelope(root, source, category, index_name, frozen=None):
    """Return an upstream multipdf workspace ready for a standard shapes card.

    The frozen bias workspace, when supplied, defines fit starts, member order
    and corrections. The source handoff supplies the observable, blinded data
    and rate object with their original pre-card names.
    """
    name = f"CMS_hzg_{category}_13_13p6TeV_bkgshape"
    data = source.data("roohist_data_mass_" + category)
    mass = source.var(data.get().first().GetName())
    source_pdf = source.pdf(name)
    reference = frozen.pdf("shapeBkg_bkg_mass_" + category) if frozen else source_pdf
    state = frozen if frozen else source
    index = state.cat(index_name)
    nominal = index.getIndex()
    ref_mass = state.var(mass.GetName())
    obs = root.RooArgSet(mass)
    ref_obs = root.RooArgSet(ref_mass)
    builder = root.PdfModelBuilder()
    builder.setObsVar(mass)
    by_name = {source_pdf.getPdf(i).GetName() for i in range(source_pdf.getNumPdfs())}
    members = root.RooArgList()
    corrections = []
    try:
        for i in range(reference.getNumPdfs()):
            old = reference.getPdf(i)
            pdf_name = old.GetName()
            match = re.fullmatch(r"(.+)_(bernerf|experf|powerf|lauerf|expbernerf|expcheberf)(\d+)", pdf_name)
            if not match or pdf_name not in by_name:
                raise ValueError(f"Unsupported or missing source member {pdf_name}")
            prefix, family, display = match.groups()
            order = int(display)
            if family == "bernerf":
                order += 1
            elif family in ("experf", "powerf"):
                order = 2 * order - 1
            method, cls = FAMILIES[family]
            pdf = root.getPdf(builder, method, order, root.std.string(), prefix, category, False, False)
            if not pdf or pdf.GetName() != pdf_name or pdf.ClassName() != cls:
                raise ValueError(f"Factory mismatch for {pdf_name}")
            params = {v.GetName(): v for v in pdf.getParameters(obs)}
            ref = {v.GetName(): v for v in old.getParameters(ref_obs)}
            if set(params) != set(ref):
                raise ValueError(f"{pdf_name}: parameter mismatch {set(params) ^ set(ref)}")
            for key, var in params.items():
                copy_parameter(ref[key], var)
            # Persist the precision used to qualify these native candidates.
            # Reuse the factory-selected method, including adaptive integration.
            config = root.RooNumIntConfig(pdf.getIntegratorConfig())
            config.setEpsAbs(1e-12)
            config.setEpsRel(1e-12)
            pdf.setIntegratorConfig(config)
            members.add(pdf)
            index.setIndex(i)
            corrections.append(reference.getCorrection())
    finally:
        index.setIndex(nominal)
    idx = root.RooCategory(index_name, index_name)
    multi = root.RooMultiPdf(name, name, idx, members)
    for i, expected in enumerate(corrections):
        idx.setIndex(i)
        if multi.getCorrection() != expected or idx.getLabel() != index.lookupName(i):
            raise ValueError(f"{category}: discrete correction or label mismatch at {i}")
    idx.setIndex(nominal)
    workspace = root.RooWorkspace("multipdf")
    for obj in (multi, source.var(name + "_norm"), data):
        if getattr(workspace, "import")(obj):
            raise ValueError(f"Failed to import {obj.GetName()}")
    return workspace


def scoped_envelope(root, source, category, suffix):
    """Give independent category shape graphs unique names before combination.

    Historical campaigns reuse ``env_pdf_catN`` across different category
    schemes. Combining those files without scoping can identify independent
    parameters by name. Keep the public mass, index, MultiPdf, data and rate
    names; scope every member graph and its shape variables. This changes no
    values, bounds, flags, member ordering or discrete correction.
    """
    name = f"CMS_hzg_{category}_13_13p6TeV_bkgshape"
    index_name = f"pdfindex_{category}_13_13p6TeV"
    data = source.data("roohist_data_mass_" + category)
    mass_name = data.get().first().GetName()
    original = source.pdf(name)
    original_index = source.cat(index_name)
    nominal = original_index.getIndex()
    workspace = root.RooWorkspace("multipdf")
    members = root.RooArgList()
    corrections = []
    try:
        for i in range(original.getNumPdfs()):
            pdf = original.getPdf(i)
            if getattr(workspace, "import")(
                pdf,
                root.RooFit.RenameAllNodes(suffix),
                root.RooFit.RenameAllVariablesExcept(suffix, mass_name),
                root.RooFit.RecycleConflictNodes(),
                root.RooFit.Silence(),
            ):
                raise ValueError(f"Cannot scope {pdf.GetName()}")
            members.add(workspace.pdf(pdf.GetName() + "_" + suffix))
            original_index.setIndex(i)
            corrections.append(float(original.getCorrection()))
    finally:
        original_index.setIndex(nominal)
    index = root.RooCategory(index_name, index_name)
    multi = root.RooMultiPdf(name, name, index, members)
    for i, correction in enumerate(corrections):
        index.setIndex(i)
        if float(multi.getCorrection()) != correction:
            raise ValueError(f"Scoping changed correction for {category}/{i}")
    index.setIndex(nominal)
    for obj in (multi, source.var(name + "_norm")):
        if getattr(workspace, "import")(
            obj, root.RooFit.RecycleConflictNodes(), root.RooFit.Silence()
        ):
            raise ValueError(f"Cannot import {obj.GetName()}")
    if getattr(workspace, "import")(data, root.RooFit.Silence()):
        raise ValueError(f"Cannot import {data.GetName()}")
    return workspace

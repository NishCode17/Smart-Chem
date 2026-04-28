"""
backend/admet.py
================
Comprehensive ADMET prediction using RDKit physicochemical descriptors.
No extra dependencies beyond rdkit.

Absorption  : HIA, Caco-2, oral bioavailability, water solubility
Distribution: BBB penetration, PPB, VDss
Metabolism  : CYP3A4 / CYP2D6 / CYP2C9 inhibition
Excretion   : Renal clearance, half-life
Toxicity    : hERG, AMES, DILI, skin sensitisation

Rules from: Lipinski 1997, Ertl 2000, Veber 2002, Gleeson 2008, Di 2013.
"""

import math
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _sig(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))

def _clamp(v: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, v))

def _gauss(x: float, mu: float, sigma: float) -> float:
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ---------------------------------------------------------------------------
# Descriptor extraction (called once, reused by all sub-predictors)
# ---------------------------------------------------------------------------
def _desc(mol) -> dict:
    return dict(
        mw    = Descriptors.ExactMolWt(mol),
        logp  = Descriptors.MolLogP(mol),
        tpsa  = Descriptors.TPSA(mol),
        hbd   = Descriptors.NumHDonors(mol),
        hba   = Descriptors.NumHAcceptors(mol),
        rb    = Descriptors.NumRotatableBonds(mol),
        rings = rdMolDescriptors.CalcNumRings(mol),
        ar    = rdMolDescriptors.CalcNumAromaticRings(mol),
        ha    = mol.GetNumHeavyAtoms(),
        fsp3  = rdMolDescriptors.CalcFractionCSP3(mol),
        charge= Chem.GetFormalCharge(mol),
    )


# ---------------------------------------------------------------------------
# A — Absorption
# ---------------------------------------------------------------------------
def _absorption(d: dict) -> dict:
    logp, tpsa, mw, hbd, hba, rb, ar = (
        d["logp"], d["tpsa"], d["mw"],
        d["hbd"],  d["hba"],  d["rb"], d["ar"]
    )

    # HIA — Veber rules + TPSA model (Ertl 2000)
    hia = _clamp(
        1.0
        - _sig((tpsa - 100) / 20) * 0.6
        - (0.1 if rb > 10 else 0)
        - (0.1 if mw > 500 else 0)
    )

    # Caco-2 permeability (Artursson 1996)
    caco2 = _clamp(_gauss(logp, 2.0, 1.5) * _clamp(1 - tpsa / 180))

    # Oral bioavailability — RO5 proxy
    ro5v = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    oral = _clamp(1.0 - ro5v * 0.2 - (0.15 if tpsa > 140 else 0))

    # Water solubility — Delaney ESOL (2004)
    log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ar

    return {
        "hia_score":        round(hia, 3),
        "hia_class":        "High" if hia >= 0.7 else ("Moderate" if hia >= 0.4 else "Low"),
        "caco2_score":      round(caco2, 3),
        "caco2_class":      "High" if caco2 >= 0.6 else ("Moderate" if caco2 >= 0.3 else "Low"),
        "oral_f_score":     round(oral, 3),
        "oral_f_class":     "High (>70%)" if oral >= 0.75 else ("Moderate" if oral >= 0.4 else "Low (<30%)"),
        "log_solubility":   round(log_s, 2),
        "solubility_class": ("Highly soluble" if log_s > -2 else
                             "Soluble"         if log_s > -4 else
                             "Moderately soluble" if log_s > -6 else
                             "Poorly soluble"),
    }


# ---------------------------------------------------------------------------
# D — Distribution
# ---------------------------------------------------------------------------
def _distribution(d: dict) -> dict:
    logp, tpsa, mw, charge = d["logp"], d["tpsa"], d["mw"], d["charge"]

    # BBB — CNS drugs: LogP 1-3, TPSA<90, MW<450 (Clark 1999)
    bbb = _clamp(
        _gauss(logp, 2.0, 1.5)
        * _clamp(1 - tpsa / 140)
        * (1.0 if mw <= 450 else _clamp(1 - (mw - 450) / 200))
        * (0.7 if abs(charge) > 0 else 1.0)
    )

    # Plasma protein binding — high logP -> high PPB (Kratochwil 2002)
    ppb = _clamp(_sig((logp - 1.5) / 1.2))

    # VDss rough estimate (L/kg)
    vd_score = _clamp(logp / 8.0 + (0.1 if charge > 0 else 0))
    vd_est   = round(0.5 + vd_score * 15, 1)

    return {
        "bbb_score":   round(bbb, 3),
        "bbb_class":   "CNS Penetrant" if bbb >= 0.6 else ("Low CNS" if bbb >= 0.3 else "CNS Impermeable"),
        "ppb_score":   round(ppb, 3),
        "ppb_class":   "High (>90%)" if ppb >= 0.7 else ("Moderate" if ppb >= 0.4 else "Low (<70%)"),
        "vd_estimate": vd_est,
    }


# ---------------------------------------------------------------------------
# M — Metabolism
# ---------------------------------------------------------------------------
def _metabolism(mol, d: dict) -> dict:
    logp, mw, ar, fsp3 = d["logp"], d["mw"], d["ar"], d["fsp3"]

    # CYP3A4 — lipophilic, aromatic, large (Gleeson 2008)
    cyp3a4 = _clamp(
        0.3 * _sig((mw - 350) / 80)
      + 0.3 * _sig((logp - 3.0) / 1.0)
      + 0.2 * _sig((ar - 2) / 1.0)
      + 0.2 * _sig((4 - fsp3 * 10) / 2)
    )

    # CYP2D6 — basic nitrogen compounds
    has_basic_n = any(
        a.GetAtomicNum() == 7 and a.GetTotalValence() <= 3 and a.GetFormalCharge() == 0
        for a in mol.GetAtoms()
    )
    cyp2d6 = _clamp(
        (0.4 if has_basic_n else 0.0)
      + 0.3 * _sig((logp - 2.0) / 1.5)
      + 0.3 * _sig((ar - 1) / 1.0)
    )

    # CYP2C9 — acidic/neutral, sulfas, NSAIDs
    cyp2c9 = _clamp(
        0.4 * _sig((logp - 2.5) / 1.2)
      + 0.3 * _sig((mw - 300) / 100)
      + 0.3 * _sig((ar - 1) / 1.5)
    )

    def _cls(v):
        return "Likely Inhibitor" if v >= 0.55 else ("Possible" if v >= 0.35 else "Unlikely")

    return {
        "cyp3a4_score": round(cyp3a4, 3), "cyp3a4_class": _cls(cyp3a4),
        "cyp2d6_score": round(cyp2d6, 3), "cyp2d6_class": _cls(cyp2d6),
        "cyp2c9_score": round(cyp2c9, 3), "cyp2c9_class": _cls(cyp2c9),
    }


# ---------------------------------------------------------------------------
# E — Excretion
# ---------------------------------------------------------------------------
def _excretion(d: dict) -> dict:
    logp, mw = d["logp"], d["mw"]

    renal = _clamp(_sig((3 - logp) / 1.5) * _sig((500 - mw) / 150))
    t_half = min(1.5 + logp * 2.5 + (mw / 500) * 5, 48.0)

    return {
        "renal_clearance_score": round(renal, 3),
        "renal_clearance_class": "High" if renal >= 0.6 else ("Moderate" if renal >= 0.35 else "Low"),
        "half_life_h":           round(t_half, 1),
        "half_life_class":       "Short (<3h)" if t_half < 3 else ("Medium (3-12h)" if t_half < 12 else "Long (>12h)"),
    }


# ---------------------------------------------------------------------------
# T — Toxicity
# ---------------------------------------------------------------------------
def _toxicity(mol, d: dict) -> dict:
    logp, mw, ar = d["logp"], d["mw"], d["ar"]

    # hERG — basic amine + high logP + aromatic (Sanguinetti 2006)
    has_basic_n = any(
        a.GetAtomicNum() == 7 and a.GetFormalCharge() == 0 and a.GetTotalValence() <= 3
        for a in mol.GetAtoms()
    )
    herg = _clamp(
        0.35 * _sig((logp - 3.0) / 1.0)
      + 0.35 * (1.0 if has_basic_n else 0.0)
      + 0.2  * _sig((ar - 1) / 1.0)
      + 0.1  * _sig((mw - 350) / 100)
    )

    # AMES — nitro groups, aromatic amines, flat systems
    has_nitro = mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]"))
    has_ar_nh2 = mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]c"))
    ames = _clamp(
        0.3  * _sig((ar - 2) / 1.0)
      + (0.30 if has_nitro  else 0.0)
      + (0.25 if has_ar_nh2 else 0.0)
      + 0.15 * _sig((logp - 2.0) / 2.0)
    )

    # DILI — high logP, large MW, reactive groups (Chen 2016)
    has_reactive = any([
        mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H]")),
        mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[F,Cl,Br]")),
        mol.HasSubstructMatch(Chem.MolFromSmarts("[N;X3][N;X2]")),
    ])
    dili = _clamp(
        0.3  * _sig((logp - 3.5) / 1.0)
      + 0.25 * _sig((mw - 400) / 100)
      + (0.25 if has_reactive else 0.0)
      + 0.2  * _sig((ar - 2) / 1.5)
    )

    # Skin sensitisation — Michael acceptors (Basketter 2014)
    has_michael = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[CX3][CX3](=O)"))
    skin = _clamp(
        (0.45 if has_michael else 0.0)
      + 0.25 * _sig((logp - 2.0) / 2.0)
      + (0.20 if has_nitro  else 0.0)
      + 0.10 * _sig((ar - 1) / 1.0)
    )

    def _cls_risk(v, hi=0.55, lo=0.3):
        return "High Risk" if v >= hi else ("Moderate Risk" if v >= lo else "Low Risk")

    return {
        "herg_score":  round(herg, 3), "herg_class":  _cls_risk(herg),
        "ames_score":  round(ames, 3), "ames_class":  ("Mutagenic" if ames >= 0.5 else ("Borderline" if ames >= 0.3 else "Non-mutagenic")),
        "dili_score":  round(dili, 3), "dili_class":  _cls_risk(dili),
        "skin_score":  round(skin, 3), "skin_class":  ("Sensitiser" if skin >= 0.45 else ("Possible" if skin >= 0.25 else "Non-sensitiser")),
    }


# ---------------------------------------------------------------------------
# Summary — maps onto existing AdmetProperties model (5 keys, backward compat)
# ---------------------------------------------------------------------------
def _summary(ab, di, me, ex, to) -> dict:
    absorption   = (ab["hia_score"] + ab["caco2_score"] + ab["oral_f_score"]) / 3
    distribution = (di["bbb_score"] + (1 - di["ppb_score"]) * 0.5 + 0.5) / 2
    metabolism   = 1 - (me["cyp3a4_score"] + me["cyp2d6_score"] + me["cyp2c9_score"]) / 3
    excretion    = ex["renal_clearance_score"] * 0.6 + _clamp(1 - ex["half_life_h"] / 48) * 0.4
    toxicity     = 1 - (to["herg_score"] + to["ames_score"] + to["dili_score"] + to["skin_score"]) / 4
    return {
        "absorption":   round(_clamp(absorption),   3),
        "distribution": round(_clamp(distribution), 3),
        "metabolism":   round(_clamp(metabolism),   3),
        "excretion":    round(_clamp(excretion),    3),
        "toxicity":     round(_clamp(toxicity),     3),
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def compute_admet(mol) -> dict:
    """
    Accepts an RDKit Mol, returns a flat dict with keys:
      summary, absorption, distribution, metabolism, excretion, toxicity, drug_score
    Returns {} on any error — never raises.
    """
    if mol is None:
        return {}
    try:
        d  = _desc(mol)
        ab = _absorption(d)
        di = _distribution(d)
        me = _metabolism(mol, d)
        ex = _excretion(d)
        to = _toxicity(mol, d)
        su = _summary(ab, di, me, ex, to)

        drug_score = _clamp(
            su["absorption"]   * 0.25
          + su["distribution"] * 0.15
          + su["metabolism"]   * 0.20
          + su["excretion"]    * 0.15
          + su["toxicity"]     * 0.25
        )

        return {
            "summary":      su,
            "absorption":   ab,
            "distribution": di,
            "metabolism":   me,
            "excretion":    ex,
            "toxicity":     to,
            "drug_score":   round(drug_score, 3),
        }
    except Exception as e:
        print("[ADMET] Error: " + str(e))
        return {}

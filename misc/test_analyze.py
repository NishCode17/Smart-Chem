import requests

BASE = "http://localhost:8000"

molecules = {
    "Aspirin":    "CC(=O)Oc1ccccc1C(=O)O",
    "Benzene":    "c1ccccc1",
    "Bad SMILES": "NOT_A_SMILES",
}

for name, smi in molecules.items():
    r = requests.post(f"{BASE}/utils/analyze", json={"smiles": smi})
    if r.status_code == 200:
        d = r.json()
        print(f"[{name}] status={r.status_code}")
        print(f"  QED={d.get('qed')}  LogP={d.get('logp')}  valid={d.get('valid')}")
        print(f"  admet keys: {list(d.get('admet', {}).keys())}")
        print(f"  admet_props keys: {list(d.get('admet_props', {}).keys())}")
        print(f"  tox_alerts: {d.get('tox_alerts', {})}")
    else:
        print(f"[{name}] status={r.status_code}  error={r.json()}")
    print()

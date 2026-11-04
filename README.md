# Expanded CRI Study Repository

Questo repository contiene il framework di **Cellular Response Index (CRI)** esteso per supportare tumori analoghi a quelli EGFR‑mutati. Include:

- Moduli di modellazione (PK/PD, reaction‑diffusion, multiscale, agent‑based)
- Pipeline di test automatizzati (unit, integration, benchmark, stress)
- Documentazione e piano di progetto (vedi `PLAN.md`)
- CI configurata con GitHub Actions

## Struttura del progetto
```
expanded_cri_study/
├─ src/                 # Codice sorgente dei moduli
│   ├─ __init__.py
│   ├─ pk_pd.py
│   ├─ reaction_diffusion.py
│   └─ multiscale.py
├─ tests/               # Test suite
│   ├─ __init__.py
│   └─ test_pk_pd.py
├─ data/                # Dati di esempio (git‑ignored)
├─ PLAN.md              # Piano dettagliato di lavoro
├─ .gitignore
└─ .github/workflows/ci.yml  # CI workflow
```
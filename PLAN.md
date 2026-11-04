## Piano di lavoro per lo studio CRI ampliato

### 1. Obiettivi
- Creare un framework *generico* per calcolare il **Cellular Response Index (CRI)** su tumori con mutazioni EGFR, KRAS e altri tumori solidi analoghi.
- Garantire **verifica** e **contro‑verifica** attraverso una suite completa di test (unit, integration, benchmark, stress).
- Produrre **documentazione** reproduci­bile e preparare materiale per una eventuale sottomissione a *Nature Cancer*.

### 2. Scoping Matrix (tumori ↔ modelli)
| Tumore target | Modello principale | Riferimenti (PMID) |
|---|---|---|
| EGFR‑mutated NSCLC | PK/PD + Reaction‑Diffusion + Multiscale | 42977708, 42969431, 42966224 |
| KRAS‑mutated lung cancer | Agent‑Based + Network meta‑analysis | 42977708, 42962543 |
| Solid tumor generic (breast, gastric) | Multiscale + Biophysical modeling | 42962543, 42946142 |
| Metastatic brain tumour | Diffusion‑reaction + PK/PD | 42977708, 42936185 |

### 3. Moduli da sviluppare (src/)
- `pk_pd.py` – fitting PK/PD con **Stan** / **PyMC**.
- `reaction_diffusion.py` – solver PDE basato su **FEniCS** o **COMSOL** (wrapper python).
- `multiscale.py` – integrazione tra modello cellulare (agent‑based, **Morpheus**) e tessuto (PDE).
- `utils.py` – funzioni comuni (loader dati, metriche, visualizzazioni).

### 4. Suite di test (tests/)
| Tipo | Descrizione | File di test |
|------|-------------|--------------|
| Unit | Verifica funzioni base (ODE, PDE, regressioni) | `test_pk_pd.py` |
| Integration | Pipeline completa su dataset di esempio (TCGA lung) | `test_integration.py` |
| Benchmark | Confronto con risultati pubblicati (PMID) | `test_benchmark.py` |
| Stress | Simulazioni su 10⁶ cellule usando **dask** | `test_stress.py` |
| Reproducibility | Esecuzione su Docker/Conda env | `test_repro.py` |

### 5. CI / Workflow (GitHub Actions)
- **ci.yml** esegue tutti i test con matrix su Python 3.11, installa dipendenze (`requirements.txt`).
- Genera artefatti di coverage (`coverage.xml`) e report di benchmark (`benchmark.json`).
- Deploy automatico di **documentation** su *Read the Docs*.

### 6. Data Management
- Usa **DVC** per versionare dataset grezzi (`data/raw/`) e pre‑processati.
- Log dei modelli con **MLflow** (parametri, metriche, artefatti).

### 7. Milestones
| Milestone | Scadenza | Deliverable |
|----------|----------|-------------|
| Setup repo + CI | 2026‑11‑07 | Repo iniziale con README, .gitignore, ci.yml.
| Implementazione PK/PD | 2026‑11‑14 | `src/pk_pd.py` + unit test.
| Reaction‑Diffusion module | 2026‑11‑21 | `src/reaction_diffusion.py` + integrazione test.
| Multiscale module | 2026‑11‑28 | `src/multiscale.py` + integrazione.
| Benchmark & Stress test | 2026‑12‑05 | Suite completa, report benchmark.
| Draft manuscript & docs | 2026‑12‑12 | `docs/` popolati, draft `manuscript.md`.
| Review & final submission | 2026‑12‑19 | Pacchetto pronto per *Nature Cancer*.

### 8. Responsabilità (esempio)
- **Lead Scientist** – definizione requisiti, revisione scientifica.
- **Data Engineer** – gestione DVC, pipeline dati.
- **Software Engineer** – sviluppo moduli, CI.
- **Bio‑statistician** – definizione metriche, analisi benchmark.
- **Documentation Lead** – scrittura README, tutorials.

### 9. Come contribuire
1. Fork del repository.
2. Creare branch `feature/<nome>`.
3. Implementare/modificare codice.
4. Eseguire `git push` e aprire una Pull Request.
5. CI verifierà automaticamente tutti i test.

---
*Questo piano è salvato come artefatto `PLAN.md` ed è parte integrante del repository.*
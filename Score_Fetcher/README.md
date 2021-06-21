# Score Fetcher Utility

## `getscores.py`

Fetches various metrics of synthetic complexity/feasibility. Defines
a ScoreFetcher class which when constructed can be passed any SMILES string
to obtain one of the following metrics:

- SCScore: Synthetic Complexity Score
    - Developed by the Coley Group, SCScore is a score on a scale from 1 to 5
    determined by a neural network model trained on millions of reactions from Reaxys
    to predict the level of complexity of a molecule.
- SAScore: Synthetic Accessibility Score
    - A heuristic score between 1 and 10 that approximates the accessibility of a compound
    using a fragment-based approach.
- RAScore: Retrosynthetic Accessibility Score
    - A score between 0 and 1 learned from the CASP tool AiZynthfinder that predicts
    if a compound would have a synthetic route (1) or not (0).

Example usage:
```python
    model = ScoreFetcher()
    # List of SMILES strings
    smis = ['CC1=C(COC2=CC(OC)=C(CN3CCCC[C@H]3C(O)=O)C(OC)=C2)C=CC=C1C4=CC=CC=C4',
            'COC1=CC(OCC2=CC=CC(=C2C)C2=CC=C3OCCOC3=C2)=CC(OC)=C1CN[C@H](CO)C(O)=O',
            'COC1=CC(OCC2=CC=CC(=C2C)C2=CC=CC=C2)=CC(OC)=C1CN[C@H](CC(F)(F)F)C1=CC=CC=C1',
            'O=C1CCC2C(C[C@@H](C)[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4(C#C)O)=C1',
            'COC1=CC(OCC2=C(C)C(=CC=C2)C2=CC=CC=C2)=CC(OC)=C1CN1CCC[C@@H](C1)C(O)=O',
            'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5']
    for smi in smis:
        model.get_SCScore(smi)
        model.get_SAScore(smi)
        model.get_RAScore(smi)
```
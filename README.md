# Slater-Condon-Rules

This is simple python code for applying Slater-Condon rules on given determinats, here how to use it:

**Input**

You’ll be asked for:
Input file (.txt or .csv) containing determinants

And the code will give CSV output with each row = ⟨Dᵢ|H|Dⱼ⟩, including the excitation degree, phase, terms, and mapping.

**Determinant file format**

Must start with a label: D_1, D2, etc.

Inside parentheses: comma-separated spin-orbitals.

Each spin-orbital is written as f_{name} s, where:

name = orbital label (e.g. V, C, O1)

s = a (alpha) or b (beta

Here is and example:
D_1 = (f_{V} b, f_{C} a, f_{O1} a, f_{O2} b)
D_2 = (f_{C} a, f_{C} b, f_{O1} a, f_{O1} b)
D_3 = (f_{C} a, f_{C} b, f_{O2} b, f_{O2} a)

## Usage

```bash
python3 sc_applier.py
```

**Ref:**
Scemama, A., & Giner, E. (2013, November 25). An efficient implementation of Slater-Condon rules. arXiv.org. https://arxiv.org/abs/1311.6244

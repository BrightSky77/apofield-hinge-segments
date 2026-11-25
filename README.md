# Hinge-based segment boundaries for ApoField

**Why ApoField fails on large conformational changes, and fixing where it cuts the chain**

Course research proposal building on ApoField, a holo→apo structure prediction project in Prof. Ruishan Liu's lab (USC) that I am working on with PhD student Jaewon Choi.

## 1. Summary

ApoField predicts the unbound (apo) structure of a protein chain from its bound (holo) structure in about 0.27 s per chain. It moves the chain as 12 rigid segments. Overall it is competitive with AlphaFold 3, but on **flex** chains (holo–apo Cα RMSD > 2 Å) it trails by about 0.7 Å (2.72 Å against 2.03 Å).

ApoField sets its segment boundaries from a single chain-averaged vector, starting from equal-length blocks. It may therefore miss where the protein actually bends.

This proposal does two things.

1. Measures how much of the flex error comes from **where the segments are cut**.
2. Replaces the boundaries with **hinge positions found from an elastic network model (ANM) of the holo structure alone**.

To check whether this direction is worth pursuing, I ran a quick preliminary test on a public antibody–antigen dataset. I took 31 proteins whose bound and unbound structures differ by more than 2 Å and cut each into 12 rigid pieces in three ways. After moving each piece to best match the true apo structure, the remaining Cα RMSD was as follows.

| Boundaries | Residual |
|---|---|
| Equal-length blocks | 1.35 Å |
| Holo-only ANM boundaries | 0.78 Å |
| Oracle-optimal boundaries | 0.67 Å |

## 2. Background

### Why holo→apo prediction matters

Proteins often change shape when they bind a partner. Many downstream tasks need the **unbound (apo)** shape, because that is the state a protein is in before binding.

- **Docking and binding prediction.** Real use cases start from unbound proteins. Methods tested on bound structures look far better than they are in practice, so realistic benchmarks and training sets need apo inputs.
- **Drug discovery.** Comparing bound and unbound shapes shows how a binding site opens, closes or moves. This matters for finding binding pockets and designing binders.
- **Data scale.** Experimental apo structures are rare, and most proteins in the PDB have been solved only in complexes. Large datasets such as PINDER therefore need apo structures that are predicted, not measured.

General structure predictors such as AlphaFold 3 can produce such structures, but they take minutes per chain once the MSA search is included, and they tend to return the bound form. A fast model that turns an existing holo structure into its apo form makes apo structures available at dataset scale.

### How ApoField sets segment boundaries

ApoField is led by PhD student Jaewon Choi and advised by Prof. Ruishan Liu. It predicts apo structures for chains in PINDER complexes and proceeds in four steps.

1. It starts from a per-residue blend of the holo and AlphaFold DB structures.
2. It assigns residues softly to K = 12 segments.
3. It predicts one rigid motion per segment, then adds small per-residue corrections.
4. At inference it shrinks small predicted motions back toward the starting structure.

Membership of residue *i* in segment *k* is a Gaussian over sequence position.

```
mu_k    = (k + 0.5) / K * L + a_k * L / K
sigma_k = (L / K) * exp(b_k)
G_ik    = softmax_k( 20 * [ log pi_k + log N(i; mu_k, sigma_k) ] )
```

The parameters `a`, `b` and `pi` come from **one vector, the mean of the chain's residue embeddings**. Their output layer is zero-initialised, so training starts from exactly equal-length blocks. Residue-level information, including the ANM features the model already takes as input, influences the boundaries only through this average.

### DEV 200 results (global Cα Kabsch RMSD, median, Å)

| Method | All (200) | Flex (71) | Rigid (44) | Time per chain |
|---|---|---|---|---|
| Copy holo | 1.60 | 3.01 | 0.76 | 0 |
| AlphaFold DB | 1.18 | 2.17 | 0.59 | lookup |
| ApoField | 1.37 | 2.72 | 0.66 | 0.27 s |
| AlphaFold 3 (+MSA) | 1.33 | 2.03 | 0.70 | 49 s + MSA |

Flex means holo–apo > 2 Å; rigid means ≤ 1 Å.

## 3. Related work and gap

| Line of work | Examples | Relevance |
|---|---|---|
| Sequence-based structure prediction | AlphaFold 3, Boltz-2, Chai-1; [Saldaño et al. 2022](https://www.biorxiv.org/content/10.1101/2021.10.27.466189) | Slow. AlphaFold2 tends to predict the holo form, and its accuracy drops as conformational diversity grows |
| Fast state-conditioned prediction | [Ibex](https://arxiv.org/abs/2507.09054) | Predicts antibody apo and holo states, but underestimates very large changes |
| Quasi-rigid domains | [Petrone & Pande 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1367309), [SPECTRUS](https://iris.sissa.it/handle/20.500.11767/17197), [PiSQRD](https://people.sissa.it/~michelet/publications/abstracts/64.html) | Large changes are well approximated by a few quasi-rigid domains, which can be found from elastic-network fluctuations |
| Contiguous optimal partitioning | [TLSMD](https://doi.org/10.1107/S0021889805038987) | Optimally splits a chain into contiguous segments using B-factors. My boundary step is the same idea, with ANM in place of B-factors |
| Single-structure hinge prediction | [FlexOracle](https://doi.org/10.1186/1471-2105-8-215), [StoneHinge](https://pmc.ncbi.nlm.nih.gov/articles/PMC2708048), [HingeMaster](https://kuhnlab.natsci.msu.edu/publication_papers/pdf/Flores_et_al.Proteins08.pdf) | Predict hinge locations from one structure |

**Gap.** Hinge and domain partitioning methods and learned holo→apo models have developed separately. I found no study that uses structure-derived boundaries as the motion units of a learned holo→apo model. Learned models such as ApoField infer their boundaries from chain-level features only.

## 4. Research questions and hypotheses

- **RQ1.** How much of ApoField's flex error comes from where the segments are cut?
- **RQ2.** Do hinge-based boundaries reduce it?

I split the flex error into three parts. All three use the true apo and serve only for diagnosis.

| Part | Meaning | How it is measured |
|---|---|---|
| (1) Representation limit | What 12 rigid segments cannot express even with the best cuts | Residual after fitting each segment rigidly to apo, using oracle-optimal contiguous boundaries |
| (2) Boundary error | Error added by cutting in the wrong places | Same fit with ApoField's learned boundaries, minus (1) |
| (3) Remaining prediction error | Error from the predicted segment motions themselves | ApoField's actual error, minus the fit with its learned boundaries |

**Premise.** With the best cuts, 12 rigid segments can express most of a flex change. In the preliminary test, oracle boundaries reduced 3.45 Å to 0.67 Å. Part (1) checks this on DEV. If (1) turns out large, better boundaries cannot help.

**Hypotheses**

- **H1 (the current boundaries are wrong).** ApoField's learned boundaries are close to equal blocks, and the boundary error (2) is a substantial part of its flex error.
- **H2 (hinge boundaries fix them).** ApoField with holo-ANM boundaries has lower flex error than the current model, and its rigid error is no worse.

## 5. Methods

### Data

- **Evaluation.** The 71 flex chains of PINDER DEV.
- **Training.** Flex chains in PINDER train.
- **Preliminary test.** The public [antibody–antigen benchmark](https://github.com/piercelab/antibody_benchmark).

I use Cα only. Apo structures and predictions are Kabsch-aligned to holo.

### Diagnosis (RQ1)

- **(1)** Find oracle-optimal contiguous 12-segment boundaries by dynamic programming. Candidate cuts lie every 5 residues and at chain breaks.
- **(2)** Take argmax_k G_ik as ApoField's boundaries and repeat the fit.
- **(3)** Compare ApoField's actual RMSD with the fit using its learned boundaries.
- Results are stratified by chain length, single- vs multi-chain partner, and size of change.

### Hinge-based boundaries (RQ2)

**Boundaries from holo only.** Build an ANM on the holo Cα atoms (cutoff 13 Å) and take the 20 lowest-frequency modes, each scaled by 1/sqrt(eigenvalue). For each candidate contiguous segment, compute the part of the mode displacements that one rigid-body motion of that segment cannot explain. Dynamic programming then picks the 12 segments that move most rigidly; the boundaries between them are the predicted hinges.

**(a) Training-free comparison.** On DEV flex and PINDER train flex, compare the 12-segment upper bound under four boundary sets, namely equal blocks, ApoField's learned boundaries, ANM boundaries and oracle boundaries.

**(b) Retraining variants.**

| Variant | Change | Question |
|---|---|---|
| V1 init | Initialise μ_k and σ_k from the ANM segments instead of equal blocks, then train | Does a better starting point help? |
| V2 fixed | Fix G to the ANM boundaries | Can boundary learning be removed? |
| V3 per-residue | Add a per-residue boundary head, supervised by oracle boundaries from train, with ANM boundaries as an input feature | Can learning beat ANM? |

All other settings stay as in the current model. Success means lower flex error than the current ApoField, with rigid error no worse.

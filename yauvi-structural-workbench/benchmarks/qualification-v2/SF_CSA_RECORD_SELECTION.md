# sf-csa record selection — 16 judgments, strata derived from SCOP

Every stratum below is assigned by `sf_csa_stratum_from_scop.py`, which reads
SCOP fold and superfamily from the PDBe SIFTS mapping and applies the definition
in Russell, Saqi, Sayle, Bates and Sternberg (J Mol Biol, 1997):

    same superfamily              -> homologous_superfamily
    same fold, other superfamily  -> fold_analogy
    different fold                -> unrelated

No stratum here was assigned from judgement or recall, and none was assigned from
a TM-score. That matters: using structural similarity to set the stratum and then
testing whether the module agrees is circular, and it is how a false-positive
gate ends up passing for the wrong reason. Raw tool output is preserved in
`SF_CSA_SELECTION_EVIDENCE.txt`.

Reproduce with:

    python3 sf_csa_stratum_from_scop.py 1fdn 1dur 1ris 1csp 1mjc 1e6a \
        1fi2 1juh 1os7 2hmz 2mhr 1cgn

## Queries are AlphaFold models, not crystal structures

Settled by measurement, 2026-09-01. sf-csa refuses any query whose PDB sequence
does not exactly match its FASTA. Screening experimental entries against that
rule, only 9 of 13 candidates passed, and the failures were systematic rather
than unlucky: unresolved residues. Every larger double-stranded β-helix screened
— twelve of them — failed, because multi-superfamily folds are mostly enzymes
and enzymes have disordered loops. The requirement quietly restricts the panel to
small, fully-ordered proteins, which is not a scientific criterion.

Predicted models have no unresolved residues, and all four entries that failed as
crystal structures are sequence-exact as AlphaFold models. This is not a
workaround: sf-csa's own declared limitation says query structures are "exact
predicted monomers, not experimental assemblies or active poses", and the panel
scope draft repeats it. Using crystal structures was fighting the module's
contract; open item 3 of the earlier draft is resolved in favour of models.

**Concordance is required, not assumed.** An AlphaFold model covers the whole
UniProt protein while SCOP classified the crystallised construct, and where those
differ the model is a different molecule from the one whose fold was classified.
Two examples from this screen: 1SNC is a 149-residue construct against a
231-residue model, and 1URN is a single RRM domain against a 282-residue
multi-domain protein. Both are excluded. Every entry below was checked and the
model differs from the deposited chain by at most one residue, which absorbs an
initiator methionine and nothing larger.

Model version is read from the AlphaFold API per accession rather than hardcoded.
It is v6 today; a hardcoded v4 URL returned 404 and briefly looked like the
service being unreachable.

## The four families

Each family is a SCOP fold that contains more than one superfamily — the
precondition for a `fold_analogy` record to exist at all.

| # | family | SCOP fold | superfamilies used |
|---|---|---|---|
| 1 | ferredoxin-like | `d.58` | `d.58.1`, `d.58.7` |
| 2 | OB-fold | `b.40` | `b.40.4`, `b.40.1` |
| 3 | double-stranded β-helix | `b.82` | `b.82.1`, `b.82.2` |
| 4 | four-helical up-and-down bundle | `a.24` | `a.24.4`, `a.24.3` |

Ferredoxin-like and double-stranded β-helix are two of the three superfolds
Russell and colleagues (J Mol Biol, 1998) identify as carrying analogous members.
The TIM barrel is deliberately absent: the same paper places its supersite among
*homologous* proteins, and the sequence evidence for shared origin across its
superfamilies is why it was withdrawn from the earlier plan.

## A fourth gate: the organism must have a reference proteome

sf-csa needs a source proteome per query for the reciprocal-best-hit leg and
`run_pipeline` raises without one, so this is a hard criterion. It removed three
entries that had passed every other check:

- **2HMZ and 2MHR** (hemerythrin, myohemerythrin) come from sipunculid worms,
  *Themiste dyscrita* and *Themiste hennahi*. Neither has a UniProt proteome.
  They were the entire four-helical family, so the family was rebuilt.
- **1JUH** (quercetin 2,3-dioxygenase, *Aspergillus japonicus*) — no proteome.
  It was family 3's homolog partner.

These were ideal on fold grounds and useless in practice. Discovering that after
freezing the records rather than before would have meant a panel that could not
execute at all.

Selection therefore passes four gates, in this order: one unambiguous SCOP domain;
stratum derived from SCOP; an AlphaFold model concordant with the deposited
construct; and a UniProt reference proteome for the organism.

## The 16 judgments

Four families, each a SCOP fold containing more than one superfamily. Twelve
distinct proteins; the query structure is the AlphaFold model of the accession
named, and the sequence leg searches the organism's reference proteome.

| # | family | stratum | query → target | basis |
|---|---|---|---|---|
| 1 | ferredoxin-like | `exact` | 1FDN → 1FDN | accession identity |
| 2 | ferredoxin-like | `homologous_superfamily` | 1FDN → 1DUR | both `d.58.1.1`, 4Fe-4S ferredoxins |
| 3 | ferredoxin-like | `fold_analogy` | 1FDN → 1RIS | fold `d.58`; `d.58.1` vs `d.58.14` ribosomal S6 |
| 4 | ferredoxin-like | `unrelated` | 1FDN → 1CSP | `d.58` vs `b.40` |
| 5 | OB-fold | `exact` | 1CSP → 1CSP | accession identity |
| 6 | OB-fold | `homologous_superfamily` | 1CSP → 1MJC | both `b.40.4.5`, nucleic-acid-binding |
| 7 | OB-fold | `fold_analogy` | 1CSP → 1E6A | fold `b.40`; `b.40.4` vs `b.40.5` inorganic pyrophosphatase |
| 8 | OB-fold | `unrelated` | 1CSP → 1FI2 | `b.40` vs `b.82` |
| 9 | ds β-helix | `exact` | 1FI2 → 1FI2 | accession identity |
| 10 | ds β-helix | `homologous_superfamily` | 1FI2 → 1DZR | superfamily `b.82.1`; families `.2` vs `.1` |
| 11 | ds β-helix | `fold_analogy` | 1FI2 → 1OS7 | fold `b.82`; `b.82.1` vs `b.82.2` clavaminate synthase-like |
| 12 | ds β-helix | `unrelated` | 1FI2 → 1CGN | `b.82` vs `a.24` |
| 13 | four-helical | `exact` | 1CGN → 1CGN | accession identity |
| 14 | four-helical | `homologous_superfamily` | 1CGN → 1CPQ | both `a.24.3.2`, cytochromes, different organisms |
| 15 | four-helical | `fold_analogy` | 1CGN → 1C02 | fold `a.24`; `a.24.3` vs `a.24.10` HPT domain |
| 16 | four-helical | `unrelated` | 1CGN → 1FDN | `a.24` vs `d.58` |

### The twelve queries — all four gates passed

| entry | SCOP | UniProt | model | proteome | protein |
|---|---|---|---|---|---|
| 1FDN | `d.58.1.1` | P00198 | 56 | UP000006094 | ferredoxin |
| 1DUR | `d.58.1.1` | P00193 | 54 | UP000603820 | ferredoxin |
| 1RIS | `d.58.14.1` | P23370 | 101 | UP000217909 | ribosomal protein S6 |
| 1CSP | `b.40.4.5` | P32081 | 67 | UP000001570 | cold shock protein B |
| 1MJC | `b.40.4.5` | P0A9X9 | 70 | UP000000625 | major cold shock protein |
| 1E6A | `b.40.5.1` | P00817 | 287 | UP000002311 | inorganic pyrophosphatase |
| 1FI2 | `b.82.1.2` | P45850 | 201 | UP001057469 | germin / oxalate oxidase |
| 1DZR | `b.82.1.1` | P26394 | 183 | UP000001014 | dTDP-sugar epimerase |
| 1OS7 | `b.82.2.5` | P37610 | 283 | UP000000625 | taurine dioxygenase |
| 1CGN | `a.24.3.2` | P00138 | 127 | UP000595916 | cytochrome c3 |
| 1CPQ | `a.24.3.2` | P00147 | 129 | UP000310597 | cytochrome c3 |
| 1C02 | `a.24.10.2` | Q07688 | 167 | UP000002311 | HPT domain, YPD1 |

All twelve acquired and verified sequence-exact against their UniProt sequence;
24 artifacts with SHA-256 recorded for `SOURCE_LOCK.json`.

Screened and rejected, with the reason kept rather than dropped: 2ACY (model
+3), 1URN (+185, one domain of a multi-domain protein), 1SNC (+82, mature
construct against full precursor), 256B (+22, signal peptide), 2HMZ / 2MHR /
1JUH (no reference proteome), 1BKB (two SCOP domains).

## Four things to settle before these are frozen

**1. The periodontal table control cannot be a family.** BamA is not classified
in SCOP: PDBe serves SCOP 1.75, which predates 5D0O and 5AYW, and both return
unclassified. So `omp85_bama` cannot take a family slot without its strata being
assigned on a different basis from the other twelve, which would quietly mix two
definitions of "homologous" inside one panel.

Better, and it costs nothing: make the table control a **control**, not a family.
Controls live in the panel's separate `controls` key and are counted
independently, so the four families stay uniformly SCOP-derived while the control
still does its job — running one `omp85_bama` comparison under the defaults and
under the override and requiring the same judgment, which is the evidence the
override changed the biology it was meant to change and nothing else. Its stratum
basis is recorded as non-SCOP.

This is a change to the family decision of 2026-08-31 and needs Yuvraj's
agreement before the records are frozen. The three classic families become four.

**2. Partly resolved — SCOP classifies domains; sf-csa compares whole structures.** `run_pipeline`
refuses a query whose PDB sequence does not exactly match its FASTA, so every
entry above must be single-domain and sequence-clean. The tool already flags
multi-fold entries — 1BKB was flagged and dropped during selection for exactly
this reason.

Checked, because the tool's weakest step is picking `domains[0]` when an entry has
several: **all twelve selected entries carry exactly one SCOP domain**, with one
distinct fold and one distinct superfamily each. The heuristic is therefore not
hiding a second assignment anywhere in this selection.

Single-domain is still a weaker check than sequence-exact, so each of the twelve
needs verifying against its FASTA before it is frozen. Substitutes are available
from the 18 candidate analogy pairs.

**3. Resolved — queries are AlphaFold models.** See the section above. The
scope note and the selection now agree, and the sequence-exact requirement is
satisfied by construction rather than by restricting the panel to small proteins.

**4. `exact` is a self-match.** Records 1, 5, 9 and 13 compare a query to itself,
which `classify_hit` settles by accession before any structural reasoning. That is
a genuine control on the identity path, but it is the cheapest of the four strata
and should not be read as evidence about fold comparison at all.

## What this does not establish

The strata are as good as SCOP's superfamily boundaries, no better. A pair
recorded here as `fold_analogy` means "same fold, different superfamily in SCOP
1.75" — not "proven independent origin". Wright (Genome Biol Evol, 2025) is the
caution worth carrying: across Foldseek clusters only about 1% of strong matches
at TM ≥ 0.5 lack sequence-level homology support, and structural resemblance
alone does not establish homology in either direction. Each record carries the
classification it came from, so a reviewer can disagree with a specific call
without having to re-derive the whole panel.

---

## Resolved — 2026-09-01, Yuvraj

**Item 1: the periodontal table control is a control, not a family.** Agreed as
proposed. The four families stay uniformly SCOP-derived; `omp85_bama` moves to
the panel's `controls` key, run once under the default tables and once under the
override, requiring the same judgment. Its stratum basis is recorded as non-SCOP.
This supersedes the family decision of 2026-08-31.

Items 2 and 4 stand as written. Item 3 was already resolved in favour of
AlphaFold models.

The twelve queries still need verifying sequence-exact against their FASTA before
freezing (item 2), which is a precondition for curation, not for this decision.

## Item 2 discharged — 2026-09-01

All twelve queries verified sequence-exact against their UniProt FASTA, which is
the precondition `run_pipeline` enforces and refuses without.

    1FDN P00198  56   1DUR P00193  54   1RIS P23370 101   1CSP P32081  67
    1MJC P0A9X9  70   1E6A P00817 287   1FI2 P45850 201   1DZR P26394 183
    1OS7 P37610 283   1CGN P00138 127   1CPQ P00147 129   1C02 Q07688 167

Twelve of twelve exact, no missing artifacts. The AlphaFold v6 models and the
FASTA files are already in `sources/alphafold` and `sources/uniprot`, so the
records can be frozen without acquiring anything further — which matters, because
`execution_policy.network_access` is `forbidden` and nothing may be fetched at run
time.

Item 1 (the periodontal table control) was resolved earlier today; item 3 was
already resolved in favour of AlphaFold models; item 4 is a statement about the
`exact` stratum, not an action. **All four items in "Four things to settle before
these are frozen" are now settled.**

---

## Reselected — 2026-09-02, collection 2.7

The twelve below replace the twelve above. Everything in the earlier sections
describing the old set is superseded; it is left in place because the reason for
the change is the interesting part.

**What was wrong.** The fourth gate was written as "the organism must have a
reference proteome" and was applied by organism. It verified that a proteome
*exists*, never that the query is *in* it. Five entries were in no reference
proteome at all — P00193, P23370, P45850, P00138, P00147 — and P26394 was pointed
at UP000002695, which loads and contains zero occurrences of it. The sequence leg
silently vanished for thirteen of sixteen records, including two `exact`
self-matches, which is what exposed it. `sf_csa_screen_candidates.py` has been
corrected to read the entry's own Proteomes cross-reference.

**The fifth gate**, now stated: *the accession must appear in the proteome it
declares*. Verified 12/12 before the table was written.

| family | SCOP fold | exact | homologous_superfamily | fold_analogy | unrelated |
|---|---|---|---|---|---|
| ferredoxin-like | `d.58` | 1FDN | 1BLU `d.58.1.1` | 1APS `d.58.10.1` | 1CSP |
| OB-fold | `b.40` | 1CSP | 1MJC `b.40.4.5` | 1E6A `b.40.5.1` | 1DZR |
| ds β-helix | `b.82` | 1DZR | 1J1L `b.82.1.12` | 1OS7 `b.82.2.5` | 1C02 |
| four-helical | `a.24` | 1C02 | 2OOC `a.24.10.6` | 1JOG `a.24.16.2` | 1FDN |

    P00198 1FDN  UP000006094      P00208 1BLU  UP000001441      P00818 1APS  UP000002281
    P32081 1CSP  UP000001570      P0A9X9 1MJC  UP000000625      P00817 1E6A  UP000002311
    P26394 1DZR  UP000001014      O00625 1J1L  UP000005640      P37610 1OS7  UP000000625
    Q07688 1C02  UP000002311      Q9A980 2OOC  UP000001816      P43934 1JOG  UP000000579

Six of the twelve are retained. Retired coordinates and proteomes were **moved to
`sources/retired-2026-09-02/`, not deleted**, so the previous state is
recoverable.

**Two choices worth seeing.** 1J1L is human, which adds 145k sequences to the
search universe; it was chosen over 1EP0 because 1EP0 sits in the same SCOP
*family* as 1DZR rather than a different one, making it a closer homolog and an
easier recall test. And the four-helical family lost its cytochromes entirely —
every `a.24.3` entry is P00138 or P00147, both proteome-less — so it is now
anchored on the HPT domain, with the nucleotidyltransferase superfamily supplying
the analogy.

**Result:** 16/16 judgments located, 16/16 passing the frozen expectation, zero
promotions against a bound of zero.

**The limitation this created, stated rather than repaired.** No
`homologous_superfamily` record in this selection carries a sequence row: 1FDN→1BLU
is above both thresholds and DIAMOND still does not seed it, and the other two are
genuinely below the coverage minimum. That stratum is therefore recovered on
structural evidence alone. Repairing it would mean selecting homologs for DIAMOND
detectability, which is choosing records to make a gate pass. See Finding 6.

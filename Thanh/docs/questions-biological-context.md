# The biological meaning of the ten Track 1 questions

Prepared November 10, 2026. Numbering follows the organizer's [starter questions](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/starter-questions.md). The headings below paraphrase the prompts; explanations, suggested readings, and prioritization are our interpretation, not additional organizer instructions.

In this context, a transcriptomic signature concerns measured RNA abundance. It does not by itself establish a DNA mutation, inherited risk variant, protein activity, or disease cause.

## 1. Characterize the disease-associated expression pattern

**What it asks:** What reproducible biological processes distinguish affected tissue from the selected control tissue?

**Why it matters:** Before proposing a mechanism, establish the phenomenon to explain. A coherent change in a cellular process is usually easier to interpret than an isolated list of genes.

**Biology to learn:** Stress responses, metabolism, angiogenesis, and the functions of the sampled tissue. Use the source paper to understand why metabolic variation was investigated. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/)

**Evidence boundary:** Describe an association with a stated comparison. It cannot identify which change initiated the disease.

## 2. Examine variation among affected pregnancies

**What it asks:** Do affected samples share one pattern, or do some exhibit different cellular programs?

**Why it matters:** A common clinical diagnosis can include different biological states. An average may conceal a meaningful subgroup.

**Biology to learn:** Phenotype versus molecular state, metabolic adaptation, and relationships among disease severity, FGR, and tissue stress. The source study found a metabolic expression pattern in only a subset of EOPE specimens. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/)

**Evidence boundary:** With a small cohort, apparent subgroups are exploratory. Check whether timing or clinical composition provides a simpler explanation before naming a new subtype.

## 3. Localize changes to particular cell populations

**What it asks:** Which cell populations show the clearest within-population disease-associated differences?

**Why it matters:** A tissue signal becomes more interpretable when the responsible cell population can be identified.

**Biology to learn:** Trophoblast differentiation, stromal and endothelial functions, and fetal versus maternal macrophages. Admati's study provides a directly relevant cell-resolved disease example. [Admati et al.](https://pubmed.ncbi.nlm.nih.gov/37572658/)

**Evidence boundary:** The optional donor-level pseudobulk extension can address within-cell-type differences more directly. Merely mapping a bulk gene to a reference cell type does not show that the population changed in the disease samples.

## 4. Distinguish cell abundance from altered cell state

**What it asks:** Is a tissue-level change compatible with different numbers of cells, different expression within those cells, or both?

**Why it matters:** These explanations imply different mechanisms. More immune-cell transcripts could mean recruitment, activation, or both.

**Biology to learn:** Cell identity markers versus activation programs, differentiation, recruitment, and tissue architecture. The distinction between bulk measurements and cell-resolved measurements is central. [Organizer workflow](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/3-workflow.md)

**Evidence boundary:** Marker expression alone is not a direct cell count. Cell-resolved data, histology, or spatial measurements can help separate the explanations.

## 5. Separate fetal growth effects from maternal disease

**What it asks:** Which placental programs track FGR, preeclampsia, or their combination?

**Why it matters:** Impaired fetal growth and maternal hypertension overlap clinically but are not interchangeable outcomes.

**Biology to learn:** Placental exchange, tissue metabolism, vascular supply, and the difference between FGR and being constitutionally small. The source cohort includes FGR without preeclampsia, making the study design particularly relevant. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/)

**Evidence boundary:** The clinical groups are small and differ in other ways. An expression association with FGR does not establish deficient nutrient transfer as its cause.

## 6. Separate disease effects from developmental timing

**What it asks:** Could a proposed disease signature instead reflect pregnancy stage at delivery?

**Why it matters:** The placenta changes during normal development, and disease can cause earlier delivery. A comparison with term controls therefore mixes important biological differences.

**Biology to learn:** Placental maturation, clinical onset versus delivery timing, and the relation between development and disease. Healthy atlases across gestation provide context. [Wang et al. 2026](https://www.nature.com/articles/s41586-026-10316-x)

**Evidence boundary:** Gestational-age adjustment is not automatically sufficient: the groups need useful overlap, and gestational age can also lie on a pathway through medically indicated delivery. State what comparison you intend to estimate. [Organizer data guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/2-the-data.md)

## 7. Compare onset-associated disease patterns

**What it asks:** Are early- and late-onset presentations associated with different cellular and molecular patterns?

**Why it matters:** Similar maternal symptoms need not imply identical placental mechanisms.

**Biology to learn:** Clinical timing, placental stress, and heterogeneity. Admati and colleagues observed pronounced early-onset changes and much smaller late-onset differences in their cohort. This is evidence to evaluate, not a universal rule to impose on new data. [Admati et al.](https://pubmed.ncbi.nlm.nih.gov/37572658/)

**Evidence boundary:** The core dataset contains EOPE; the onset comparison needs the optional extension. The organizer notes only three early control donors versus ten early cases. More cells do not create more independent pregnancies. [Starter-question caveats](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/starter-questions.md)

## 8. Propose altered signaling relationships

**What it asks:** Which sending and receiving cell populations, and which molecular interactions, are plausible candidates for altered communication?

**Why it matters:** This connects expression patterns to a specific mechanistic hypothesis that can guide a future experiment.

**Biology to learn:** Ligands, receptors, molecular complexes, secreted versus membrane-bound forms, downstream responses, and anatomical proximity. The normal interface atlas demonstrates a framework for inferring candidate interactions. [Vento-Tormo et al.](https://www.nature.com/articles/s41586-018-0698-6)

**Evidence boundary:** Co-expression can nominate a potential interaction; it does not measure binding, signal transmission, or causality. A good hypothesis includes an alternative and a perturbation that could distinguish them.

## 9. Identify where a bulk signal normally originates

**What it asks:** Are candidate genes broadly expressed or enriched in particular placental cell populations?

**Why it matters:** Cell identity can help translate a tissue-level result into a biologically plausible source.

**Biology to learn:** Cell-type enrichment versus specificity; the distinction between a marker and a functional regulator. Consult the reference atlas and gene browser for anatomical and cellular context. [Vento-Tormo et al.](https://www.nature.com/articles/s41586-018-0698-6), [organizer data guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/2-the-data.md)

**Evidence boundary:** Localization in a reference is supporting context. It neither measures cell abundance in the bulk sample nor proves which population caused the observed difference.

## 10. Use normal development to interpret disease

**What it asks:** What does a healthy developmental reference suggest about the identities, programs, or interactions implicated in affected tissue?

**Why it matters:** Understanding a normal developmental role can make a later disease-associated expression pattern interpretable.

**Biology to learn:** Implantation, normal immune regulation, trophoblast differentiation, and gestational changes. Begin with the first-trimester atlas; newer spatiotemporal work extends the developmental map. [Vento-Tormo et al.](https://www.nature.com/articles/s41586-018-0698-6), [Wang et al. 2026](https://www.nature.com/articles/s41586-026-10316-x)

**Evidence boundary:** An early healthy atlas is not a matched control for a diseased delivered placenta. A developmental analogy can motivate a hypothesis, but it cannot establish that an implantation-stage defect occurred in those pregnancies.

## How to choose after learning the biology

A useful first project can combine a focused core-data question (1, 5, or 6) with cell-source interpretation (9). Questions 3, 7, and 8 generally need more cell-resolved work for stronger claims. This is our proposed prioritization, not a required track choice.

For each candidate, complete these sentences before choosing an analysis:

1. The biological uncertainty is ...
2. The cells, compartment, and pregnancy stage are ...
3. If our explanation is right, we expect ...
4. Another explanation for that observation is ...
5. The supplied data can distinguish these explanations by ...
6. The supplied data cannot establish ...
7. A follow-up experiment that could disprove our explanation is ...

Choose the question after inspecting the released sample dictionary and group sizes. Background reading should guide the hypothesis while leaving room for the data to contradict it.

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The first pass suggested Q1, Q5, or Q6 plus Q9, then deferred the choice until the released metadata could be checked. The morning review supports a narrower first pass for the new four-graph article and short video. Lead with **Q1** (EOPE without FGR vs. iPTB, **8 vs. 5**), use **Q6** for the delivery-age-adjusted and **5-vs.-5** shared-age checks, and use **Q9** for three documented gene-browser lookups. [Q1](../../Q01-placental-gene-expression.md) · [Q6](../../Q06-delivery-timing.md) · [Q9](../../Q09-baseline-gene-localization.md) · [Article guide](../../openbio-datathon-2026/docs/article-guide.md)

**Q5** is the focused alternative if the team is most interested in fetal growth. Its within-EOPE comparison retains **4 vs. 7** in the shared-age subset. For **Q3, Q4, Q7, and Q8**, let the primary result justify a targeted follow-up. The early arm has three controls and some classes have even fewer. **Q10** can use both tissues for all ten donors in its primary comparison. The healthy atlas can help interpret the delivery-stage result by placing genes in their normal cell and tissue context. [Released data dictionary](../../openbio-datathon-2026/track-01/data/track01_data_dictionary.md) · [Q10 paired plan](../../Q10-healthy-early-pregnancy-atlas.md)

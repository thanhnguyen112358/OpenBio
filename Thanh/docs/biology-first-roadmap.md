# Learning preeclampsia through cell biology and biochemistry

Prepared November 10, 2026. Starting point: comfortable with cell biology/biochemistry, new to placental biology. This roadmap follows the user's requested scientific direction; the session sequence and study exercises are proposals.

## The objective

Arrive able to explain the clinical problem, identify the cells and signals relevant to a biological question, recognize competing explanations, and specify which observations would support or weaken a hypothesis. Computational methods should follow those decisions.

Build four linked explanations:

1. What happens to the mother and fetus clinically?
2. What normally happens in the placenta and maternal tissues?
3. Which processes may become abnormal, and how could those abnormalities produce disease?
4. Which parts of that explanation can the datathon data actually examine?

## 1. Clinical grounding: what is the syndrome?

Preeclampsia usually presents after 20 weeks of pregnancy, but can first become apparent after birth. It involves hypertension with additional abnormalities; proteinuria is common but is not required when other qualifying organ dysfunction is present. The kidneys, liver, brain, and clotting system can be affected. Placental dysfunction can impair fetal growth, and severe disease may necessitate premature delivery. The exact cause remains incompletely understood. [RCOG clinical introduction](https://www.rcog.org.uk/for-the-public/browse-our-patient-information/pre-eclampsia/)

Learn to distinguish gestational hypertension, preeclampsia, eclampsia, fetal growth restriction (FGR), and preterm birth. Distinguish when disease becomes clinically visible from when its initiating abnormalities may have arisen. Include postpartum disease: delivery does not mean every clinical abnormality immediately resolves. The management problem includes balancing maternal and fetal deterioration against the consequences of prematurity. [RCOG](https://www.rcog.org.uk/for-the-public/browse-our-patient-information/pre-eclampsia/)

**Learning output:** Explain in two minutes why this is a multisystem pregnancy disorder, why the placenta matters, and why the timing of delivery is clinically consequential. Learn clinical management conceptually; this project does not require memorizing treatment doses.

## 2. Normal biology before disease mechanisms

Start with two different locations: placental villi, which support exchange, and the implantation site, where fetal trophoblasts interact with maternal decidua and vessels. Do not treat every placental cell as fetal, or every immune cell as maternal.

Build a cell map containing:

| Cell or tissue | What to learn first |
|---|---|
| Villous cytotrophoblast | Progenitor state and differentiation into other trophoblast states |
| Syncytiotrophoblast | Fused trophoblast layer; exchange and secretory functions |
| Extravillous trophoblast | Invasion into maternal tissue and participation in vascular remodeling |
| Decidual stromal cells | Pregnancy-adapted maternal uterine tissue and local regulation |
| Decidual NK cells and maternal macrophages | Immune regulation and developmental interactions at implantation |
| Hofbauer cells | Fetal macrophages within placental villi |
| Endothelial cells | Distinguish fetal placental vessels from maternal vessels |

Use the [Vento-Tormo et al. 2018 atlas](https://www.nature.com/articles/s41586-018-0698-6) to learn cellular identities and proposed interactions in early pregnancy. The [Wang et al. 2026 atlas](https://www.nature.com/articles/s41586-026-10316-x) adds a reference across gestation, spatial organization, and trophoblast–stromal regulation. Healthy reference maps describe normal development; they are not themselves disease case-control experiments.

**Learning output:** Draw the two compartments, label maternal versus fetal cells, and explain which cells can plausibly interact in each location.

## 3. Mechanisms: separate origin, messenger, and injury

Use the following as a provisional organizing model, not a universal sequence proven in every patient:

**Abnormal placental development and/or maternal susceptibility → placental stress → altered circulating and local signals → maternal vascular dysfunction and clinical manifestations.**

There is no single established initiating cause that explains all preeclampsia. Keep developmental, immune, vascular, and metabolic explanations distinct while studying how they could interact. Differences among patients and between onset groups are part of the problem to explain, not noise to discard. [RCOG](https://www.rcog.org.uk/for-the-public/browse-our-patient-information/pre-eclampsia/), [Ackerman et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/), [Admati et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37572658/)

Begin with one experimentally supported biochemical example. sFlt1 is a soluble form of the VEGF receptor FLT1 that can bind VEGF and PlGF, reducing their availability. Maynard and colleagues combined human observations with endothelial experiments and a rat model to support its contribution to vascular dysfunction. This illustrates how a placental signal can affect distant maternal organs; it does not establish sFlt1 as the sole initiating cause. [Maynard et al. 2003](https://www.jci.org/articles/view/17189)

**Learning output:** For one mechanism, write the source cell, molecular product, target, expected cellular response, clinical consequence, and strongest alternative explanation.

## 4. What maternal–fetal communication means

Study at least three kinds of communication:

- **Contact-dependent recognition:** a cell-surface molecule on one cell interacting with a receptor on another. Maternal KIR–fetal HLA-C genetics is a useful example connecting immune recognition to placentation risk. The genetic association is not something ordinary bulk expression counts can reproduce. [Hiby et al. 2004](https://pubmed.ncbi.nlm.nih.gov/15477349/)
- **Local signaling:** secreted ligands and receptors within the decidual or placental neighborhood. The early interface atlas offers candidate interactions and a framework for identifying likely sending and receiving cells. [Vento-Tormo et al. 2018](https://www.nature.com/articles/s41586-018-0698-6)
- **Systemic signaling:** placental products acting through the maternal circulation, exemplified by the sFlt1–VEGF/PlGF axis. Local trophoblast–immune interactions and distant maternal vascular effects are related but separate questions. [Maynard et al. 2003](https://www.jci.org/articles/view/17189)

For every proposed interaction, ask: Are sender and receiver in the relevant tissue and developmental stage? Is there evidence for the protein, secretion or surface localization, receptor availability, and a downstream response? What would a perturbation experiment add?

**Learning output:** Make three communication cards using the same fields. Label each statement as observed, experimentally supported, inferred, or unresolved.

## 5. Where the research landscape stands

This is a selected orientation to the literature, not an exhaustive systematic review.

| Study | Why read it | Evidence boundary |
|---|---|---|
| [Maynard et al., JCI, 2003](https://www.jci.org/articles/view/17189) | A biochemical bridge from placental factors to endothelial dysfunction | Experimental support for an important mechanism, not a complete explanation of all disease |
| [Vento-Tormo et al., Nature, 2018](https://www.nature.com/articles/s41586-018-0698-6) | Normal early interface cell identities and inferred signaling relationships | Early healthy pregnancy; predicted communication requires further validation |
| [Ackerman et al., Hypertension, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/) | The source study behind the core dataset; EOPE, FGR, and metabolic heterogeneity | Tissue collected at delivery; expression patterns do not establish initiating events |
| [Admati et al., Med, 2023](https://pubmed.ncbi.nlm.nih.gov/37572658/) | Cell-resolved differences in early- and late-onset disease; directly relevant to the optional extension | Onset-related differences in this cohort should not be turned into a universal two-category explanation |
| [Wang et al., Nature, April 2026](https://www.nature.com/articles/s41586-026-10316-x) | Normal development across gestation, spatial niches, and multimodal cell-state maps | A healthy reference is not a gestationally matched disease-control cohort |
| [Sanchez-Corrales et al., Science Advances, July 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13398505/) | Severe disease across multiple interface tissues, with spatial context and gestational-age adjustment | Tissue-specific molecular signatures and candidate interactions still require causal testing |

The working landscape is moving from cataloging altered genes toward identifying the relevant cell state, anatomical neighborhood, pregnancy stage, and maternal or fetal contribution. The 2026 studies help frame those questions; they are not additional required datathon datasets.

## 6. Understand the dataset's biological context

Read the Ackerman paper before choosing a question. It studied 30 pregnancies and separated villous tissue from decidua basalis. Its five groups comprised EOPE with FGR (7), normotensive preterm FGR (5), EOPE without FGR (8), spontaneous preterm birth (5), and uncomplicated term birth (5). The datathon guide specifies the villous-tissue bulk core; the released dictionary must determine the actual usable fields. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/), [organizer data guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/2-the-data.md)

Ask what a comparison means biologically. A term control provides a different comparison from a preterm control. FGR can occur without maternal preeclampsia. A delivered placenta is a late observation of a process that may have begun much earlier. The paper's metabolic findings occurred in a subset of cases, supporting careful attention to heterogeneity. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/)

Two methodological implications follow from the supplied gene-level bulk counts:

- A signal can reflect cell abundance, expression within cells, or both.
- An FLT1 gene count does not by itself measure the soluble sFlt1 isoform, circulating protein, or receptor signaling activity. That distinction follows from the assay's resolution and sFlt1 biology. [Data guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/2-the-data.md), [Maynard et al.](https://www.jci.org/articles/view/17189)

**Learning output:** Explain what was sampled, when, from whom, and which biological quantities remain unmeasured.

## A proposed study week

The event starts Saturday, November 14. There are four calendar preparation days from November 10 through 13, so front-load the foundations. Adjust session length to available time.

| Date | Main focus | Concrete output |
|---|---|---|
| Tue Nov 10 | Clinical syndrome, consequences, and basic placental orientation | One-page clinical explanation and vocabulary list |
| Wed Nov 11 | Normal anatomy, cell identities, implantation, and vascular remodeling | Labeled maternal–fetal cell map |
| Thu Nov 12 | Mechanisms and communication; Maynard plus selected atlas figures | Three communication cards, each with an evidence label |
| Fri Nov 13 | Ackerman dataset context and all ten question explanations | Two candidate questions with biological rationale and competing explanations |
| Sat Nov 14 | Kickoff: reconcile the released material with the biological model | One justified question and comparison selected using actual metadata |
| Sun Nov 15 | Interpret results, challenge explanations, and complete the submission | A supported hypothesis, alternatives, and a falsification experiment |
| Mon Nov 16 | Deeper reading of the 2026 studies and reflection | Identify what the datathon could not resolve and what study would come next |

The organizer's dated event website supports the November 14–15 format; the November 16 activity is optional follow-up. [Event overview](https://www.openbiolab.org/openbio-datathon)

## How to read each paper

Record six things: clinical question; tissue and pregnancy stage; comparison groups; proposed mechanism; evidence supporting it; and the claim the experiment cannot establish. Separate the paper's conclusion from your interpretation. Read methods and figures before borrowing an analysis or asserting a mechanism.

The source screening for this roadmap used organizer guides, the RCOG clinical page, the Maynard article, and abstracts/metadata for the selected atlas and disease studies. The reading assignments call for deeper figure- and methods-level study; that study is not represented as completed here.

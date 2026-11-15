# Zeisel/Admati browser lookups: how the plots were read

The browser was consulted on 2026-11-15 for four genes: `FLT1`, `NDRG1`, `HSD17B1` and `HK2`. The recorded observations are in [`browser_lookups.csv`](browser_lookups.csv).

**The browser page has no written legend.** The reading below was inferred from the plot layout and should be stated as an inference.

- Each cell class has four violins.
- A grey bar runs under the violins. Light grey covers violins 1–2 and dark grey covers 3–4. The bar image is named `early_late_bar.png`, so light grey was read as **early onset** and dark grey as **late onset**.
- A row of alternating open and filled squares runs under the violins. Open was read as **control** and filled as **PE**. The marker shapes follow the same pattern: circle = early control, star = early PE, diamond = late control, square = late PE.
- Asterisks sit above each control–PE pair. Their colors match the relative position of the medians. **Red** marks the pair where the PE median is higher and **blue** the pair where it is lower. For example, `HSD17B1` has blue asterisks where the early PE median is below the early control median, and a red asterisk in the late SCT pair where the late PE median is higher.

**What these lookups can and cannot support.** The browser is a different cohort, analyzed with different methods. Following the organizer guide, use it as supporting context for where genes are expressed, not as validation, and do not present screenshots as quantitative evidence.

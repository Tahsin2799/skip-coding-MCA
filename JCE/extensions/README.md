# Post hoc extensions: runbook

Built 10 October 2026 for the open items of the 7 October editor review. All of
it is post hoc: it writes nothing into `methodology/v2` or
`outputs/methodology/production_v03`, sources the production generator, arms
and metric code unchanged, and changes none of the 132 prespecified verdicts.

| Study | Configurations | Question (review item) |
|---|---|---|
| `sweep` | advice loading on *V* of 0, .25, .5, .75, 1 x contact and diagnosis filters (10) | Where does imputation start to beat dropping the items as advice carries more information about need? (B4) |
| `devsize` | development samples of 2,000 and 4,000 under the diagnosis filter (2) | Are the imputers' losses under the diagnosis filter small-sample failures of the item models? (A6) |
| `devcontact` | development contact of 20 and 80 percent x both filters (4) | Is the loss of the skip codings a development-to-target shift, or a property of the target's own contact level? (A1) |

All at the Q3 reference cell (*r²* .6, no recognition age slope). Every
configuration scores the ten production methods plus missing fuzzy average,
four principal component indices on indicator variables (the extension of the
theory to asset-type indices: skip as 0, mean imputation, own indicators,
items dropped) and regularized iterative MCA (Josse et al. 2012, `missMDA`).
Diagnosis-filter configurations also score imputation with Firth-penalized
item models and intercept correction (FLIC; Puhr et al. 2017, `logistf`).

Configurations with a development sample of 1,000 reuse the production seed
streams of the reference cell. The loading-.25, 50-percent-contact
configurations therefore reproduce production exactly for the ten production
methods (checked by `--verify`), and every other loading or development
contact level sees the same people and the same random draws.

**Packages.** `logistf` 1.26.1 and `testthat` 3.3.2 are in the main library.
`missMDA` 1.23 requires FactoMineR 2.17 or later, while production ran 2.13 and
FactoMineR 2.17 moves the MCA arms in the fifth or sixth decimal place, enough
to break exact reproduction. So the main library keeps FactoMineR 2.13, and
`~/R/missMDA-lib` holds FactoMineR 2.17 and `missMDA`; only the regularized
iterative MCA pass puts that library first on the path (`--lib`). Each output
file records the FactoMineR version it ran under.

**Status, 11 October 2026.** Both passes ran on 10 and 11 October (640 chunks,
no chunk error, 18 failed FLIC rows in 3 replicates at 20% development contact
under the diagnosis filter) and reproduce production exactly at the reference
configuration over all 1,000 replicates. Results are in section 4.4, Figure 4,
Tables S11 to S14 and appendix A7 of the paper.

## Commands

From the repository root. `caffeinate -i` keeps the Mac awake; both passes
resume from the last completed 50-replicate chunk if interrupted.

```sh
# 0. Tests and the production reproduction check (about 5 minutes)
Rscript JCE/extensions/test_extensions.R
Rscript JCE/extensions/run_extensions.R --verify

# 1. Main pass: 16 configurations, 1,000 replicates (about 28 hours on six workers)
caffeinate -i Rscript JCE/extensions/run_extensions.R --reps 1000 \
  --studies sweep,devsize,devcontact --workers 6 \
  --out outputs/methodology/extensions_v1 > outputs/methodology/extensions_v1.log 2>&1

# 2. Regularized iterative MCA pass, same configurations (about 14 hours)
caffeinate -i Rscript JCE/extensions/run_extensions.R --reps 1000 \
  --studies sweep,devsize,devcontact --arms rimca --lib ~/R/missMDA-lib --workers 6 \
  --out outputs/methodology/extensions_v1 > outputs/methodology/extensions_v1_rimca.log 2>&1

# 3. Summary (a minute)
Rscript JCE/extensions/extensions_summary.R outputs/methodology/extensions_v1
```

Run the two passes one after the other; each uses six workers. Progress: the
main pass writes 320 files and the second pass another 320
(`ls outputs/methodology/extensions_v1 | wc -l`), and each log line gives a
chunk's minutes and error count.

**Measured cost** (six workers, 50 replicates): 6.6 minutes for a
diagnosis-filter configuration, 3.6 for a contact-filter configuration, 2.6
for the regularized iterative MCA pass. Totals at 1,000 replicates: about 20
hours for the nine diagnosis-filter configurations, 8 for the seven
contact-filter configurations, 14 for the second pass.

## Files

| File | Role |
|---|---|
| `extensions.R` | grid, generator wrapper, FLIC item models, regularized iterative MCA, principal component indices and their closed-form diagnostics, one replicate |
| `run_extensions.R` | parallel runner, chunked and resumable; `--verify` |
| `extensions_summary.R` | paired contrasts with Monte Carlo standard errors; writes `extensions_summary.rds` |
| `test_extensions.R` | generator identity with production, common draws across configurations, FLIC fits, principal component closed forms, every arm scores |

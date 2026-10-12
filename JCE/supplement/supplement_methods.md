## Technical appendix

### A1. Donor records and analysis pools

The DHS records supply covariates only. Disease, contact, diagnosis and advice
are simulated; no survey advice response is used as a simulation outcome.
Eligible records are adults with observed water, sanitation, cooking fuel, BMI,
education and work categories. Country comparisons use ages 18 to 49. The
Bangladesh adult pool extends to age 95. The Nepal complete-input pool contains
women only. Table S8 reports the eligible and complete-input counts and the
composition of every analysis pool.

The four structural pools use fixed region-by-urban/rural membership from
`quadrant_cells.csv`. They are predominantly Indian but include Albanian and
Bangladeshi strata. Q1 has low deprivation and low coherence; Q2 has low
deprivation and high coherence; Q3 has high deprivation and low coherence; Q4
has high deprivation and high coherence. Deprivation is the mean of normalized
water, sanitation, fuel and education scores. Coherence is the mean of the 15
bias-corrected pairwise Cramér's V statistics among the six inputs. Membership
is read from the frozen file, rather than recalculated in each replicate.
Earlier descriptions calling these pools exclusively Indian were incorrect;
Table S8 describes the membership used in production, whose hash matches the
run manifest.

The four structural pools, the Bangladesh adult pool and four country pools
overlap through their source surveys. They test performance across covariate
structures, not nine independent national replications. Development and target
samples are separate draws from the same pool. A source cluster can occur in
both draws or more than once within a draw; each draw of a cluster receives a
new synthetic cluster identifier and new simulated variables. This is transport
across simulated access levels within a pool, not training in one country and
testing in another.

Clusters are sampled with equal probability and replacement. The generator adds sampled clusters whole, in order, and stops before the one that would exceed the nominal size.
If the first cluster alone exceeds that size, it is truncated. Consequently,
the nominal development and target sizes are 1,000 and 20,000, but actual sizes
usually fall slightly below them. The actual target size determines every
denominator and budget. Records are not survey-weighted, and the results concern
the finite empirical pools rather than national populations.

### A2. Harmonization and the deprivation count

Water, sanitation and fuel come from PR variables `hv201`, `hv205`/`hv225` and
`hv226`; education comes from `hv109`. BMI uses `sbbm` in Bangladesh and
`ha40`/`hb40` elsewhere, divided by 100. Values at or above 6000 before division
are treated as missing. BMI categories are <18.5, 18.5 to <25, 25 to <30, and
at least 30 kg/m². Work comes from Bangladesh PR `sb308`, or IR `v717` and MR
`mv717` joined to PR records by cluster, household and person identifiers.
Age and sex use `hv105` and `hv104`, the age in years and the sex (male or
female) of each household member as recorded in the household roster.

Most categories are harmonized from value labels because numeric codes differ
across surveys. Education uses the standardized `hv109` codes 0 to 5. The
Bangladesh work mapping uses the explicit country-specific map below. All
unmatched values remain missing and are excluded from the complete-input pool.
The release includes the ordered matching rules and membership file.

| Input | Harmonized categories in generator order | Count adds one for |
|---|---|---|
| Water | on premises; other improved; unimproved; surface | unimproved or surface |
| Sanitation | improved private; improved shared; unimproved; open defecation | shared, unimproved or open |
| Fuel | clean; kerosene/coal; charcoal; biomass | any category except clean |
| Education | higher; secondary; incomplete secondary; primary; incomplete primary; none | incomplete primary or none |
| Work | not working; services; manual; agriculture | manual or agriculture |
| BMI | underweight; normal; overweight; obese | underweight or obese |

Water labels for rivers, dams, lakes, ponds, streams, canals and irrigation map
to surface; unprotected, open or dug wells and “other” map to unimproved; piped
water into a dwelling, yard or plot maps to on premises; remaining listed
improved supplies include protected sources, tubewells, boreholes, rainwater,
bottled water, tankers and carts. Sanitation labels for no facility, bush or
field map to open; pits without slabs, hanging or bucket toilets, unspecified
discharge and other listed unimproved facilities map to unimproved. Improved
facilities are divided by reported sharing. Unknown sharing is not separately
excluded in the implemented map.

Clean fuel includes electricity, solar, LPG, natural gas, biogas and alcohol
fuels. Kerosene, coal and related liquid fossil fuels share one rung; biomass
includes wood, crop residues and dung. “No food cooked” and “other” fuel records
are excluded. Work label matching gives agriculture precedence over service or
manual labels. Bangladesh `sb308` maps 0/61/62 to not working, 11/41/51/52 to
services, 16/21/22/23/31 to manual, and 12/13/14/15 to agriculture.

The count is the sum of six binary indicators, ranging from 0 to 6 before the
shared orientation step. Water, sanitation, fuel and schooling adapt the global
MPI definitions (Alkire, Kanagaratnam and Suppa, 2018), rather than reproducing
the full household MPI. Travel time for water is not used, education is
individual, and the combined kerosene/coal rung is counted as deprived. BMI
uses the WHO category thresholds (WHO, 2000), with underweight and obesity
flagged but overweight unflagged. The work indicator is a study design choice.
The score is therefore a prespecified six-domain count, not a validated MPI or
an externally validated measure of healthcare need.

### A3. Imposed mechanisms

Let \(X_{ij}\) denote the six categories for person \(i\) in synthetic cluster
\(j\). For each non-BMI variable with \(K\) canonical categories in the order
above, the category score is \((k-1)/(K-1)\). Empty categories retain their
positions. BMI scores are .70, 0, .50 and 1 for underweight, normal, overweight
and obese. In the order water, sanitation, fuel, education, work and BMI, the
loadings are \(\lambda=(.8,.7,.6,.6,.5,.5)\).

With \(S_{ij}=\sum_k\lambda_k s_k(X_{ij})\), vulnerability is

\[
V_{ij}=\sqrt{r^2}\frac{S_{ij}-\mu_{S,p}}{\sigma_{S,p}}
       +\sqrt{1-r^2}\,\epsilon_{ij},\qquad \epsilon_{ij}\sim N(0,1).
\]

The center and scale are computed from the full donor pool for the resolved
specification and reused in development and target samples. The core varies
\(r^2\) over .4, .6 and .8. Both the construct and its weights are design choices.

Disease is drawn independently conditional on the generated predictors:

\[
Y_{ij}\sim\operatorname{Bernoulli}\left[\operatorname{expit}\left\{
\alpha_p+.8V_{ij}+u_j+\log(1.74)(a_{ij}-\bar a)/10
 +.1\,1(\mathrm{male})\right\}\right],
\quad u_j\sim N\left(0,\frac{.05}{.95}\frac{\pi^2}{3}\right).
\]

Here \(\bar a\) is the sample mean age. The age coefficient comes from Abba et
al. (2021). The .05 ICC is a residual latent-scale value, conditional on
\(V\), age and sex, set by calibration; observed disease clustering is higher
because \(V\) is itself clustered. The sex and vulnerability
coefficients are design choices. Population-specific intercepts target roughly
25% disease prevalence.

For each development or target draw, let \(Z_j\) be the standardized mean
\(V\) of cluster \(j\). Cluster access supply is
\(L_j=\rho Z_j+\sqrt{1-\rho^2}\eta_j\), where \(\eta_j\sim N(0,1)\) and
\(\rho=-.25\). A finite draw realizes this correlation only approximately.
With independent \(e_{ij}\sim N(0,1)\), contact under setting \(s\) is

\[
C_{ij,s}=1\left[U_{C,ij}<\operatorname{expit}
 (\gamma_s+.4e_{ij}+.3V_{ij}+.736L_j)\right].
\]

The baseline contact intercepts are −1.5831718277591, .0030644116036 and
1.5892429868911 for the 20%, 50% and 80% settings. The .736 access coefficient
was calibrated against a Nepal contact proxy; its value and the -.25 design
correlation are not direct estimates of these latent mechanisms. The individual
term \(.3V_{ij}\) pushes the other way, so the reference design is nearly
access-neutral at cluster level. In probe replicates (25 per cell, Q3 and
Bangladesh 18 to 95 at \(r^2=.6\), `JCE/analysis/editor_review_checks.R`), the
cluster-level Spearman correlation between contact rate and mean \(V\) was
−.04 to −.05 in Q3 and −.02 in Bangladesh 18 to 95 at every contact level.
Under the diagnosis filter the cluster eligibility rate correlated .07 to .13
with cluster need at 20% contact and .22 to .28 at 80%. The access-correlation
−.5 variation is the inverse-care scenario. Probe values are diagnostics, not
production estimates.

Let \(z(a)\) standardize age within the draw. Diagnosis is

\[
q_{ij}=\begin{cases}
\operatorname{expit}\{\operatorname{logit}(.50)+\delta z(a_{ij})\},&Y_{ij}=1,\\
.05,&Y_{ij}=0,
\end{cases}\qquad
R_{ij,s}=C_{ij,s}1(U_{R,ij}<q_{ij}),
\]

where \(\delta\) is 0 or .5. The direction of the age association is motivated
by Foti and Hardy (2025); its magnitude is a design choice. The .50 and .05
values are design anchors, not estimates of diagnosis conditional on contact
from population awareness or self-report validation studies.

For advice item \(k\), with baseline probabilities
\(b_k=(.10,.10,.05,.15)\), factual receipt and observation are

\[
A^*_{ijk,s}=C_{ij,s}1\left[U_{A,ijk}<\operatorname{expit}
 \{\operatorname{logit}(b_k)+\log(16)R_{ij,s}+.25V_{ij}\}\right],
\]
\[
G_{ij,s}=\begin{cases}C_{ij,s},&\text{contact filter},\\
R_{ij,s},&\text{diagnosis filter},\end{cases}\qquad
A_{ijk,s}=\begin{cases}A^*_{ijk,s},&G_{ij,s}=1,\\
\text{not asked},&G_{ij,s}=0.\end{cases}
\]

Because \(A^*_{ijk,s}=0\) whenever \(C_{ij,s}=0\), a skip under the contact
filter is a structural zero, and coding it as "no" reproduces factual advice
receipt. Under the diagnosis filter, contacted but undiagnosed respondents can
have \(A^*_{ijk,s}=1\), so the same coding misreads them.

The odds ratio of 16 relating diagnosis to advice and the advice loading .25 are
design choices. The motivating survey cannot estimate that contrast because it
does not ask undiagnosed people about advice. There is no ordinary item nonresponse
among eligible respondents. Contact, diagnosis and advice collapse different
survey time windows into one simulated period. No intervention effect on disease
or vulnerability is generated.

Figure S2 shows the imposed structure. It omits independent noise terms and
empirical correlations among donor covariates. The arrow from vulnerability to
cluster access supply represents the imposed inverse-care association, and the
two dashed arrows into eligibility are the alternative filters.

**Figure S2. The imposed data-generating structure.** The six covariates, age
and sex are resampled from the surveys; every other node is simulated. Dashed
arrows mark the two alternative filters that decide who is asked about advice.
File: `JCE/figures/figS2.pdf`.

The outcome intercepts for Q1, Q2, Q3, Q4, Bangladesh adults, India, Nepal,
Albania and Bangladesh ages 18 to 49 are, respectively, −1.4018739424817,
−1.3874000585407, −1.3979337882166, −1.3897578168946, −1.5156911575790,
−1.3925435620646, −1.3550533785203, −1.3936592700210 and −1.3962307474406.

The stress settings each replace one parameter. They are diagnosis given disease
.20 or .80; false-positive diagnosis .15; access correlation 0 or −.5; access
coefficient 0; two or eight advice items; all six loadings .65; nominal
development size 500; and 20 imputations. Two-item advice uses the first two baseline probabilities;
eight-item advice repeats the four-item vector twice. The low and high
diagnosis settings are also run under the diagnosis filter. Other stress
settings and all mechanism checks use the contact filter at the Q3 reference
configuration (\(r^2=.6,\delta=0\)).

The mechanism checks set both diagnosis probabilities to .1625; replace advice
by four always-observed proxies with probability
\(\operatorname{expit}(1.2V)\); increase the advice loading from .25 to 1;
or set BMI scores to 0, 1/3, 2/3 and 1. When contact calibration changes, the
20%/50%/80% intercepts are:

| Variant | 20% | 50% | 80% |
|---|---:|---:|---:|
| Access correlation 0 | −1.5940221939971 | .0042973493238 | 1.6014012231968 |
| Access correlation −.5 | −1.5978230961461 | −.0233018728466 | 1.5523706500620 |
| Access coefficient 0 | −1.4580612563225 | .0002043784744 | 1.4584134522060 |
| Equal loadings | −1.6008579919599 | −.0143568055093 | 1.5725634087540 |
| Monotone BMI | −1.6009325417957 | −.0146981182331 | 1.5705678407246 |

### A4. Fitting, projection and allocation

All methods receive the same generated data within a replicate. The development
sample is built at 50% contact and reused for all target contact settings. The
two filters share the same populations and uniforms. Per-purpose seed streams
separate generation, development fitting and target scoring. Uniforms for
individual and cluster budget ties are drawn before any method is fitted and
do not depend on truth, disease or diagnosis.

MCA uses one dimension. Each development score is split at its unique empirical
quintile boundaries. If disease prevalence in the top tier is below that in the
bottom tier, its sign is reversed and the cuts are recomputed. This orientation
also applies to the count, oracle and supervised reference. Scores with ties
can produce fewer than five tiers. Target scores use the development orientation
and cut-points; values outside the development range are clamped to its bounds
for tier assignment. Budget rankings use the projected score before clamping.

Each imputer fits one logistic model per item among eligible respondents using
the six inputs, age and sex. Constant terms are omitted. For item prediction,
an unseen covariate category is replaced by the modal category among eligible
fitting respondents; the count of replacements is retained. Five sets of
Bernoulli responses are drawn from fixed fitted probabilities, keeping observed
responses. Parameter uncertainty and residual dependence between items are not
propagated. Each completed development dataset yields an MCA axis, oriented and standardized on its development coordinates, and each person's five standardized scores are averaged into one before allocation. The target-model variant
refits item models in the target; both variants retain the development MCA axes.

Under the diagnosis filter the development item models rest on few people. In the probe replicates of A3 each was fitted on 83 to 86 diagnosed respondents. With 19 to 21 parameters and 19 to 21 people in the smaller outcome class, that is about one event per parameter, and R reported fitted probabilities of 0 or 1. Under the contact filter the models had about 490 respondents and about three events per parameter.
In production, against dropping the items under the diagnosis filter, the mean
differences in cluster rank correlation were −.087 for development models and
−.055 for target models, the medians −.028 and −.009, and 23.9 and 14.4 percent
of replicates trailed by more than .10. The mean deficit comes from a tail of
replicates with large deficits; the typical target-model replicate was within .01 of dropping.

Disease measures are the development and target top-tier odds ratios, the
development within-cluster and between-cluster slopes, and the target AUCs of
the index and of the development-fitted tier, age and sex model. Each is
evaluated on the people its method scores. Regression results flagged for
separation are excluded as the metric specification states (A5, Table S6).

The supervised reference is logistic regression on all categorical inputs and
zero-coded advice, without age and sex in its score. Its development scores use
five folds separated by synthetic cluster; a model fitted on all development
rows scores the target. Age and sex enter the downstream tier validation models
for every method. MCA and supervised projection leave target rows unscored if
they contain categories absent from the relevant development fit. Subset MCA excludes the not-asked categories but keeps the full-data masses of the others. Its target rows are placed by the transition formula, because FactoMineR's `predict.MCA` ignores the exclusion.

Individual budgets select the highest-scoring people up to
\(\lfloor bN\rfloor\), for \(b=.05,.10,.20\). Unscoreable people are not
selected. Cluster budgets select up to \(\lfloor .10J\rfloor\) clusters using
the mean score of scored residents. Clusters with no scored resident are never selected; for
whole-target rank correlation they tie below every scored cluster. Slots remain
unfilled when too few people or clusters can be scored. Equal numbers of
clusters need not mean equal numbers of residents or equal delivery costs.
Table S5 reports the share of budget slots used and the residents covered.

The high-need set is the top fifth of *V* in the whole target, fixed before any
method is evaluated. Undiagnosed disease means \(Y=1,R=0\). Cluster truth uses
all residents, including unscored residents. Inclusion error is the fraction
of selected clusters below median truth, whereas top-decile precision is the
fraction belonging to the true top decile. The high-need share selected is the selected
high-need count divided by the whole-target high-need count. The corresponding
shortfall per 1,000 is
\(1000[\min(\lfloor .10N\rfloor,n_H)-n_{H,selected}]/N\).

### A5. Reporting correction and sensitivity analysis

The original calculation for Supplementary Figure S1 (Figure 3 of the
27 September draft) omitted the separation exclusions stated in
the metric specification. The corrected calculation masks only affected disease
metrics, preserving the targeting outcomes on the same rows. Development
separation masks the development top-tier odds ratio and within/between slopes;
target separation masks the target top-tier odds ratio and the two target AUCs,
as the frozen metadata specify. The operational flag is any fitted tier-model
probability below \(10^{-8}\) or above \(1-10^{-8}\); it is not a formal proof
of separation. Table S6 gives exclusions, and `correlation_correction.csv`
retains the original and corrected medians.

Within each core scenario, correlations pool the ten methods and 1,000
replicates, using pairwise available values after exclusions. Inclusion error
and shortfall are sign-reversed. The resulting 54 correlations include dependent
and algebraically redundant endpoints; they are descriptive comparisons, not
independent hypothesis tests. Supplementary Figure S1 summarizes each pairing by its median
and range across 72 scenarios. Its bars are not confidence intervals.

A post-review sensitivity first averages valid results for each method within
each scenario, then correlates the ten method means. Table S7 reports the six
disease measures against the three principal targeting metrics. This check was
added after production and was not part of the frozen decision rule.

The correction is implemented in `JCE/analysis/review_analysis.R`. It checks the
production source hashes, donor cache and quadrant membership, then reads the
saved core results. It does not alter the simulation, frozen summarizer,
production records, or the 132 paired targeting verdicts. Protocol development
followed the 17 to 19 September pilots; the amended design was fixed on
20 September, the separate sizing pilot ran on 21 September, and production ran
on 23 to 25 September 2026. Proposed margins were confirmed unchanged on
25 September. Neither the amendment nor this reporting correction is a
preregistration.

The separate 100-replicate sizing pilot projected 124 of the 132 targeting
contrasts to be decided at 1,000 replicates, with no further gain at
2,000. This projection used the pilot's paired spreads and held its mean
differences fixed. Pilot results are not pooled with production, which used a
different random stream and manifest. In production 122 were
decided, including one equivalence.

### A6. Proofs and verification of the theory

The theory, its tests and missing fuzzy average were derived after
production and are post hoc. They change no prespecified verdict.

**Notation.** There are *n* respondents and *Q* = *P* + *k* categorical items,
the substantive items \(X_1,\dots,X_P\) and the filtered items
\(A_1,\dots,A_k\). The filtered items are all asked of the same eligible
respondents, \(G_i=1\), with \(p=\bar G\). A coding turns item
*q* into indicator columns \(z_j\), \(j\in J_q\). Eligible respondents receive
the crisp indicator of their answer, and a skipped respondent receives a profile
\(v_i\) over \(J_q\), crisp or fuzzy. Subset analysis keeps a set \(S\) of
analysed columns with their full-data masses. Write \(n_j=\sum_i z_{ij}\),
\(f_j=n_j/n\), \(f_{j\mid G}\) for the answer distribution among the eligible
(zero for a column the eligible never use, such as a skip category), and
\(\bar v_j\) for the mean profile given to the skipped, so that
\(f_j=pf_{j\mid G}+(1-p)\bar v_j\). For a centered score \(h\) with
\(n^{-1}\sum_ih_i^2=1\), the squared correlation ratio of item *q* is

\[
\eta_q^2(h)=\frac1n\sum_{j\in J_q\cap S}\frac{\left(\sum_iz_{ij}h_i\right)^2}{n_j},
\]

which for a crisp item is \(\sum_jf_jm_j(h)^2\), with \(m_j\) the mean of \(h\)
in category *j*. The eligibility score is \(s_G=(G-p)/\sqrt{p(1-p)}\).

**Lemma A6.0 (the MCA criterion).** With
\(B=Q^{-1}\sum_q\sum_{j\in J_q\cap S}(z_j-f_j\mathbf 1)(z_j-f_j\mathbf 1)^\top/n_j\),
every centered unit-variance score satisfies
\(h^\top Bh/n=\lambda(h)=Q^{-1}\sum_q\eta_q^2(h)\). The first axis maximizes
\(\lambda(h)\), and \(\lambda_1=\max_h\lambda(h)\). *Proof.* Column *j* of the
standardized residual matrix is \((z_j-f_j\mathbf 1)/\sqrt{Qn_j}\), and
\((z_j-f_j\mathbf 1)^\top h=z_j^\top h\) for centered \(h\). Subset analysis keeps
only the columns in \(S\) (Greenacre and Pardo 2006). For complete crisp data
this is the homogeneity form of MCA (Le Roux and Rouanet 2010). ∎

**Theorem 1 (squared correlation ratio with eligibility).** The squared correlation ratio of item *q* with eligibility, its discrimination measure on the eligibility score in the language of homogeneity analysis (Michailidis and de Leeuw 1998), is

\[
\eta_q^2(s_G)=p(1-p)\sum_{j\in J_q\cap S}\frac{(\bar v_j-f_{j\mid G})^2}{f_j}.
\]

*Proof.* \(\sum_iz_{ij}s_{G,i}=\big[np(1-p)f_{j\mid G}-n(1-p)p\,\bar v_j\big]/\sqrt{p(1-p)}
=n\sqrt{p(1-p)}\,(f_{j\mid G}-\bar v_j)\). Square, divide by \(n\,n_j=n^2f_j\) and
sum over the analysed columns. ∎ When every column is analysed, \(\eta_q^2(s_G)\) is
the mean-square contingency between eligibility and the coded item.

**Corollary A6.1 (the codings).** For a binary item that the eligible answer
"yes" with probability \(\pi\):

(a) Missing single, the skip as its own category, with \(\bar v=(0,0,1)\) on (no, yes,
skip), gives \(\eta_q^2(s_G)=p(1-p)\,[(1-\pi)/p+\pi/p+1/(1-p)]=1\), whatever the number of answer
categories.

(b) Subset MCA, with the skip column excluded, gives
\(\eta_q^2(s_G)=p(1-p)\,[(1-\pi)/p+\pi/p]=1-p\). Subset MCA gives the same solution as
missing passive modified margin (Josse et al. 2012), so modified margin keeps
\(1-p\) as well. The "multiple" variant of subset MCA for incomplete data, with
one missing category per cell, has the same analysed columns and the same
value.

(c) The skip coded "no", with \(\bar v=(1,0)\), has \(f_{\text{no}}=1-p\pi\) and
\(f_{\text{yes}}=p\pi\), so
\(\eta_q^2(s_G)=p(1-p)\,[\pi^2/(1-p\pi)+\pi/p]=\pi(1-p)/(1-p\pi)\). It is below the
subset value exactly when \(\pi<1/(1+p)\).

(d) Missing fuzzy average, with \(\bar v=f_{\cdot\mid G}\), gives \(\eta_q^2(s_G)=0\). Since every
term of Theorem 1 is non-negative, \(\eta_q^2(s_G)=0\) if and only if
\(\bar v_j=f_{j\mid G}\) for every analysed column with \(f_j>0\). Among codings
that analyse every column, it is the only mean profile that carries no filter.

(e) Under imputation, \(\bar v\) is the mean imputed profile, and \(\eta_q^2(s_G)\) is near 0
when that profile is near the answers of those asked.

**Corollary A6.2 (a bound under missing fuzzy average).** Under missing fuzzy average,
\(\eta_q^2(h)\le p\cdot n^{-1}\sum_{i:G_i=1}h_i^2\le p\) for every centered
unit-variance \(h\). The filtered block analysed alone therefore has
\(\lambda_A\le p\), and \(\kappa\le kp/(P\lambda_X)\). *Proof.* The skipped rows of
each column equal \(f_j\), so with \(S_j\) the sum of \(h\) over the eligible in
category *j* and \(H=\sum_jS_j\), \(\sum_iz_{ij}h_i=S_j-f_jH\) because \(h\) is
centered, and \(\eta_q^2=n^{-2}\big[\sum_jS_j^2/f_j-H^2\big]\). The eligible in
category *j* number \(npf_j\), so by the Cauchy-Schwarz inequality
\(S_j^2\le npf_j\sum_{i\in E_j}h_i^2\); sum over *j*. ∎ Van der Heijden and
Escofier (2003) note that missing fuzzy average gives the lowest total inertia;
the corollary bounds the inertia it leaves to the items.

**Proposition A6.3 (a lower bound).**
\(\lambda_1\ge\lambda(s_G)=Q^{-1}\sum_{q}\eta_q^2(s_G)\ge
Q^{-1}\sum_{q\le k}\eta_q^2(s_G)\), indexing the filtered items first. Under missing
single, \(\lambda_1\ge k/Q\), and the filtered block analysed alone
has eigenvalue exactly 1, so the filter is an exact dimension. In the core design
\(k/Q=.40\); in probe replicates the first eigenvalue was .404 to .407, with an
axis correlating .995 to .998 with eligibility.

**Proposition A6.4 (block dominance, the inertia ratio and instability).** Let
\(h_X\) and \(h_A\) be the first axes of the substantive block and of the coded
filtered block analysed alone, with eigenvalues \(\lambda_X\) and \(\lambda_A\).
On the span of \(\{h_X,h_A\}\) the criterion is the generalized eigenproblem
\(Mw=\mu Cw\), with \(M_{ab}=Q^{-1}\sum_q\eta_q(h_a,h_b)\), the cross form of the
squared correlation ratio, and \(C\) the Gram matrix of the two scores. Its
leading value is a Rayleigh-Ritz lower bound for \(\lambda_1\), exact when the
first axis lies in that span. The diagonal entries are
\(Q^{-1}[P\lambda_X+\sum_{q\le k}\eta_q^2(h_X)]\) and
\(Q^{-1}[k\lambda_A+\sum_{q>k}\eta_q^2(h_A)]\). The leading term of each diagonal entry, \(P\lambda_X\) or \(k\lambda_A\), is the block's number of items times its maximum axial inertia, the first eigenvalue of its separate analysis. It is what the block contributes along its own first axis in the joint analysis. When the blocks are weakly associated the cross terms are small and the leading direction follows the larger diagonal, so the filtered block dominates as the inertia ratio \(\kappa=k\lambda_A/(P\lambda_X)\) passes 1. Multiple factor analysis divides each group by its first
eigenvalue so that no group dominates (Escofier and Pagès 1994); \(\kappa\)
measures the dominance an undivided filtered block exerts. By A6.3,
\(k\lambda_A\ge\sum_{q\le k}\eta_q^2(s_G)\), so Theorem 1 bounds \(\kappa\) from below, with
\(\kappa\ge k/(P\lambda_X)\) under missing single and
\(\kappa\ge k(1-p)/(P\lambda_X)\) under subset MCA.
For binary items, \(k\lambda_A\) is the largest eigenvalue of their correlation
matrix, at least \(1+(k-1)\bar r\). In probe replicates with
\(|\kappa-1|>.25\), the span of \(\{h_X,h_A\}\) explained a median .997 of
the first axis (5th percentile .954).

By the Davis-Kahan theorem in the form of Yu, Wang and Samworth (2015), the
angle between a sample axis and its population counterpart is bounded by twice
the sampling perturbation of \(B\) divided by the eigengap
\(\lambda_1-\lambda_2\). In the two-block reduction the gap closes as
\(\kappa\to1\), so development samples near the threshold give axes that swing
between the blocks.

**Proposition A6.5 (cluster ranking by a weighted score, and the spread of eligibility).** Write
the deployed score as \(s=\alpha h_X+\beta h_A\) plus a remainder. At cluster
level in the target, with between-cluster standard deviations
\(\sigma_X,\sigma_A\), correlations \(\rho_X,\rho_A\) of the cluster means with
true cluster need and \(r\) between them,

\[
\tau=\frac{\alpha\sigma_X\rho_X+\beta\sigma_A\rho_A}
{\sqrt{\alpha^2\sigma_X^2+\beta^2\sigma_A^2+2\alpha\beta\sigma_X\sigma_Ar}},
\]

exact for the Pearson correlation of cluster means when the remainder is zero.
As \(\beta/\alpha\to\infty\), \(\tau\to\rho_A\). When the filtered block dominates and carries most of the filter (missing
single, or subset MCA with few people asked), \(h_A\) is close to eligibility,
and as \(\beta\sigma_A\) outgrows \(\alpha\sigma_X\) the index ranks communities
as the cluster eligibility rate does. A block can dominate without carrying
much of the filter, as subset MCA does when developed at 80% contact (A7). The weights are fixed in development, but \(\sigma_A\) is
a property of the target. If cluster eligibility rates vary with logit-scale
variance \(\sigma_L^2\) around *p*, the delta method and binomial sampling within
clusters of mean size \(\bar n\) give
\(\sigma_A^2\approx[p(1-p)]^2\sigma_L^2+p(1-p)/\bar n\), largest at
\(p=\tfrac12\). Under the contact filter the target's *p* is its contact rate,
so the filter's pull peaks at 50% contact. Under the diagnosis filter *p* stays
below one half, and both \(\sigma_A\) and \(\rho_A\) rise with contact. With few
people eligible, \(\sigma_A\) is small and the substantive part still orders
many clusters, so the index can rank them better than the eligibility rate
does (A6.11). In probe
replicates, the development-time two-direction reduction predicted target
cluster targeting with r = .94 and mean absolute error .039 where
\(|\kappa-1|>.25\).

**Proposition A6.6 (a budget below the share asked).** If the index ranks every eligible
person above every skipped person and the budget fraction is \(b\le p\), every
selected person is eligible, and the selected share of any group outside the
filter is zero. If \(b>p\), at most \((b-p)/b\) of the slots go outside it.
*Proof.* Counting. ∎ The skip coded "no" meets the premise only approximately, since
an eligible respondent who answered "no" to every item has the same coded
items as a skipped one, which is why it selected up to .00002 rather than 0.

**Proposition A6.7 (structural zeros and identification).** In the generator,
factual advice is \(A^*_{ijk,s}=C_{ij,s}\cdot1[\cdot]\) (A3). Under the contact
filter a skipped respondent has \(A^*=0\) with certainty, so the skip coded "no"
reproduces factual advice exactly; any harm comes from the item itself, since at a given need a person's chance of having been advised depends on their access. Under the diagnosis filter, contacted
but undiagnosed respondents may have been advised, and observed advice comes
only from the diagnosed, so their factual value is not identified from the
observed data and every completion is an assumption. Missing fuzzy average
assigns no factual value. An item then informs the index only through
differences among those asked.

**A6.8 Principal component indices.** For binary items analysed by principal
components of their correlation matrix, the first eigenvalue is
\(\max_h\sum_j\operatorname{corr}(h,x_j)^2\). An item coded 0 for the skipped,
\(x=G\cdot Y\) with \(Y\) the answer, has
\(\operatorname{cov}(G,x)=p\pi(1-p)\), \(\operatorname{var}(x)=p\pi(1-p\pi)\)
and so \(\operatorname{corr}(s_G,x)^2=\pi(1-p)/(1-p\pi)\), the value of
Corollary A6.1(c). Propositions A6.3 to A6.6 carry over with correlation-matrix
eigenvalues. A7 tests this extension by simulation.

**A6.9 Verification on the production replicates.**
`JCE/theory/theory_production_link.R` regenerates each production replicate
from its seed stream. Regeneration reproduced the saved true \(V\) and the
production dropped-items rank correlation exactly in all 52 configurations. In
each development sample it computes, for missing single, subset MCA and the
skip coded "no", the empirical squared correlation ratio with eligibility and the value of
Theorem 1, \(\lambda_X\), \(\lambda_A\), \(\kappa\), the first two eigenvalues,
and the squared correlations of the index's first axis with \(h_A\) (block
R²), with eligibility (eligibility R²) and with \(h_X\). Missing fuzzy average
is fitted as correspondence analysis of the fuzzy indicator matrix
(FactoMineR `CA`); development and target rows are placed by the transition
formula on the development column standard coordinates, which reproduced the
fitted row coordinates to \(5\times10^{-14}\). In the target, the skipped
receive the answer distribution of the target's own eligible respondents, so no
model is fitted, and rows with a category the development sample never saw are
left unscored, as for every MCA method. It is scored by the production metric
code unchanged. The full run (1,000 replicates of 52 configurations, 51,000
with filtered items) took 2 hours 51 minutes on six workers on 8 October 2026;
`JCE/theory/theory_link_summary.R` produces Table 2, Figure 1, Table S10 and the
numbers below.

- *Theorem 1.* The maximum difference between the empirical value and the
  formula over 204,000 coding-by-replicate rows was \(6\times10^{-15}\), and the
  skip coded "no" carried less filter than subset MCA in 51,000 of 51,000 replicates.
- *Inertia ratio.* Over 153 configuration-coding pairs, its Spearman
  correlation with the loss against dropping the items was .83 (.88 for the block R², .79
  replicate by replicate). Above \(\kappa=1\) the block took on average .99 of
  the axis under missing single, .74 under subset MCA and .86 under
  the skip coded "no"; the two pairs below 1 took .15 and .24 and lost −.005 and .000.
- *Instability.* The within-configuration replicate standard deviation of targeting was .12 for \(|\kappa-1|\) up to .15, .16 from .15 to .30, .08 from .30 to .60 and .06 above .60, with mean eigengaps .036, .038, .099 and .195 in the same bands. The Spearman
  correlation between \(|\kappa-1|\) and the standard deviation was −.64.
- *Budget below the share asked.* In every diagnosis-filter cell whose diagnosed share exceeded the
  10% budget, the three codings selected 0 to .00001 of the people with
  undiagnosed disease, against .17 for dropping the items; below the budget
  they selected .03 to .12.
- *Corollary A6.2.* On probe replicates (`JCE/theory/theory_checks.R`, 12 per
  cell), the fuzzy-average block eigenvalue reached at most .64 of *p* in all 12
  settings, and its inertia ratio was .40 to .63 under the contact filter (.89
  with eight items) and .03 to .09 under the diagnosis filter.

**A6.10 A published simulation.** Josse et al. (2012, Table 3, rows 7 and 8)
report that modified margin "can really crash" when 60 of 100 individuals miss
variables 1 to 3 together and 60 miss variables 9 and 10 together, without a
stated reason. Theorem 1 gives modified margin a squared correlation ratio with the missingness pattern of .6 per item
there and missing fuzzy average zero. `JCE/theory/josse2012_check.R` rebuilds
the design with 500 replications per cell and reproduces all six published
medians of the modified RV coefficient within .01. For weak and strong
structure they are .30 and .45 for missing single, .42 and .71 for modified
margin, and .77 and .91 for missing fuzzy average. The missingness indicators explained .50, .48 and .01 of the
three methods' first two axes under weak structure, and .50, .23 and .01 under
strong structure. The collapse is the missingness pattern dominating the first axes.

**A6.11 Cluster ranking under a dominant filter.** Production saved person-level
data for two replicates of every configuration, 24 per filter in the core
design; `JCE/analysis/filter_ranking_check.R` reads them. Under the contact
filter missing single ranked clusters about as the cluster contact rate did.
Their mean Spearman correlations with cluster need were .04 and −.04 at 20%
contact, .02 and −.03 at 50% and .05 and −.03 at 80%. Under the diagnosis
filter missing single ranked clusters better than the diagnosis rate, at .41,
.33 and .36 against .10, .18 and .25. Only 4%, 9% and 14% of the target were
eligible, so eligibility varied little between clusters (A6.5), and 54%, 23%
and 10% of clusters had no diagnosed resident. Among those clusters the rest
of the score still ordered need, with a correlation of about .83 at every
contact level. At 80% contact under the diagnosis filter, everyone that
missing single, subset MCA and the skip coded "no" selected at the 10%
individual budget had been diagnosed, against .22 for dropping the items.

### A7. Additional simulations (post hoc)

Three studies, designed after the theory and run on 10 and 11 October 2026,
take up three questions the production design left open. They ask how
informative the advice must be before imputing it helps, whether the imputers' losses under the
diagnosis filter come from small item models, and whether the loss of the
skip codings depends on the development sample's contact level or on the
target's. All use the reference scenario of pooled stratum Q3 (\(r^2=.6\), no
diagnosis age slope), a nominal target of 20,000 evaluated at 20%, 50% and 80% contact, and 1,000 replicates per configuration.

| Study | Configurations | What varies |
|---|---:|---|
| Advice loading | 10 | the coefficient of \(V\) in the advice model (A3), 0, .25, .5, .75 or 1, under each filter |
| Development size | 2 | nominal development samples of 2,000 and 4,000, diagnosis filter, loading .25 |
| Development contact | 4 | contact in the development sample 20% or 80%, under each filter, loading .25 |

Configurations with a development sample of 1,000 reuse the production seed
stream of the reference configuration, so the people and the uniforms are the
same in every one of them and only the advice, or the development sample's
contact, changes. At loading .25 and 50% development contact the configuration
is production itself. Over all 1,000 replicates under both filters, the ten
production methods reproduced the saved production metrics exactly (60,000
rows, 81 metrics, largest difference 0). The larger development samples draw
their own streams.

Every configuration scores the ten production methods and missing fuzzy
average (A6.9), and four further kinds of index.

- *Imputation with FLIC item models*, under the diagnosis filter only, where
  the development item models rest on about 84 diagnosed respondents. Each
  item model's slopes are estimated by Firth's penalized likelihood and its
  intercept is then re-estimated so that the mean predicted probability equals
  the observed rate (Puhr et al. 2017; `logistf` 1.26.1). Development and
  target models are used as in the production imputers.
- *Regularized iterative MCA* (Josse et al. 2012). The skipped cells of the
  indicator matrix are imputed by a regularized reconstruction of rank two
  (`missMDA` 1.23, `imputeMCA`), and the index is the first axis of the MCA of
  the completed fuzzy indicator matrix. The target is imputed on its own data
  and placed by the transition formula. `missMDA` requires FactoMineR 2.17,
  which moves the production MCA methods in the fifth or sixth decimal place,
  so this method ran in a separate pass under 2.17 and every other method under
  2.13. Its filter diagnostics were not computed.
- *Principal component indices* on 0/1 indicators of every input category, the
  score being the first principal component of their correlation matrix, as
  for asset indices. A filtered item enters as its yes indicator with the skip
  coded 0, as the same indicator with the skip replaced by the share of those
  asked answering yes (mean imputation, the counterpart of missing fuzzy
  average), as yes, no and skipped indicators, or not at all. A6.8 gives the
  closed forms checked in every development sample.

The two passes took 17.3 and 7.9 hours of computation on six workers in 640
chunks of 50 replicates, with no chunk error. Of 822,000 method rows, 18
failed. They were both FLIC methods in 3 of 1,000 replicates with 20% development
contact under the diagnosis filter, where development held 27 to 31 diagnosed
respondents for 17 to 20 coefficients per item model and the penalized fit
stopped on an information matrix that was not positive definite. No table uses
those rows. `JCE/extensions/extensions_summary.R` produces Tables S11 to S14
and Figure 4.

- *Advice loading* (Table S11, Figure 4A). Under the contact filter the
  imputers overtook dropping the items at loadings of .16 (target models) and
  .24 (development models), against .25 in production. Under the diagnosis
  filter target models crossed at .67 and development models stayed behind at
  every loading (−.019 at 1); with FLIC item models the crossings were .47 and
  .81. Regularized iterative MCA led dropping by .003 to .041 under the contact
  filter and stayed within .002 of it under the diagnosis filter; missing fuzzy
  average led by .002 to .032 and stayed within .001. As the loading grew,
  subset MCA rose from .342 to .833 and the skip coded "no" from .525 to .791
  under the contact filter; under the diagnosis filter "no" rose from .396 to
  .517 and subset MCA stayed at .334 to .336.
- *Development size* (Table S12, Figure 4B). Under the diagnosis filter the
  imputers' mean deficits against dropping were .107 and .161 with 1,000
  development respondents (84 diagnosed, 1.4 events per item-model
  parameter) and .012 and .029 with 2,000 (2.7); with 4,000 (5.4) target
  models led by .004 and development models trailed by .001. With 1,000, FLIC
  item models cut the deficits to .033 and .061 and the share of replicates
  trailing dropping by more than .10 from .293 and .470 to .101 and .181.
  Dropping the items itself rose with development size, from .859 to .882 and
  .894.
- *Development contact* (Table S13, Figure 4C). Missing single, subset MCA and
  the skip coded "no" ranked worst at 50% target contact in 17 of the 18
  combinations of filter, coding and development contact; the exception,
  subset MCA developed at 80% contact under the contact filter, fell from .848
  to .801 as target contact rose. The development sample's contact level, more than the target's, set how well these codings ranked clusters. Under the
  contact filter, as development contact went from 20% to 50% and 80%, subset
  MCA's inertia ratio fell from 1.93 to 1.20 and 1.09 and its filtered block
  still took .992, .706 and .742 of the first axis, but eligibility explained
  .992, .625 and .002 of it, and its cluster rank correlation, over the three target contact levels, was .007 to .042, .434 to .520 and .801 to .848. With 80% asked, the block keeps
  \(1-p=.2\) of the filter (Theorem 1) and its axis follows the answers of
  those asked. Missing single kept \(\kappa=2.42\), eligibility explained .996
  of its axis at every level, and it stayed below .05. Under the diagnosis
  filter the minimum at 50% target contact also appeared; there the spread of
  cluster eligibility and its correlation with need both rise with contact
  (A6.5), and this design does not separate the two.
- *Principal component indices* (Table S14). The closed forms of A6.8 held in
  all 16,000 development samples, to \(3\times10^{-13}\). With yes, no and skipped indicators, the counterpart of missing single, the inertia ratio was 2.73 under the contact filter and 3.00 under the
  diagnosis filter, eligibility explained .991 and .997 of the first
  component, and the index's cluster rank correlation was .019 to .050 and .302
  to .304.
  With the skip coded 0 the ratio was .81 under the contact filter, where
  eligibility explained .058 and the index reached .865 to .894, against .862
  for dropping the items, and 1.13 under the diagnosis filter, where
  eligibility explained .678 and the index reached .744 to .757. Mean
  imputation has zero squared correlation with eligibility and reached .866
  to .905 and .861 to .873.

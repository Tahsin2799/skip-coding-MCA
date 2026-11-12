## Technical appendix

### A1. Donor records and analysis pools

The DHS records supply covariates only. Disease, contact, diagnosis and advice
are simulated, and no survey answer about advice is used as a simulation
outcome. Eligible records are adults with observed water, sanitation, cooking
fuel, BMI, education and work categories. Country comparisons use ages 18 to 49.
The Bangladesh adult pool extends to age 95. The Nepal pool of complete records
contains women only. Table S8 reports the counts of eligible records and of
records complete on all six inputs, and the composition of every analysis pool.

The four structural pools are fixed sets of strata, each a region and its urban
or rural part, and their membership does not change between replicates. They
are mainly Indian but include Albanian and Bangladeshi strata (Table S8c). Q1
has low deprivation and low coherence, Q2 low deprivation and high coherence,
Q3 high deprivation and low coherence, and Q4 high deprivation and high
coherence. Deprivation is the mean of normalized water, sanitation, fuel and
education scores. Coherence is the mean of the 15 pairwise Cramér's V
statistics among the six inputs, each corrected for bias.

The four structural pools, the Bangladesh adult pool and four country pools
overlap through their source surveys. They test performance across covariate
structures, not nine independent national replications. Development and target
samples are separate draws from the same pool. A source cluster can occur in
both draws or more than once within a draw, and each draw of a cluster receives
a new synthetic cluster identifier and new simulated variables. This is
transport across simulated access levels within a pool, not training in one
country and testing in another.

Clusters are sampled with equal probability and replacement. The generator adds sampled clusters whole, in order, and stops before the one that would exceed the nominal size.
If the first cluster alone exceeds that size, it is truncated. Consequently,
the nominal development and target sizes are 1,000 and 20,000, but actual sizes
usually fall slightly below them. The actual target size determines every
denominator and budget. Records are not weighted by their survey weights, and
the results concern the finite empirical pools rather than national
populations.

### A2. Harmonization and the deprivation count

Water, sanitation and fuel come from the variables `hv201`, `hv205` with
`hv225`, and `hv226` of the DHS household member recode (PR), and education
from `hv109`. BMI uses `sbbm` in Bangladesh and `ha40` or `hb40` elsewhere,
divided by 100. Values at or above 6000 before division
are treated as missing. BMI categories are <18.5, 18.5 to <25, 25 to <30, and
at least 30 kg/m². Work comes from `sb308` of the Bangladesh PR file, or from `v717` and
`mv717` of the women's (IR) and men's (MR) recodes, joined to PR records by
cluster, household and person identifiers.
Age and sex use `hv105` and `hv104`, the age in years and the sex (male or
female) of each household member as recorded in the household roster.

Most categories are harmonized from value labels because numeric codes differ
across surveys. Education uses the standardized `hv109` codes 0 to 5. The
Bangladesh work mapping uses the explicit map given below. All unmatched
values remain missing and are excluded from the pool of complete records. The
code archive contains the matching rules in order and the stratum membership.

| Input | Harmonized categories in generator order | Count adds one for |
|---|---|---|
| Water | on premises, other improved, unimproved, surface | unimproved or surface |
| Sanitation | improved private, improved shared, unimproved, open defecation | shared, unimproved or open |
| Fuel | clean, kerosene or coal, charcoal, biomass | any category except clean |
| Education | higher, secondary, incomplete secondary, primary, incomplete primary, none | incomplete primary or none |
| Work | not working, services, manual, agriculture | manual or agriculture |
| BMI | underweight, normal, overweight, obese | underweight or obese |

Water labels for rivers, dams, lakes, ponds, streams, canals and irrigation map
to surface. Unprotected, open or dug wells and “other” map to unimproved. Piped
water into a dwelling, yard or plot maps to on premises, and the remaining
improved supplies include protected sources, tubewells, boreholes, rainwater,
bottled water, tankers and carts. Sanitation labels for no facility, bush or
field map to open, and pits without slabs, hanging or bucket toilets,
unspecified discharge and other unimproved facilities map to unimproved.
Improved facilities are divided by reported sharing. Unknown sharing is not
excluded separately.

Clean fuel includes electricity, solar, LPG, natural gas, biogas and alcohol
fuels. Kerosene, coal and related liquid fossil fuels share one rung, and
biomass includes wood, crop residues and dung. “No food cooked” and “other” fuel records
are excluded. Work label matching gives agriculture precedence over service or
manual labels. In Bangladesh, `sb308` codes 0, 61 and 62 map to not working,
11, 41, 51 and 52 to services, 16, 21, 22, 23 and 31 to manual, and 12, 13, 14
and 15 to agriculture.

The count is the sum of six binary indicators, ranging from 0 to 6 before the
shared orientation step. Water, sanitation, fuel and schooling adapt the global
MPI definitions (Alkire, Kanagaratnam and Suppa, 2018), rather than reproducing
the full household MPI. Travel time for water is not used, education is
individual, and the combined kerosene and coal rung is counted as deprived. BMI
uses the WHO category thresholds (WHO, 2000), with underweight and obesity
flagged but overweight unflagged. The work indicator is a study design choice.
The score is therefore a prespecified count over six domains, not a validated MPI or
an externally validated measure of healthcare need. None of its six indicators
is filtered. The global MPI itself counts a household without the eligible
population for an indicator, such as one with no child of school age for
school attendance, as not deprived (Alkire and Santos 2014), which is the skip
coded "no" of the main text.

### A3. Imposed mechanisms

Let \(X_{ij}\) denote the six categories for person \(i\) in synthetic cluster
\(j\). For each variable other than BMI, with \(K\) canonical categories in the order
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
al. (2021). The .05 ICC is a residual value on the latent scale, conditional
on \(V\), age and sex, and set by calibration. Observed disease clustering is
higher because \(V\) is itself clustered. The sex and vulnerability
coefficients are design choices. An intercept for each population targets
roughly 25% disease prevalence.

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
was calibrated against a contact proxy in the Nepal survey. Its value and the
design correlation of −.25 are not direct estimates of these latent mechanisms.
The individual term \(.3V_{ij}\) pushes the other way, so the reference design
is nearly neutral to access at cluster level. In diagnostic replicates,
generated separately from the main simulation with the same generator (25 per
scenario, pooled stratum Q3 and Bangladesh 18 to 95 at \(r^2=.6\)), the
Spearman correlation between cluster contact rate and cluster mean \(V\) was
−.04 to −.05 in Q3 and −.02 in Bangladesh 18 to 95 at every contact level.
Under the diagnosis filter the cluster eligibility rate correlated .07 to .13
with cluster need at 20% contact and .22 to .28 at 80%. The variation with an
access correlation of −.5 represents the inverse care law. These diagnostic
values are not estimates from the main simulation.

Let \(z(a)\) standardize age within the draw. Diagnosis is

\[
q_{ij}=\begin{cases}
\operatorname{expit}\{\operatorname{logit}(.50)+\delta z(a_{ij})\},&Y_{ij}=1,\\
.05,&Y_{ij}=0,
\end{cases}\qquad
R_{ij,s}=C_{ij,s}1(U_{R,ij}<q_{ij}),
\]

where \(\delta\) is 0 or .5. The direction of the age association is motivated
by Foti and Hardy (2025), and its magnitude is a design choice. The .50 and .05
values are design anchors, not estimates of diagnosis conditional on contact
from studies of population awareness or of the validity of self reports.

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
design choices. BDHS 2017-18 cannot estimate that contrast because it
does not ask undiagnosed people about advice. There is no ordinary item nonresponse
among eligible respondents. Contact, diagnosis and advice collapse different
survey time windows into one simulated period. No intervention effect on disease
or vulnerability is generated.

Figure S2 shows the imposed structure. It omits independent noise terms and
empirical correlations among donor covariates. The arrow from vulnerability to
cluster access supply represents the imposed inverse care association, and the
two dashed arrows into eligibility are the alternative filters.

**Figure S2. The imposed structure of the simulated data.** The six covariates,
age and sex are resampled from the surveys, and every other node is simulated. Dashed
arrows mark the two alternative filters that decide who is asked about advice.
File: `JCE/figures/figS2.pdf`.

The outcome intercepts for Q1, Q2, Q3, Q4, Bangladesh adults, India, Nepal,
Albania and Bangladesh ages 18 to 49 are, respectively, −1.4018739424817,
−1.3874000585407, −1.3979337882166, −1.3897578168946, −1.5156911575790,
−1.3925435620646, −1.3550533785203, −1.3936592700210 and −1.3962307474406.

The variations each replace one parameter. They are diagnosis given disease
.20 or .80, diagnosis without disease .15, access correlation 0 or −.5, access
coefficient 0, two or eight advice items, all six loadings .65, nominal
development size 500, and 20 imputations. With two advice items the first two
baseline probabilities are used, and with eight the four are repeated. The low
and high diagnosis settings are also run under the diagnosis filter. The other
variations and all mechanism checks use the contact filter at the reference
configuration of pooled stratum Q3 (\(r^2=.6,\delta=0\)).

The mechanism checks set both diagnosis probabilities to .1625, replace advice
by four proxies observed for everyone with probability
\(\operatorname{expit}(1.2V)\), increase the advice loading from .25 to 1,
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
two filters share the same populations and uniforms. Separate random number
streams serve generation, development fitting and target scoring. Uniforms for
individual and cluster budget ties are drawn before any method is fitted and
do not depend on truth, disease or diagnosis.

MCA uses one dimension. Each development score is split at its unique empirical
quintile boundaries. If disease prevalence in the top tier is below that in the
bottom tier, its sign is reversed and the cuts are recomputed. This orientation
also applies to the count, oracle and supervised reference. Scores with ties
can produce fewer than five tiers. Target scores use the development orientation
and cut points, and values outside the development range are clamped to its
bounds for tier assignment. Budget rankings use the projected score before clamping.

Each imputer fits one logistic model per item among eligible respondents using
the six inputs, age and sex. Constant terms are omitted. For item prediction,
an unseen covariate category is replaced by the modal category among the
eligible respondents used in fitting, and the number of replacements is recorded. Five sets of
Bernoulli responses are drawn from fixed fitted probabilities, keeping observed
responses. Parameter uncertainty and residual dependence between items are not
propagated. Each completed development dataset yields an MCA axis, oriented and standardized on its development coordinates, and each person's five standardized scores are averaged into one before allocation. The variant with
target models refits the item models in the target, and both variants keep the
development MCA axes.

Under the diagnosis filter the development item models rest on few people. In the diagnostic replicates of A3 each was fitted on 83 to 86 diagnosed respondents. With 19 to 21 parameters and 19 to 21 people in the smaller outcome class, that is about one event per parameter, and fitted probabilities of 0 or 1 occurred. Under the contact filter the models had about 490 respondents and about three events per parameter.
In the main simulation, against dropping the items under the diagnosis filter,
the mean differences in cluster rank correlation were −.087 for development
models and −.055 for target models, the medians −.028 and −.009, and 23.9 and
14.4 percent of replicates trailed by more than .10. The mean deficit comes from
a tail of replicates with large deficits. The median replicate with target
models was within .01 of dropping.

Disease measures are the odds ratios of the top tier in development and in the
target, the development slopes within and between clusters, and the target AUCs
of the index and of the model of tier, age and sex fitted in development. Each
is evaluated on the people its method scores. Regression results flagged for
separation are excluded, as the specification of the performance measures
states (A5, Table S6).

The supervised reference is logistic regression on all categorical inputs and
the advice items with the skip coded 0, without age and sex in its score. Its
development scores use five folds separated by synthetic cluster, and a model
fitted on the whole development sample scores the target. Age and sex enter the downstream tier validation models
for every method. MCA and supervised projection leave target rows unscored if
they contain categories absent from the relevant development fit. Subset MCA excludes the categories of those not asked but keeps the masses of the others from the full data. Its target rows are placed by the transition formula, because the predict.MCA function of FactoMineR ignores the exclusion.

Individual budgets select the people with the highest scores, up to
\(\lfloor bN\rfloor\), for \(b=.05,.10,.20\). Unscoreable people are not
selected. Cluster budgets select up to \(\lfloor .10J\rfloor\) clusters using
the mean score of scored residents. Clusters with no scored resident are never
selected, and in the rank correlation over the whole target they tie below every
scored cluster. Slots remain
unfilled when too few people or clusters can be scored. Equal numbers of
clusters need not mean equal numbers of residents or equal delivery costs.
Table S5 reports the share of budget slots used and the residents covered.

The most vulnerable fifth is the top fifth of *V* in the whole target, fixed
before any method is evaluated. Undiagnosed disease means \(Y=1,R=0\). Cluster
truth uses all residents, including unscored residents. Inclusion error is the
fraction of selected clusters below median truth, whereas precision is the
fraction belonging to the true top decile. The share of the most vulnerable
fifth selected is the number of them selected divided by their number in the
whole target, \(n_H\). The corresponding shortfall per 1,000 is
\(1000[\min(\lfloor .10N\rfloor,n_H)-n_{H,selected}]/N\).

### A5. Prespecification, separation exclusions and a sensitivity analysis

The protocol was developed after exploratory pilot runs on 17 to 19 September
2026. The amended design, the performance measures, the decision rule and its
proposed margins were fixed on 20 September. A separate pilot of 100 replicates
on 21 September sized the main simulation, which ran on 23 to 25 September
2026, and the proposed margins were confirmed unchanged after it. None of these
steps is a preregistration.

The sizing pilot projected 124 of the 132 targeting contrasts to be decided at
1,000 replicates, with no further gain at 2,000. The projection used the
pilot's paired spreads and held its mean differences fixed. Pilot results are
not pooled with the main simulation, which used different random numbers. In
the main simulation 122 contrasts were decided, including one equivalence.

The specification of the performance measures excludes a disease measure when
its fit is flagged for separation. The prespecified summary did not apply these
exclusions to the correlations of Supplementary Figure S1, so they were applied
in a separate analysis that leaves the targeting results and the 132 verdicts
unchanged. Separation in development masks the development odds ratio of the
top tier and the slopes within and between clusters. Separation in the target
masks the target odds ratio of the top tier and the two target AUCs. The flag is
any fitted probability of the tier model below \(10^{-8}\) or above
\(1-10^{-8}\), an operational rule rather than a formal proof of separation.
Table S6 gives the exclusions, and the code archive keeps the medians computed
without them.

Within each core scenario, correlations pool the ten methods and 1,000
replicates, using the values available for each pair after exclusions.
Inclusion error and shortfall are reversed in sign. The 54 correlations include
dependent and algebraically redundant end points, so they are descriptive
comparisons, not independent hypothesis tests. Supplementary Figure S1
summarizes each pairing by its median and range across 72 scenarios, and its
bars are not confidence intervals.

A sensitivity analysis added after the main simulation first averages the valid
results of each method within each scenario, then correlates the ten method
means. Table S7 reports the six disease measures against the three principal
targeting measures. It is not part of the prespecified decision rule.

### A6. Proofs and verification of the theory

The theory, its tests and missing fuzzy average were derived after the main
simulation and are post hoc. They change no prespecified verdict.

**Notation.** There are *n* respondents and *Q* = *P* + *k* categorical items,
the substantive items \(X_1,\dots,X_P\) and the filtered items
\(A_1,\dots,A_k\). The filtered items are all asked of the same eligible
respondents, \(G_i=1\), with \(p=\bar G\). A coding turns item
*q* into indicator columns \(z_j\), \(j\in J_q\). Eligible respondents receive
the crisp indicator of their answer, and a skipped respondent receives a profile
\(v_i\) over \(J_q\), crisp or fuzzy. Subset analysis keeps a set \(S\) of
analysed columns with their masses from the full data. Write \(n_j=\sum_i z_{ij}\),
\(f_j=n_j/n\), \(f_{j\mid G}\) for the answer distribution among the eligible
(zero for a column the eligible never use, such as a skip category), and
\(\bar v_j\) for the mean profile given to the skipped, so that
\(f_j=pf_{j\mid G}+(1-p)\bar v_j\). For a centered score \(h\) with
\(n^{-1}\sum_ih_i^2=1\), the squared correlation ratio of item *q* is

\[
\eta_q^2(h)=\frac1n\sum_{j\in J_q\cap S}\frac{\left(\sum_iz_{ij}h_i\right)^2}{n_j},
\]

which for a crisp item is \(\sum_jf_jm_j(h)^2\), with \(m_j\) the mean of \(h\)
in category *j*. For fuzzy codings and subset analyses it is the item's
contribution to the criterion, and we keep the name. The eligibility score is \(s_G=(G-p)/\sqrt{p(1-p)}\).

**Lemma A6.0 (the MCA criterion).** With
\(B=Q^{-1}\sum_q\sum_{j\in J_q\cap S}(z_j-f_j\mathbf 1)(z_j-f_j\mathbf 1)^\top/n_j\),
every centered score with unit variance satisfies
\(h^\top Bh/n=\lambda(h)=Q^{-1}\sum_q\eta_q^2(h)\). The first axis maximizes
\(\lambda(h)\), and \(\lambda_1=\max_h\lambda(h)\). *Proof.* Column *j* of the
standardized residual matrix \(T\) is \((z_j-f_j\mathbf 1)/\sqrt{Qn_j}\), so \(B=TT^\top\), and
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
the mean square contingency between eligibility and the coded item.

**Corollary A6.1 (the codings).** For a binary item that the eligible answer
"yes" with probability \(\pi\):

(a) Missing single, the skip as its own category, with \(\bar v=(0,0,1)\) on (no, yes,
skip), gives \(\eta_q^2(s_G)=p(1-p)\,[(1-\pi)/p+\pi/p+1/(1-p)]=1\), whatever the number of answer
categories.

(b) Subset MCA, with the skip column excluded, gives
\(\eta_q^2(s_G)=p(1-p)\,[(1-\pi)/p+\pi/p]=1-p\). Subset MCA gives the same solution as
missing passive modified margin (Josse et al. 2012), so modified margin keeps
\(1-p\) as well. The "multiple" variant of subset MCA for incomplete data
(Nienkemper-Swanepoel et al. 2024), with one missing category per cell, has the
same analysed columns and the same value. In a simulation, Hendry (2015) found
that subset correspondence analysis performed well provided less than 30% of any
variable was missing, which is consistent with a residual of \(1-p\) per item.

(c) The skip coded "no", with \(\bar v=(1,0)\), has \(f_{\text{no}}=1-p\pi\) and
\(f_{\text{yes}}=p\pi\), so
\(\eta_q^2(s_G)=p(1-p)\,[\pi^2/(1-p\pi)+\pi/p]=\pi(1-p)/(1-p\pi)\). It is below the
subset value exactly when \(\pi<1/(1+p)\).

(d) Missing fuzzy average, with \(\bar v=f_{\cdot\mid G}\), gives \(\eta_q^2(s_G)=0\). Since every
term of Theorem 1 is nonnegative, \(\eta_q^2(s_G)=0\) if and only if
\(\bar v_j=f_{j\mid G}\) for every analysed column with \(f_j>0\). Among codings
that analyse every column, it is the only mean profile that carries no filter.

(e) Under imputation, \(\bar v\) is the mean imputed profile, and \(\eta_q^2(s_G)\) is near 0
when that profile is near the answers of those asked.

**Corollary A6.2 (a bound under missing fuzzy average).** Under missing fuzzy average,
\(\eta_q^2(h)\le p\cdot n^{-1}\sum_{i:G_i=1}h_i^2\le p\) for every centered
\(h\) with unit variance. The filtered block analysed alone therefore has
\(\lambda_A\le p\), and \(\kappa\le kp/(P\lambda_X)\). *Proof.* The skipped rows of
each column equal \(f_j\), so with \(S_j\) the sum of \(h\) over the eligible in
category *j* and \(H=\sum_jS_j\), \(\sum_iz_{ij}h_i=S_j-f_jH\) because \(h\) is
centered, and \(\eta_q^2=n^{-2}\big[\sum_jS_j^2/f_j-H^2\big]\). The eligible in
category *j* number \(npf_j\), so by the Cauchy-Schwarz inequality
\(S_j^2\le npf_j\sum_{i\in E_j}h_i^2\), and summing over *j* gives the bound. ∎ The bound is reached by
any score that is a centered function of the item's answer among the eligible
and zero for the skipped. Under missing fuzzy average the skipped sit at each
item's centroid and add no inertia (Josse et al. 2012). The corollary bounds
the inertia it leaves to the items, and A6.3 gives \(\lambda_A\) exactly.

**Proposition A6.3 (a lower bound, and the eigenvalues of the filtered block).**
\(\lambda_1\ge\lambda(s_G)=Q^{-1}\sum_{q}\eta_q^2(s_G)\ge
Q^{-1}\sum_{q\le k}\eta_q^2(s_G)\), indexing the filtered items first. Let
\(\lambda_E\) be the first eigenvalue of the MCA of the filtered items among the
eligible alone. Analysed alone, the coded filtered block has the eligibility
score \(s_G\) as an eigenvector, with eigenvalue 1 under missing single,
\(1-p\) under subset MCA and 0 under missing fuzzy average. Its other
eigenvalues are those of the analysis of the eligible alone, multiplied by \(p\) under
missing fuzzy average, so \(\lambda_A\) is 1, \(\max(1-p,\lambda_E)\) and
\(p\lambda_E\). *Proof.* For an answer column,
\(z_j-f_j\mathbf 1=u_j+f_{j\mid G}(G-p\mathbf 1)\) under missing single and
subset MCA, where \(u_j\) equals \(z_{ij}-f_{j\mid G}\) for the eligible and 0
for the skipped. Under missing fuzzy average \(z_j-f_j\mathbf 1=u_j\). The
\(u_j\) of an item sum to zero and each is orthogonal to \(G-p\mathbf 1\), so the
cross terms vanish. The block's matrix is then the matrix of the eligible alone,
multiplied by \(p\) under missing fuzzy average, plus
\((G-p\mathbf 1)(G-p\mathbf 1)^\top\) divided by \(np\) under subset MCA, and by
\(np(1-p)\) once the skip column is added under missing single. ∎ Under missing
single, \(\lambda_1\ge k/Q\), and the filter is an exact dimension. In the core
design \(k/Q=.40\), and in diagnostic replicates the first eigenvalue was .404
to .407, with an axis correlating .995 to .998 with eligibility. In the 51,000
development samples of the main simulation the block eigenvalue was 1 under
missing single in every sample. Under subset MCA it equalled \(1-p\) in every sample under the
diagnosis filter and in 83% under the contact filter, and was larger in the
rest, where \(\lambda_E\) exceeds \(1-p\).

**Theorem 2 (the joint analysis as an update of rank one).** Let \(B_X\) be the
matrix of Lemma A6.0 for the substantive items analysed alone, \(B_E\) that of
the filtered items analysed among the eligible alone, with zeros for the
skipped, and \(u=s_G/\sqrt n\), a unit vector. For any data,

\[
QB=PB_X+c\,kB_E+k\theta\,uu^\top,
\]

with \((c,\theta)=(1,1)\) under missing single, \((1,1-p)\) under subset MCA
and \((p,0)\) under missing fuzzy average. Dropping the items leaves \(PB_X\).
*Proof.* \(QB\) is the sum over items of
\(\sum_{j\in J_q\cap S}(z_j-f_j\mathbf 1)(z_j-f_j\mathbf 1)^\top/n_j\). By the
proof of A6.3 this sum is, for a filtered item, its counterpart among the
eligible alone multiplied by \(c\), plus \(\theta uu^\top\), because
\((G-p\mathbf 1)(G-p\mathbf 1)^\top=np(1-p)\,uu^\top\). ∎ The three codings
therefore give one analysis, apart from a term along eligibility of weight
\(k\theta\) and the factor \(c\) on the answers of those asked. Since \(u\) is
orthogonal to every column of \(B_E\), each filtered item contributes exactly
\(\theta\) to the criterion at eligibility, which is Theorem 1 for these
codings.

Write \(B_0=(PB_X+c\,kB_E)/Q\), with eigenvalues \(\mu_1\ge\mu_2\ge\dots\ge0\)
and unit eigenvectors \(w_i\), \(\rho=k\theta/Q\) and \(c_i=u^\top w_i\). Then
\(B=B_0+\rho uu^\top\) is a modification of \(B_0\) of rank one (Golub 1973). Its
eigenvalues are the roots of

\[
1=\rho\sum_i\frac{c_i^2}{\lambda-\mu_i},
\]

together with every \(\mu_i\) whose \(c_i\) is zero. When \(c_1\ne0\) the first
eigenvalue is the largest root. It lies between
\(\max(\mu_1,\rho+u^\top B_0u)\) and \(\mu_1+\rho\), and the squared
correlation of the first axis with eligibility is
\(\big[\rho^2\sum_ic_i^2/(\lambda_1-\mu_i)^2\big]^{-1}\). Here
\(u^\top B_0u=Q^{-1}\sum_{q>k}\eta_q^2(s_G)\), the substantive items' own
squared correlation ratios with eligibility. *Proof.* If \(Bv=\lambda v\) and
\(\lambda\) is not an eigenvalue of \(B_0\), then
\(v=\rho(u^\top v)(\lambda I-B_0)^{-1}u\). Multiplying by \(u^\top\) gives the
equation, and normalizing \(v\) gives the correlation. The bounds are Weyl's
inequality and the value of the criterion at \(u\). ∎

Two consequences follow. First, if the substantive items are uncorrelated with
eligibility in the sample, in that every centered indicator of a substantive
item is orthogonal to \(G-p\mathbf 1\), then \(c_i=0\) whenever \(\mu_i>0\).
The eigenvalues of \(B\) are then those of \(B_0\) together with \(\rho\), and
eligibility is the first axis exactly when \(k\theta>Q\mu_1\). This asks less
than the independence of the two blocks used in A6.4, since the answers of
those asked may be associated with the substantive items, and \(Q\mu_1\) is at
least \(\max(P\lambda_X,c\,k\lambda_E)\). Second, three quantities compete for
the first axis, \(P\lambda_X\) for the substantive items, \(c\,k\lambda_E\) for
the answers of those asked and \(k\theta\) for the filter. The inertia ratio of
A6.4 compares the first with the larger of the other two. It says when the
filtered block leads, and the filter leads only if \(\theta\) also exceeds
\(c\lambda_E\).

The skip coded "no" has no such split. For a binary item its matrix is the
projector onto \(a_q=u_q+\pi_q(G-p\mathbf 1)\), where \(u_q\) holds the centered
answers of those asked and zeros for the skipped. Each item is then a single
direction that mixes the answers with eligibility, at a squared cosine with
eligibility equal to the value of Corollary A6.1(c).

**Proposition A6.4 (block dominance, the inertia ratio and instability).** Let
\(h_X\) and \(h_A\) be the first axes of the substantive block and of the coded
filtered block analysed alone, with eigenvalues \(\lambda_X\) and \(\lambda_A\).
On the span of \(\{h_X,h_A\}\) the criterion is the generalized eigenproblem
\(Mw=\mu Cw\), with \(M_{ab}=Q^{-1}\sum_q\eta_q(h_a,h_b)\), the cross form of the
squared correlation ratio, and \(C\) the Gram matrix of the two scores. Its
leading value is a Rayleigh-Ritz lower bound for \(\lambda_1\), exact when the
first axis lies in that span. The diagonal entries are
\(Q^{-1}[P\lambda_X+\sum_{q\le k}\eta_q^2(h_X)]\) and
\(Q^{-1}[k\lambda_A+\sum_{q>k}\eta_q^2(h_A)]\). The leading term of each diagonal entry, \(P\lambda_X\) or \(k\lambda_A\), is the block's number of items times its maximum axial inertia, the first eigenvalue of its separate analysis. It is what the block contributes along its own first axis in the joint analysis. When the blocks are weakly associated the cross terms are small and the leading direction follows the larger diagonal, so the filtered block dominates as the inertia ratio \(\kappa=k\lambda_A/(P\lambda_X)\) passes 1. \(B\) is the sum of the two blocks' matrices weighted by \(P/Q\) and \(k/Q\) for any data. If the blocks are independent in the sample, in that every centered indicator of one is orthogonal to every centered indicator of the other, the two weighted matrices have orthogonal ranges, the eigenvalues of \(B\) are theirs so weighted, and \(\kappa=1\) is the exact threshold. Multiple factor analysis weights each group by the inverse of its first
eigenvalue so that no group dominates (Escofier and Pagès 1994), whereas
\(\kappa\) measures the dominance an unweighted filtered block exerts. The weighting
equalizes the blocks and does not remove the filter. Under missing single the
filtered block's first eigenvalue is 1 with eligibility as its axis (A6.3), so
the weighted analysis gives eligibility and the substantive axis the same
leading weight. In 300 diagnostic development samples (A6.9) eligibility explained
.64 of the first weighted axis on average (range .22 to .98) and the
substantive axis .39. By A6.3,
\(\kappa\) equals \(k/(P\lambda_X)\) under missing single,
\(k\max(1-p,\lambda_E)/(P\lambda_X)\) under subset MCA and
\(kp\lambda_E/(P\lambda_X)\) under missing fuzzy average. Theorem 1 gives the
criterion at eligibility, \(\sum_{q\le k}\eta_q^2(s_G)\), which is \(k\),
\(k(1-p)\) and 0 under these codings. It equals \(k\lambda_A\) under missing
single, and under subset MCA whenever \(1-p\ge\lambda_E\).
For items coded in two categories, as under the skip coded "no",
\(k\lambda_A\) is the largest eigenvalue of their correlation
matrix, at least \(1+(k-1)\bar r\) with \(\bar r\) their mean correlation. In diagnostic replicates with
\(|\kappa-1|>.25\), the span of \(\{h_X,h_A\}\) explained a median .997 of
the first axis (5th percentile .954).

By the Davis-Kahan theorem in the form of Yu, Wang and Samworth (2015), the
sine of the angle between a sample axis and its population counterpart is at
most twice the sampling perturbation of \(T^\top T\), in operator norm, divided
by the population eigengap \(\lambda_1-\lambda_2\). \(T^\top T\) is the matrix
over the analysed categories with the same nonzero eigenvalues as
\(B=TT^\top\), and unlike \(B\) its size does not grow with \(n\). The theorem is
an upper bound. In the reduction to two directions the gap closes as \(\kappa\to1\),
where the bound no longer limits how far a sample axis can turn between the
blocks, and A6.9 shows that it does turn. Under subset MCA a second crossing
lies inside the filtered block, between \(1-p\) and \(\lambda_E\) (Theorem 2).

**Proposition A6.5 (cluster ranking by a weighted score, and the spread of eligibility).** Write
the deployed score as \(s=\alpha h_X+\beta h_A\) plus a remainder. At cluster
level in the target, with standard deviations between clusters
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
a property of the target. If cluster eligibility rates vary around *p* with
variance \(\sigma_L^2\) on the logit scale, the delta method and binomial sampling
within clusters of mean size \(\bar n\) put the variance of the eligibility rate
between clusters at about \([p(1-p)]^2\sigma_L^2+p(1-p)/\bar n\), largest at
\(p=\tfrac12\). With \(h_A\) close to eligibility and scaled in development,
\(\sigma_A\) is proportional to its square root. Under the contact filter the target's *p* is its contact rate,
so the filter's pull peaks at 50% contact. Under the diagnosis filter *p* stays
below one half, and both \(\sigma_A\) and \(\rho_A\) rise with contact. With few
people eligible, \(\sigma_A\) is small and the substantive part still orders
many clusters, so the index can rank them better than the eligibility rate
does (A6.11). In diagnostic replicates, the reduction to two directions,
computed in development, predicted cluster targeting in the target with
r = .94 and mean absolute error .039 where \(|\kappa-1|>.25\).

**Proposition A6.6 (a budget below the share asked).** If the index ranks every eligible
person above every skipped person and the budget fraction is \(b\le p\), every
selected person is eligible, and the selected share of any group outside the
filter is zero. If \(b>p\), at most \((b-p)/b\) of the slots go outside it.
*Proof.* Counting. ∎ The skip coded "no" meets the premise only approximately, since
an eligible respondent who answered "no" to every item has the same coded
items as a skipped one, which is why it selected .00001 of the people with
undiagnosed disease on average, and at most .00004 in any configuration, rather
than none.

Missing fuzzy average is neutral in mean, not in spread. The skipped sit at
each item's centroid, so the filtered items add nothing to their scores, while
those asked keep the items' variance and reach both ends of the score. When
the filtered block carries weight on the axis, those asked therefore fill more
of a tight budget than their share. With eight advice items under the contact
filter, where its inertia ratio averaged .89 in diagnostic replicates, the
share of the most vulnerable fifth selected at the 10% budget fell from .314 to
.286 and .272 as contact rose from 20 to 50 and 80%, against .322 for dropping
the items. In the core design, with inertia ratios of about .45 to .54 under the
contact filter and .05 to .06 under the diagnosis filter, it stayed within .012
of dropping on this share (Table S10).

**Proposition A6.7 (structural zeros and identification).** In the generator,
factual advice is \(A^*_{ijk,s}=C_{ij,s}\cdot1[\cdot]\) (A3). Under the contact
filter a skipped respondent has \(A^*=0\) with certainty, so the skip coded "no"
reproduces factual advice exactly. Any harm comes from the item itself, since at a given need a person's chance of having been advised depends on their access. Under the diagnosis filter, contacted
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
Corollary A6.1(c). The DHS wealth index codes answers that do not apply as 0 in
this way, for instance whether a toilet is shared (Rutstein, no date). The bound
of A6.3 and Propositions A6.4 to A6.6 carry over with eigenvalues of correlation
matrices. Since such an eigenvalue already sums over columns, the inertia ratio
is \(\lambda_A/\lambda_X\). Corollary A6.2 does not carry over. Replacing
the skip by the share of those asked answering yes (mean imputation) leaves the
column uncorrelated with eligibility, but the imputed columns keep the
correlation matrix of those asked, so their block is not shrunk by *p*. Its
inertia ratio was .51 under the diagnosis filter (Table S14), against .06 for
missing fuzzy average at 50% development contact (Table S13). A7 tests this
extension by simulation.

**A6.9 The theory in the main simulation.** Every replicate of the main
simulation was regenerated from its random seeds, which reproduced the true
\(V\) and the rank correlation of dropping the items exactly in all 52
configurations. In each development sample we computed, for missing single,
subset MCA and the skip coded "no", the empirical squared correlation ratio with
eligibility and its value by Theorem 1, \(\lambda_X\), \(\lambda_A\),
\(\kappa\), the first two eigenvalues, and the squared correlations of the
index's first axis with \(h_A\) (block R²), with eligibility (eligibility R²)
and with \(h_X\). Missing fuzzy average is fitted as correspondence analysis of
the fuzzy indicator matrix. Development and target rows are placed by the
transition formula on the development column standard coordinates, which
reproduced the fitted row coordinates to \(5\times10^{-14}\). In the target,
the skipped receive the answer distribution of the target's own eligible
respondents, so no model is fitted, and rows with a category the development
sample never saw are left unscored, as for every MCA method. Missing fuzzy
average was evaluated on the same performance measures as the prespecified
methods, in all 1,000 replicates of the 52 configurations, 51,000 of which have
filtered items. Table 2, Figure 1, Table S10 and the figures below come from
this analysis.

- *Theorem 1.* Over the 153,000 combinations of development sample and
  coding, the largest difference between the empirical value and the formula was
  \(6\times10^{-15}\). Missing fuzzy average is 0 by Corollary A6.1(d) and was
  not recomputed. Theorem 1 is an identity, so this checks that each coding was
  applied as defined. The skip coded "no" carried less filter than subset MCA in
  51,000 of 51,000 replicates.
- *Inertia ratio.* Over 153 pairs of configuration and coding, its Spearman
  correlation with the loss against dropping the items was .83 (.88 for the block R², .79
  replicate by replicate). For each coding separately it was .74 for subset MCA, .86 for the
  skip coded "no" and .23 for missing single, whose inertia ratio was above 1.2
  in every configuration. Above \(\kappa=1\) the block took on average .99 of
  the axis under missing single, .74 under subset MCA and .86 under
  the skip coded "no". The two pairs below 1 took .15 and .24 and lost −.005 and .000.
- *The inertia ratio in single samples.* Take eligibility explaining more than
  half the first axis as the outcome. In the development samples of the main simulation
  \(\kappa>1\) agreed with it in every sample under missing single and under
  the diagnosis filter. Under subset MCA and the contact filter (33,000
  samples) the outcome occurred in .505 and \(\kappa>1\) in .827, with
  sensitivity .998, specificity .346 and accuracy .675. Requiring also that the
  block's own first axis be eligibility (block eigenvalue \(1-p\), A6.3) raised
  the accuracy to .785. Where \(\kappa>1\) and the block's axis was the answers
  of those asked, the block explained .72 of the first axis and eligibility
  .09. In the core configurations under the contact filter, \(P\lambda_X\), \(k(1-p)\) and
  \(k\lambda_A\) averaged 1.83, 1.99 and 2.00, so three directions were nearly
  tied (Theorem 2). Under the skip coded "no" and the contact filter \(\kappa\)
  exceeded 1 in .863 of samples and eligibility never explained half the axis.
  That coding's losses come through the advice answers (A6.7). The inertia
  ratio says when the filtered block leads. Whether the filter leads is read
  from the fitted index (Table 5, step 3) or from Theorem 2.
- *Theorem 2.* We built the indicator matrices directly in 300 diagnostic
  development samples, 50 in each of six scenarios (the two core populations
  under each filter, and two and eight advice items under the contact filter). Under missing single and subset MCA the two sides of the
  decomposition differed by at most \(1\times10^{-16}\). The largest root of
  the secular equation equalled the first eigenvalue to \(1\times10^{-15}\) and
  FactoMineR's to \(2\times10^{-14}\), and the squared correlation of the first
  axis with eligibility equalled its predicted value to \(2\times10^{-13}\).
  These are identities and check the algebra and its implementation. Eligibility
  explained more than half the first axis in .825 of the 600 pairs of sample
  and coding. The rule \(k\theta>Q\mu_1\), exact when the substantive items are
  unrelated to eligibility, agreed with that in .975 and \(\kappa>1\) in .94.
  For subset MCA with four items under the contact filter (100 samples), where
  \(P\lambda_X\) was 1.67 and 1.97 in the two populations, \(k(1-p)\) 1.98 and
  2.00 and \(k\lambda_E\) 1.82 and 1.81, the agreement was .85 and .64. The
  multiple factor analysis weighting of A6.4 was checked against the MFA
  function of FactoMineR on one sample.
- *Instability.* The standard deviation of targeting across replicates within a configuration was .12 for \(|\kappa-1|\) up to .15, .16 from .15 to .30, .08 from .30 to .60 and .06 above .60, with mean eigengaps .036, .038, .099 and .195 in the same bands. The Spearman
  correlation between \(|\kappa-1|\) and the standard deviation was −.64. The
  pattern is that of subset MCA, whose standard deviation was .21 within .3 of
  \(\kappa=1\) and .06 beyond. For the skip coded "no", which carries little
  filter near \(\kappa=1\), it was .08 and .07.
- *Budget below the share asked.* In the 18 scenarios under the diagnosis
  filter whose diagnosed share exceeded the 10% budget, 17 at 80% contact and
  one at 50%, missing single and subset MCA selected none of the people with
  undiagnosed disease and the skip coded "no" at most .00004, against .17 for
  dropping the items. Below the budget they selected .03 to .12.
- *Corollary A6.2.* In diagnostic replicates (12 per scenario), the block
  eigenvalue under missing fuzzy average averaged at most .64 of *p* in the
  12 scenarios (by A6.3 this ratio is \(\lambda_E\)), and its inertia ratio
  averaged .40 to .63 under the contact filter (.89 with eight items) and .03 to
  .09 under the diagnosis filter.

**A6.10 A published simulation.** Josse et al. (2012, Table 3, rows 7 and 8)
report that modified margin "can really crash" when 60 of 100 individuals miss
variables 1 to 3 together and 60 miss variables 9 and 10 together, without a
stated reason. Theorem 1 gives modified margin a squared correlation ratio with the missingness pattern of .6 per item
there and missing fuzzy average zero. We rebuilt the design with 500
replications per scenario and reproduced all six published medians of the
modified RV coefficient within .01. The reproduced medians for
weak and strong structure are .30 and .45 for missing single, .42 and .71 for
modified margin (published .43 and .71), and .77 and .91 for missing fuzzy
average. The missingness indicators explained .50, .48 and .01 of the
three methods' first two axes under weak structure, and .50, .23 and .01 under
strong structure. The collapse is the missingness pattern dominating the first axes.

**A6.11 Cluster ranking under a dominant filter.** Data on individual people
were kept for two replicates of every configuration, 24 per filter in the core
design, and the figures below come from them. Under the contact
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

Three studies, designed after the theory, take up three questions the main
design left open. They ask how
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

Configurations with a development sample of 1,000 reuse the random numbers of
the reference configuration of the main simulation, so the people and the
uniforms are the same in every one of them and only the advice, or the
development sample's contact, changes. At loading .25 and 50% development
contact the configuration is the reference configuration itself, and the ten
prespecified methods reproduced its results exactly in all 1,000 replicates
under both filters. The larger development samples use their own random
numbers.

Every configuration scores the ten prespecified methods and missing fuzzy
average (A6.9), and four further kinds of index.

- *Imputation with FLIC item models*, under the diagnosis filter only, where
  the development item models rest on about 84 diagnosed respondents. Each
  item model's slopes are estimated by Firth's penalized likelihood and its
  intercept is then estimated again so that the mean predicted probability
  equals the observed rate (Puhr et al. 2017), with the R package logistf
  1.26.1. Development and target models are used as in the prespecified
  imputers.
- *Regularized iterative MCA* (Josse et al. 2012). The skipped cells of the
  indicator matrix are imputed by a regularized reconstruction of rank two
  (R package missMDA 1.23), and the index is the first axis of the MCA of the
  completed fuzzy indicator matrix. The target is imputed on its own data and
  placed by the transition formula. Because missMDA requires FactoMineR 2.17,
  this method was fitted with that version and every other method with 2.13.
  The two versions give MCA results that differ in the fifth or sixth decimal
  place. The filter diagnostics of this method were not computed.
- *Principal component indices* on 0/1 indicators of every input category, the
  score being the first principal component of their correlation matrix, as
  for asset indices. A filtered item enters as its yes indicator with the skip
  coded 0, as the same indicator with the skip replaced by the share of those
  asked answering yes (mean imputation, the counterpart of missing fuzzy
  average but without the shrinkage of Corollary A6.2, A6.8), as yes, no and
  skipped indicators, or not at all. A6.8 gives the
  closed forms checked in every development sample.

Of 822,000 method results, 18 failed, all from the two FLIC methods in 3 of
the 1,000 replicates with 20% development contact under the diagnosis filter.
There the development sample held 27 to 31 diagnosed respondents for 17 to 20
coefficients per item model, and the penalized fit stopped on an information
matrix that was not positive definite. No table uses those results.

- *Advice loading* (Table S11, Figure 4A). Under the contact filter the
  imputers overtook dropping the items at loadings of .16 (target models) and
  .24 (development models), against the base loading of .25. Under the diagnosis
  filter target models crossed at .67 and development models stayed behind at
  every loading (−.019 at 1). With FLIC item models the crossings were .47 and
  .81. Regularized iterative MCA led dropping by .003 to .041 under the contact
  filter and stayed within .002 of it under the diagnosis filter. Missing fuzzy
  average led by .002 to .032 and stayed within .001. As the loading grew,
  subset MCA rose from .342 to .833 and the skip coded "no" from .525 to .791
  under the contact filter. Under the diagnosis filter "no" rose from .396 to
  .517 and subset MCA stayed at .334 to .336.
- *Development size* (Table S12, Figure 4B). Under the diagnosis filter the
  imputers' mean deficits against dropping were .107 and .161 with 1,000
  development respondents (84 diagnosed, 1.4 events per parameter of the item
  models) and .012 and .029 with 2,000 (2.7). With 4,000 (5.4) target
  models led by .004 and development models trailed by .001. With 1,000, FLIC
  item models cut the deficits to .033 and .061 and the share of replicates
  trailing dropping by more than .10 from .293 and .470 to .101 and .181.
  Dropping the items itself rose with development size, from .859 to .882 and
  .894.
- *Development contact* (Table S13, Figure 4C). Missing single, subset MCA and
  the skip coded "no" ranked worst at 50% target contact in 17 of the 18
  combinations of filter, coding and development contact. The exception,
  subset MCA developed at 80% contact under the contact filter, fell from .848
  to .801 as target contact rose. The development sample's contact level, more than the target's, set how well these codings ranked clusters. Under the
  contact filter, as development contact went from 20% to 50% and 80%, subset
  MCA's inertia ratio fell from 1.93 to 1.20 and 1.09 and its filtered block
  still took .992, .706 and .742 of the first axis, but eligibility explained
  .992, .625 and .002 of it, and its cluster rank correlation, over the three target contact levels, was .007 to .042, .434 to .520 and .801 to .848. With 80% asked, the block keeps
  \(1-p=.2\) of the filter (Theorem 1). As this is below \(\lambda_E\), its axis
  follows the answers of those asked (A6.3). Missing single kept \(\kappa=2.42\), eligibility explained .996
  of its axis at every level, and it stayed below .05. Under the diagnosis
  filter the minimum at 50% target contact also appeared. There the spread of
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

### A8. Filtered items in published indices

We examined how indices built from questionnaire data treat items that are
asked only after a filter, in the documentation of widely used indices and
instruments and in open access papers that build an index by multiple
correspondence analysis. This is a scoping search, not a systematic review.

| Index and how it is built | Filtered items | Treatment of the skip | In the terms of this paper |
|---|---|---|---|
| DHS wealth index, principal components (Rutstein, no date) | toilet sharing, asked of households with a facility | "0 if not shared or question is not applicable" | "no" |
| Global Multidimensional Poverty Index, weighted count (Alkire and Santos 2014, Dotter and Klasen 2017) | school attendance, nutrition and child mortality, defined for households with an eligible member | counted as not deprived, and over a third of households had no child of school age in the three countries examined | "no" |
| Household Food Insecurity Access Scale, sum (Coates et al. 2007) | nine frequency questions, asked when the occurrence answer is yes | frequency coded 0 "for all cases where the answer to the corresponding occurrence question was 'no'" | "no" |
| Harmonized asset index, principal components (Varghese et al. 2021) | assets not asked in one survey wave | imputed as zero | "no" |
| Female Sexual Function Index, sum of domain scores (Meyer-Bahlburg and Dolezal 2007, Boehmer et al. 2012, Forbes et al. 2014) | twelve items with a response for no sexual activity and three for no attempt at intercourse | scored 0, the low point of the scale, and "zero scoring biases sexually inactive women towards greater measures of sexual dysfunction" | "no" |
| International Index of Erectile Function, sum (Yule et al. 2011, Forbes et al. 2014) | items with a response for no sexual activity | scored 0, below the least functional answer | "no" |
| SWPER women's empowerment index, principal components (Ewerling et al. 2017) | items asked only of women in a union, and age at first birth | restricted to women in a union, which left out 34% of women on average (6% to 62% by country), and hot deck imputation for women without children | complete case, imputation |
| Women's bargaining power index, MCA (Tchakounté Tchuimi and Kamga 2020) | items on decision making | 6,617 of 14,677 women, those not in a union, excluded | complete case |
| Maternal empowerment index, MCA (Mulamba 2026) | items on decision making | married mothers only | complete case |
| Antenatal care quality, count of services (Arsenault et al. 2018) | services received, asked of women who attended | reported among attenders, beside coverage | complete case |
| Antenatal care quality index, MCA (Mchenga et al. 2019) | services received, asked of women who attended | analysed among women who used antenatal care | complete case |
| Women's autonomy score, MCA (Rai et al. 2020) | items on decision making, in an analysis of all women aged 15 to 49 | not stated | not reported |
| Kansas City Cardiomyopathy Questionnaire, scaled means (Spertus and Jones 2015) | activities not done for reasons other than heart failure | scored as missing | items set aside |
| WHODAS 2.0, sum (Üstün et al. 2010) | four work or school items | "32 items being used if the respondent is not in gainful employment" | items dropped for the skipped |

For the Multidimensional Poverty Index the rule is in section 3.3 of the
working paper version of Alkire and Santos (OPHI Working Paper 59, 2013).
Questions on decision making in the Demographic and Health Surveys have been asked
only of currently married women since about 2004 (ICF, no date).

Three of the treatments compared in this paper are in use, the skip coded "no",
restriction to those asked and imputation. Two instruments instead score the
skipped on their remaining items. Coding the skip "no" is the rule in the two
most used development indices and in several clinical scores. It is a
structural zero where the skip implies the answer, as a condition that never
occurred has a frequency of zero, and benign where passing the filter is part
of the construct, as having a facility is part of wealth. It has been
criticized where neither holds, as for the Multidimensional Poverty Index and
the sexual function indices. Restriction to those asked is usual for
empowerment indices built from Demographic and Health Surveys. We found no
index that gives the skipped the answers of those asked, none that reports how
much of a filter its score carries, and no applied index that states it coded a
skip as its own category or set that category aside in a subset analysis. None
of these indices was audited here.

Europe PMC, searched on 12 November 2026, held 169 open access papers with full
text that name multiple correspondence analysis and a wealth, asset,
socioeconomic, composite, vulnerability, deprivation, empowerment, knowledge or
quality index in their methods. In 124 the analysis was used to build an index
or score. In 101 of these it was a household wealth or socioeconomic index. Of
the 23 others, at least four included items asked only of a subgroup, on
decision making or on the content of antenatal care. Three restricted the
analysis to those asked and one did not say how the others were scored (table).
The search covers open access papers only, and it cannot show how often a skip
is coded without comment.

### References cited only in the supplement

Arsenault C, Jordan K, Lee D, Dinsa G, Manzi F, Marchant T, et al. Equity in
antenatal care quality: an analysis of 91 national household surveys. Lancet
Glob Health 2018;6:e1186-95. https://doi.org/10.1016/S2214-109X(18)30389-9.

Coates J, Swindale A, Bilinsky P. Household Food Insecurity Access Scale (HFIAS)
for measurement of food access: indicator guide. Version 3. Washington, DC:
Food and Nutrition Technical Assistance Project; 2007.

Dotter C, Klasen S. The Multidimensional Poverty Index: achievements,
conceptual and empirical issues. Courant Research Centre Poverty, Equity and
Growth Discussion Paper 233. Göttingen: University of Göttingen; 2017.

Forbes MK, Baillie AJ, Schniering CA. Critical flaws in the Female Sexual
Function Index and the International Index of Erectile Function. J Sex Res
2014;51:485-91. https://doi.org/10.1080/00224499.2013.876607.

Hendry GM. The management of missing categorical data: comparison of multiple
imputation and subset correspondence analysis [PhD thesis]. Durban: University
of KwaZulu-Natal; 2015. http://hdl.handle.net/10413/15643.

ICF. Guide to DHS Statistics: participation in decision making. Rockville, MD:
The DHS Program; no date.
https://dhsprogram.com/data/Guide-to-DHS-Statistics/Participation_in_Decision_Making.htm
(accessed 12 November 2026).

Mchenga M, Burger R, von Fintel D. Examining the impact of WHO's Focused
Antenatal Care policy on early access, underutilisation and quality of
antenatal care services in Malawi: a retrospective study. BMC Health Serv Res
2019;19:295. https://doi.org/10.1186/s12913-019-4130-1.

Meyer-Bahlburg HFL, Dolezal C. The Female Sexual Function Index: a
methodological critique and suggestions for improvement. J Sex Marital Ther
2007;33:217-24. https://doi.org/10.1080/00926230701267852.

Mulamba KC. Maternal empowerment and child mortality in the Democratic Republic
of the Congo: evidence from multilevel modelling of DHS data. BMC Public Health
2026;26:2907. https://doi.org/10.1186/s12889-026-28504-9.

Nienkemper-Swanepoel J, le Roux NJ, Gardner-Lubbe S. On subset multiple
correspondence analysis for incomplete multivariate categorical data. Commun
Stat Simul Comput 2024;53:5229-41. https://doi.org/10.1080/03610918.2023.2173779.

Rai A, Khan MN, Thapa S. Trends and determinants of anaemia in women of Nepal:
a multilevel analysis. Matern Child Nutr 2020;16:e13044.
https://doi.org/10.1111/mcn.13044.

Spertus JA, Jones PG. Development and validation of a short version of the
Kansas City Cardiomyopathy Questionnaire. Circ Cardiovasc Qual Outcomes
2015;8:469-76. https://doi.org/10.1161/CIRCOUTCOMES.115.001958.

Tchakounté Tchuimi D, Kamga BF. The effect of women's bargaining power within
couples on contraceptive use in Cameroon. Gates Open Res 2020;4:20.
https://doi.org/10.12688/gatesopenres.13100.3.

Üstün TB, Kostanjsek N, Chatterji S, Rehm J, editors. Measuring health and
disability: manual for WHO Disability Assessment Schedule (WHODAS 2.0). Geneva:
World Health Organization; 2010.

Varghese JS, Maluccio JA, Cunningham SA, Ramirez-Zea M, Stein AD. Development
of a temporally harmonized asset index: evidence from across 50 years of follow
up of a birth cohort in Guatemala. BMC Med Res Methodol 2021;21:85.
https://doi.org/10.1186/s12874-021-01263-4.

Yule M, Davison J, Brotto L. The International Index of Erectile Function: a
methodological critique and suggestions for improvement. J Sex Marital Ther
2011;37:255-69. https://doi.org/10.1080/0092623X.2011.582431.

Every other source in this supplement is in the reference list of the main
text.

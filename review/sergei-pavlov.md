# Review of NHANES sample weights for BMI and fasting glucose

Reviewer: Sergei Pavlov. Verdicts refer to the original implementation. I ran all nine cases with the proposed fixes and checked selected standard errors with a delete-one-PSU jackknife.

## Verdicts

| Question | Cases | Verdict | Reason | Reference |
|---|---|---|---|---|
| 1. Weight choice | 1, 3 | Correct | The fasting subsample weight applies to the joint analysis of BMI and fasting glucose. | [CDC weighting], [GLU_J] |
| 2. Strata and PSUs | 2 | Correct | Taylor linearization compares PSU influence totals within strata. All 30 PSUs contribute to the design with 15 strata. | [CDC Taylor linearization] |
| 3. Rescaling | 4, 5, 7 | Correct | The common scaling factor cancels in both the correlation and its linearized SE. | [CDC weighting] |
| 4. Survey choice | 5–7 | Correct | Starting in 2003 simplifies weighting. NCHS permits pooling the later cycles, and the pooled result describes the two selected periods. | [CDC weighting], [CDC release] |
| 5. Filters | 8, 9 | Correct | The domain analysis preserves the survey design for SE estimation. | [CDC subgroup analysis] |
| 6. Inference | Weighted cases | Debatable | Both Wald approximations are reasonable. The original interval and test use different scales. | [CDC design df] |
| 7. Missing checks | 1, 5 | Debatable | The naive comparison changes the eligible sample. Pooled glucose measurements need a documented common scale for substantive interpretation. | [GLU_J], [GLU_L] |

## Findings and proposed fixes

**Q5 · Minor inference correction.** The domain analysis is correct. The Case 9 p-value uses the subgroup degrees of freedom (22 PSUs − 13 strata = 9). The same df is used for the confidence interval. r and SE are unchanged. [CDC subgroup df]

**Q6 · Minor.** Use the Fisher-z statistic (`z / se_z`) for the two-sided p-value so the test and 95% CI agree at the 5% threshold. In Case 1, p moves from 0.000073 to 0.000096 with the same r, SE and CI. With the subgroup df correction, Case 9 has p ≈ 0.00049.

**Q7 · Minor for the demonstration.** It would be useful to compare the weighted and naive correlations using the same participants with positive fasting subsample weights.

**Q7 · Minor for the demonstration.** In addition, it would be useful to compare the pooled correlation in Case 5 before and after bringing fasting glucose measurements from all eight cycles onto a common scale. In 2015, NHANES switched from the Roche Cobas C501 to the C311 ([GLU_I]). The published 2015–2016 results from samples measured on the C501 were already adjusted to the C311 scale ([GLU_L]).

## Rationale for the code changes

### Domain degrees of freedom (Q5, Case 9)

Taylor SE uses the original survey design. A PSU with no eligible subgroup observations contributes zero, so a stratum with one represented PSU still contributes to the variance estimate through the comparison with the other PSU. For inference, NCHS counts the PSUs and strata represented in the subgroup. Case 9 has 22 such PSUs in 13 strata, giving df = 9. The smaller df widens the interval while retaining the same correlation and SE. [CDC subgroup analysis] · [CDC subgroup df]

### Consistent scale for the interval and test (Q6)

The interval uses Fisher’s transformation, `z = arctanh(r)`, and the delta-method standard error `se_z = se / (1 − r²)`. The original test used `r / se`. Computing the two-sided p-value from `z / se_z` with the same df makes the test agree with the interval about whether zero is excluded at the 5% level. The change affects the reported p-value. Both calculations still use the survey-based SE and df.

### Degrees of freedom for the difference between cycles (Case 6)

Under the assumption of independent cycles, the estimated variance of the difference is the sum of the two squared SEs. The two cycles have 15 df each, but the earlier cycle contributes roughly two-thirds of this variance. The Welch–Satterthwaite approximation accounts for these unequal contributions and gives effective df ≈ 27.2 instead of 30. The estimated difference remains 0.034 with SE 0.0447. The resulting 95% CI is [−0.058, 0.126] and p ≈ 0.45. [NIST]

[CDC weighting]: https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx
[CDC Taylor linearization]: https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx#:~:text=Taylor%20Linearization%20Procedures
[CDC subgroup analysis]: https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx#:~:text=Analyzing%20Subgroups%20in%20NHANES
[CDC design df]: https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx#:~:text=Degrees%20of%20Freedom%20for%20a%20Complex%20Survey
[CDC subgroup df]: https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx#:~:text=Degrees%20of%20Freedom%20when%20Analyzing%20Subgroups%20in%20NHANES
[CDC release]: https://wwwn.cdc.gov/nchs/nhanes/ContinuousNhanes/releasenotes.aspx?Cycle=2021-2023
[GLU_J]: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/GLU_J.htm
[GLU_I]: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/GLU_I.htm
[GLU_L]: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/GLU_L.htm
[NIST]: https://www.itl.nist.gov/div898/handbook/mpc/section5/mpc571.htm

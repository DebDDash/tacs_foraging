# Variation in metrics pre and post stimulation during decision making while playing game

We first performed **epoching** on EEG recordings for each participant around 'Trigger#3' taking 1000 ms data.This corresponds to time taken to make stay or leave decision in the foraging game.

## Feature Extraction

For every participant and condition, we computed the following EEG metrics across all epochs:

| Metric | Description |
|:--------|:-------------|
| **Theta/Beta Ratio** | Indicator of attentional control and arousal balance |
| **Alpha/Beta Ratio** | Reflects relaxation vs alertness states |
| **PAC (Phase–Amplitude Coupling)** | Measures coupling between theta phase and gamma amplitude |
| **Average Alpha Power** | Reflects cortical idling or relaxed attentional state |
| **Average Beta Power** | Linked with motor engagement or alertness |
| **Average Gamma Power** | Associated with higher cognitive processing |

All metrics were calculated for both **pre** and **post** conditions and saved in a combined CSV file (`EEG_metrics_all_participants.csv`).

---

## Statistical Analysis

We then conducted **paired-sample t-tests** (and **Wilcoxon signed-rank tests** as nonparametric alternatives) to assess whether each metric significantly changed from pre to post stimulation across participants.  
The effect size (**Cohen’s d**) was also computed to evaluate the magnitude of change.

---

## Results Summary

| Metric | PreMean | PostMean | tStat | pValue | Cohen’s d | Wilcoxon p |
|:--------|---------:|----------:|-------:|--------:|-----------:|-------------:|
| **ThetaBeta** | 2.8316 | 2.5468 | 1.5955 | 0.1451 | -0.5046 | 0.3223 |
| **AlphaBeta** | 1.0964 | 1.4121 | -2.0785 | 0.0674 | 0.6573 | **0.0488** |
| **PAC** | 0.00022 | 0.00013 | 1.0226 | 0.3332 | -0.3234 | 0.6953 |
| **Alpha** | 21.148 | 24.400 | -2.2012 | 0.0552 | 0.6961 | 0.0645 |
| **Beta** | 19.357 | 17.674 | 0.9236 | 0.3798 | -0.2921 | 0.1934 |
| **Gamma** | 9.0388 | 6.6221 | 1.3809 | 0.2006 | -0.4367 | 0.1309 |

---

## Interpretation

- **Alpha/Beta Ratio:**  
  Showed a **significant increase**, indicating a possible **shift toward higher alpha dominance** after stimulation.  
  This suggests participants may have transitioned to a **more relaxed or internally focused state** post-stimulation.

- **Alpha Power:**  
  Showed a **trend-level increase**, consistent with enhanced alpha activity, often associated with **reduced cortical excitability** or **restful alertness**.

- **Other Metrics (Theta/Beta, Beta, Gamma, PAC):**  
  No statistically significant changes were detected.


---

## Files

- `EEG_metrics_all_participants.csv` — extracted EEG metrics per participant  
- `EEG_pre_post_stats.csv` — results of paired t-tests, p-values, and effect sizes  
- `analyze_EEG_stats.m` — MATLAB script performing the statistical tests  
- `band_power.m` and `compute_PAC.m` — helper functions for feature extraction  


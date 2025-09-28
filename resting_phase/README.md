# Analysis of data extracted during resting phase 
## analysis.m
Computes statistical metrics for pre vs. post stimulation:
- Cohen’s d – Positive = post > pre, Negative = post < pre.
- Effect size r – Sign shows direction (positive = increase, negative = decrease).
- Mean Difference – Absolute change in original units (e.g., µV²).
- Percent Change – Relative change compared to pre.
- p-value – < 0.05 indicates statistical significance.
- Z-value – Magnitude shows strength of evidence, sign shows direction.
Analyses performed for theta, alpha, beta, and gamma bands.
## features.m
Extracts EEG features including:
- PAC (Phase-Amplitude Coupling): Quantifies cross-frequency coupling.
- PLV (Phase Locking Value): Measures synchronization between signals.
## comodulogram.m
Generates comodulograms to visualize cross-frequency coupling.
## results
Has analysis and outputs of above codes



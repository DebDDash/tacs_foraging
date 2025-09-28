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

## Key Findings (Resting Phase)  

- **Significant Effects:**  
  - **theta–gamma PAC (F5–Fz):** Significant decrease (*p = 0.0443*).  

- **Slightly Significant / Trend-Level Effects:**  
  - **theta–theta PLV (F5–AFz):** Increase (*p = 0.116*).  
  - **theta–theta PAC/PLV (F5–FCz):** Weak trends (*p = 0.101–0.131*).  
  - **theta–theta PAC/PLV (F5–Fz):** Weak trends (*p = 0.13–0.14*).  

We see a **decrease in theta–gamma coupling** at the F5–Fz electrode pair, with additional weaker trends in **theta–theta connectivity**. 


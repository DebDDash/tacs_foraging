# Epoching

This pipeline processes EEG data in EDF format, extracts 1-second pre-trigger epochs around Trigger#3, and saves the data in a MATLAB/EEGLAB-friendly format with annotations for stay and leave decisions.

Outputs include 
- EEGLAB .set file for MATLAB analysis.
- Metadata CSV with epoch index, decision, reaction time, and reward.
- Visualization of ERPs

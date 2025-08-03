# Enhancing Empathic Accuracy: Penalized Functional Alignment Method to Correct Temporal Misalignment in Real-time Emotional Perception

Empathic accuracy (EA) is the ability to accurately understand another person's thoughts and feelings, which is crucial for social and psychological interactions. Traditionally, EA is assessed by comparing a perceiver’s moment-to-moment ratings of a target’s emotional state with the target’s own self-reported ratings at corresponding time points. However, misalignments between these two sequences are common due to the complexity of emotional interpretation and individual differences in behavioral responses. Conventional methods often ignore or oversimplify these misalignments, for instance by assuming a fixed time lag, which can introduce bias into EA estimates. 

To address this, we propose a novel alignment approach that captures a wide range of misalignment patterns. Our method leverages the square-root velocity framework to decompose emotional rating trajectories into amplitude and phase components. To ensure realistic alignment, we introduce a regularization constraint that limits temporal shifts to ranges consistent with human perceptual capabilities. We validate our method through simulations and real-world applications involving video and music datasets.

# Reproducing Simulation and Case Studies
1. Simulation studies \
- run `simulation_1_main_*.R`, `simulation_2_main.R`, and `simulation_3_main_*.R` for simulation studies.
- run `simulation_1_summary.R`, `simulation_2_summary.R`, and `simulation_3_summary.R` for summaries.
  
2. Case study 1 \
- run “alignment.R” to produce “Delvinetal_allcorrelations_08012025.Rdata”
- run “postAlignment-analyses.R” to produce all plots in the folder plots
- run “covariates-analysis.R” to produce Table 3 in the main ms
3. Case study 2 \
- run “alignment.R” to produce “Tabaketal_allMeasures_06012025.Rdata” (this will call the file “casestudy2-ancillary-func.R)
- run “Tabaketal_plots_threshold6s” to produce all plots in plots/Tabak_threshold6s.pdf
- run “Tabaketal_plots_threshold8s” to produce all plots in plots/Tabak_threshold8s.pdf
- run “Tabaketal_plots_threshold10s” to produce all plots in plots/Tabak_threshold10s.pdf

* For running `*_l1p.R`, use `fdasrvf` version >= 2.0.0. For the others, use the attached `fdasrvf_1.9.8.tar.gz`.

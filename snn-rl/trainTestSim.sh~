#!/bin/bash
# Script to train and test the simulation with running it from the command line

# train
RUN_TIME_VAR=55010
POS_REINF_VAR=1.15519
NEG_REINF_VAR=5.44256
RANDOMIZATION_VAR="0.0713-0.9546"
SHOW_PLOT_VAR=False

python gupta_paper_further_formulas.py evaluateClassifier=False runTime=$RUN_TIME_VAR posReinf=$POS_REINF_VAR negReinf=$NEG_REINF_VAR standardPrint=False verbosePrint=False randomization=$RANDOMIZATION_VAR showPlot=$SHOW_PLOT_VAR

# test
RUN_TIME_VAR=7900
SHOW_PLOT_VAR=True

python gupta_paper_further_formulas.py evaluateClassifier=True runTime=$RUN_TIME_VAR optResultsFile="optimizationResults.txt" posReinf=$POS_REINF_VAR negReinf=$NEG_REINF_VAR standardPrint=False verbosePrint=False randomization=$RANDOMIZATION_VAR showPlot=$SHOW_PLOT_VAR


RWR With eta parameter 0.01 
===

Drug Target Prediction by Random walk with restart 
There are 5 scripts as follows:
1. Transition.R - This script contains the codes to generate four trasition matrices.
2. rwr_predict.R - This core read the full M matrix which is the TT,TD, DT, DD matrices into a single matrixruns the prediction algorithm.  
3. log_score.R- This script generated the log score of the drug target and target drug matrix and averages it and generates the scores.It also generates raw score matrix.
4. dtp.R - Predicts the chembl dataset which is extracted from chembl verision 15. There are 544 drugs and 467 target proteins arranged into a matrix. 
5. Ranking.R - Takes the predicted network and finds the ranks of the targets and writes the csv file.Also calculated the rate of fraction of interactions.
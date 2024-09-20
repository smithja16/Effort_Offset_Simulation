This code evaluates various approaches to including sampling effort in GLMs, including offset terms. There is an accompanying maniscript in draft.

The code works in three stages: 
1) generates 'data scenarios', which are count data determined by three variables, with flexibility to apply collinearity among these variables as well as change the magnitude and shape of some effects
2) fits a varity of GLMs to these data scenarios
3) measures performance of each model for each data scenario, to gain insight into offset terms and variously structured effort covariates

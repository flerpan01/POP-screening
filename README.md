# Project summary

Environmental contaminants are present in human tissues, and many persistent organic pollutants (POPs) with endocrine disruptive properties are detectable in human serum.<br>
In this study the [OECD test guideline #456 steroidogenesis assay](https://www.oecd-ilibrary.org/environment/test-no-456-h295r-steroidogenesis-assay_9789264122642-en) was scaled up to screen 24 POPs in 96-well microplates for their effects on testosterone and estradiol synthesis, as well as cell viability, using H295R cells. The chemical test library included six polyfluoroalkyl substances (PFAS), five organochlorine pesticides, ten polychlorinated biphenyls (PCBs) and three polybrominated diphenyl ethers (PBDEs), and all compounds were tested at eight concentrations (1 nM to 10 µM). The data was analyzed using a linear mixed-effects model, and Dunnett’s multiple comparison was applied as a post-hoc test.

>The code was used in this paper: [Screening persistent organic pollutants for effects on testosterone and estrogen synthesis at human-relevant concentrations using human H295R cells in 96-well plates](http://dx.doi.org/10.1007/s10565-024-09902-4)

## Project structure

The dataset is composed of human _H295R_ cells exposed to 24 POPs compounds with measures collected for activity of estradiol, testosteron and MTT assay. The molar concentration of the compounds range from 0 (control) to 1e-5 mol/L. Experiments were ran in 3 batches. The dataset is curated in an excel file, `doc/steroidogenesis_dataset.xlsx`. `code/analysis.R` contains the code for statistical analysis (linear mixed-effects model, ANOVA and post-hoc) and will output (1) plots of all compounds (`img/plots.pdf`) and (2) the statistical output (`data/fit_linear_mixed_model_output.xlsx`). All analysis were run in [R](https://www.r-project.org/)

```
project/
├─- code
│   └-- analysis.R
├-- data
│   └-- fit_linear_mixed_model_output.xlsx
├-- doc
│   └-- steroidogenesis_dataset.xlsx
├-- img
│   └-- plots.pdf
└-- README.md
```

>This project was made by the [Karlsson Laboratory](https://karlssonlab.se/) at Stockholm University, Sweden
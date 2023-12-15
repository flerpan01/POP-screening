# Project summary

Environmental contaminants are present in human tissues, and many persistent organic pollutants (POPs) with endocrine disruptive properties are detectable in human serum.<br>
In this study the [OECD test guideline #456 steroidogenesis assay](https://www.oecd-ilibrary.org/environment/test-no-456-h295r-steroidogenesis-assay_9789264122642-en) was scaled up to screen 24 POPs in 96-well microplates for their effects on testosterone and estradiol synthesis, as well as cell viability, using H295R cells. The chemical test library included six polyfluoroalkyl substances (PFAS), five organochlorine pesticides, ten polychlorinated biphenyls (PCBs) and three polybrominated diphenyl ethers (PBDEs), and all compounds were tested at eight concentrations (1 nM to 10 µM). The data was analyzed using a linear mixed-effects model, and Dunnett’s multiple comparison was applied as a post-hoc test.

>The code was used in this paper: [Large scale in vitro screen for steroidogenic disruption by persistent organic pollutants found in blood reveal effects on testosterone, estradiol and viability at low concentration in H295R cells](https://doi.org/)

## Project structure

The dataset is curated in an excel file, `misc/steroidogenesis\ R\ 2.0.xlsx`. `code/analysis.R` contains the code which runs the statistical analysis (linear mixed-effects model), produces an excel, `data/fit_linear_mixed_model_output.xlsx`, with the model outputs and also produces plots for each compound, `img/plots_w_sign_indicator.pdf`

```sh
proj/
|-- code
|   `-- analysis.R
|-- data
|   `-- fit_linear_mixed_model_output.xlsx
|-- img
|   `-- plots_w_sign_indicator.pdf
|-- misc
|   `-- steroidogenesis\ R\ 2.0.xlsx
`-- readme.md
```
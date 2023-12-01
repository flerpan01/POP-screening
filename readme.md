# Project summary

Målet är att screena mixer av kemikalier i blod efter hormonstörande effekter. Vi använder oss av in vitro assays för att mäta effekter på steroidogenes av steroidsexhormonen estradiol och testosteron som är rekommenderade  av OECD. Metoden har skalats upp och validerats.
Kemikarlier vi valt att fokusera på är sk. POPs, "persistent organic pollutants" som uppmätts i blod i en kohort från Västerbotten, och av dessa valde vi de 24 som återfanns i högst koncentration i ett dataset med nära 800 individer.

För att undersöka mix-effekter och jämföra mixer från olika individer med varandra har vi först undersökt alla ämnen var för sig, och sedan återskapat mixer såsom uppmätts i olika individer och applicerat dessa på celler för att kvantifiera effekt på steroidogenes.

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
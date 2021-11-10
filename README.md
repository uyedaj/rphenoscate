# Phylogenetic Comparative Analysis of Integrated Anatomical Traits: A Short Course

This R package encapsulates the material for the short course on _Phylogenetic Comparative Analysis of Integrated Anatomical Traits_.

The purpose of doing so is to leverage the R package installation mechanism for ensuring that all R dependencies needed for the course are installed and available, simply as a result of installing this package.

The material itself is in the form of vignettes written in RMarkdown. Users are encouraged to modify the code chunks and re-execute the vignettes as a means of exploring the material.

### Installation

You can download the course material to your computer directly from [Github](http://github.com/phenoscape/scate-shortcourse) as a Zip archive. Click the _Download Zip_ link that appears under the _Clone or Download_, move the Zip archive to a directory of your choice, then unpack it. 
![](https://i.imgur.com/qL9NZ7L.png)

You can then open the folder containing the material as a project in RStudio, and render ("knit") and modify the vignette(s) (see `vignette/` subdirectory).

The RPhenoscape vignette depends only on the [RPhenoscape package](http://rphenoscape.phenoscape.org), which can be installed directly from Github as per its documentation. (RPhenoscape has itself a number of recursive dependencies; installing the package will install those, too.)

The PARAMO vignette has a number of dependencies. These can be installed automatically by installing the workshop package, either from the command line (use `R CMD build`, followed by `R CMD INSTALL`), or from within RStudio using the `install_github()` function in the `remotes` package (which can be installed from CRAN using `install.packages()`):

```r
remotes::install_github("phenoscape/scate-shortcourse")
```

This command will likely trigger a (possibly long) list of packages that are required to be installed and upgraded, or suggested for upgrade, with a prompt for choosing how to select which packages to upgrade or install. Unless you happen to have installed an earlier version of RPhenoscape that fails to meet the required version, we recommend to use the option for installing and upgrading _from CRAN only_. (Typically, this will be option 2.)

Note that the above command will use the default options for build arguments, which will leave out building and installing the vignettes. Since the best use of the vignettes is to go through (and potentially modify) them within RStudio, this shouldn't be an issue.

#### Installation on Unix / Linux

If your operating system is some flavor of Linux, many of the package dependencies that include compiled code will end up being installed from source. Some of these have dependencies external to R, and thus installing them from source will fail if these external dependencies are not pre-installed.

Currently, the following packages are required on top of a Ubuntu 18.04 / R 5.3.3 environment:
- libmagick++-dev
- libgmp-dev
- libmpfr-dev
- pandoc

Other base systems may require additional (or fewer) external packages.

### When and Where
- [SSB 2020] Workshop in Gainesville, FL
- January 6, 2020 (Monday)
- Venue: Hotel Indigo

### Schedule
* 8:30am	Participants arrive and get settled
* 9:00	Welcome and intros, goals
* 9:15	Background concepts: anatomy ontologies, phenotype annotations and character dependencies
    - In parallel: tech check troubleshooting
* 9:30	Using [Phenoscape Knowledgebase]
* 10:00 _Break_
* 10:30	Using [RPhenoscape]
* 10:45	Constructing stochastic character maps on a phylogeny
    - Using the PARAMO pipeline with RevBayes
    - Obtaining a presence-absence dependency matrix
* 11:30	Measuring semantic similarity
    - In RPhenoscape
    - Similarity between phenotype profiles in the KB
* 11:45	How to manage your own data
    - Demo of Phenex and Phenoscape annotation guidelines
    
### Instructors
- [Jim Balhoff](https://orcid.org/0000-0002-8688-6599) Renaissance Computing Institute
- [Wasila Dahdul](https://scholar.google.com/citations?user=qHfrfGwAAAAJ&hl=en), U South Dakota
- [Hilmar Lapp](https://scholars.duke.edu/person/Hilmar.Lapp), Duke U
- [Paula Mabee](https://www.usd.edu/faculty-and-staff/Paula-Mabee), U South Dakota
- [Josef Uyeda](https://www.uyedalab.com/), Virginia Tech
- [Todd Vision](https://orcid.org/0000-0002-6133-2581), UNC Chapel Hill

### Recommended papers and resources
- Dececchi et al. Toward synthesizing our knowledge of morphology: using ontologies and machine reasoning to extract presence/absence evolutionary phenotypes across studies. Systematic Biology 64, 936. 2015. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4604830/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4604830/)
- Jackson et al. (2018) Automated integration of trees and traits: a case study using paired fin loss across teleost fishes. Systematic Biology, 67(4):559–575, [https://doi.org/10.1093/sysbio/syx098](https://doi.org/10.1093/sysbio/syx098)
- Tarasov et al. (2019). PARAMO: A pipeline for reconstructing ancestral anatomies using ontologies and stochastic mapping. Insect Systematics and Diversity, Volume 3, Issue 6, [https://doi.org/10.1093/isd/ixz009](https://doi.org/10.1093/isd/ixz009)
- Guide to Character Annotation [https://wiki.phenoscape.org/wiki/Guide_to_Character_Annotation](https://wiki.phenoscape.org/wiki/Guide_to_Character_Annotation)
- Dahdul et al. (2018) Phenoscape guide to character annotation. figshare. Fileset. [https://doi.org/10.6084/m9.figshare.1210738.v2](https://doi.org/10.6084/m9.figshare.1210738.v2)


### Funding
[SCATE] is funded by the US National Science Foundation (NSF) as collaborative awards [1661456] (Duke University), [1661529] (Virginia Tech), [1661516] (University of South Dakota), and [1661356] (UNC Chapel Hill and RENCI) within the Advances in Biological Informatics (ABI) program.

## Terms of reuse

To the extent possible under law, the [SCATE] Project Team has waived all copyright and related or neighboring rights to this work. See the [CC Zero] public domain waiver for details.

If you reuse or make derivatives of this work, we request that you cite or otherwise give proper credit to this work.

[Phenoscape Knowledgebase]: http://kb.phenoscape.org/#/home
[SCATE]: http://scate.phenoscape.org
[CC Zero]: https://creativecommons.org/publicdomain/zero/1.0/
[RPhenoscape]: http://rphenoscape.phenoscape.org/
[Evolution Meetings]: https://www.evolutionmeetings.org/evolution-2019---providence.html
[1661456]: https://nsf.gov/awardsearch/showAward?AWD_ID=1661456
[1661529]: https://nsf.gov/awardsearch/showAward?AWD_ID=1661529
[1661356]: https://nsf.gov/awardsearch/showAward?AWD_ID=1661356
[1661516]: https://nsf.gov/awardsearch/showAward?AWD_ID=1661516
[SSB 2020]: https://systbiol.github.io/ssb2020/index.html

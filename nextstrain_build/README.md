### Nextstrain builds for ORF8 KO
In this manuscript, we built 5 different SARS-CoV-2 nextstrain trees:
- [WA_20k](https://nextstrain.org/groups/blab/ncov/WA/20k), which is enriched for WA sequences and evenly sampled across time.
- [Alpha](https://nextstrain.org/groups/blab/ncov/WA/Alpha), an Alpha specific phylogeny enriched for potential ORF8 knockouts in WA.
- [Delta](https://nextstrain.org/groups/blab/ncov/WA/Delta), a Delta specific phylogeny enriched for potential ORF8 knockouts in WA.
- [Omicron](https://nextstrain.org/groups/blab/ncov/WA/Omicron), an Omicron specific phylogeny enriched for potential ORF8 knockouts in WA.
- [WA_other](https://nextstrain.org/groups/blab/ncov/WA/other), non-Delta, non-Omicron, non-Alpha phylogeny enriched for potential ORF8 knockouts in WA.

`orf8ko_profiles` contains the code to build these trees using the [Nextstrain ncov workflow](https://github.com/nextstrain/ncov) v13. 
To build the trees, [install the ncov workflow](https://docs.nextstrain.org/projects/ncov/en/latest/tutorial/setup.html) and copy the `orf8ko_profiles` folder into the `ncov` directory.
Run the workflow from inside the `ncov` directory using the below command:
```
nextstrain build . --configfile orf8ko_profiles/orf8ko_builds.yaml
```

`results` contains the knockout clusters for each gene inferred for each tree. 

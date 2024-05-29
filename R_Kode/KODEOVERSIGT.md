# Readthrough_analysis R-kode oversigt
Under projektforløbet har vi bidraget til at generalisere og automatisere det readthrough_analysis script, som Søren Lykke Andersen (fremadrettet SLA) tidligere har udarbejdet.
Vi har brudt hele den samlede readthrough_analysis kode op i mindre dele (funktioner) for at gøre scriptet mere overskueligt. Dette er ikke af praktiske årsager, og det er ikke hensigten, at scriptet skal køres på denne måde.

### Bidrag
Nedenstående tabel indeholder en kort oversigt over projektets kode-bidrag. Under tabellen fremgår en mere dybdegående oversigt, hvor de specifikke kodestykkers funktionalitet er forklaret.  

| Funktion      | Fælles indsats| SLA           | 
| ------------- | ------------- | ------------- |
| check_upstream_signal  | X  |   | 
| hmm_and_ds_fitting  | X  |   | 
| load_GOI  | X  |  | 
| log_2_transform  | X  |   | 
| normalize_to_gene_body_signal  | X  |   | 
| subtract_ctrl_signal  | X  |   | 
| wrapper  | X  |   | 
| remove_batch_effects  |   | X  | 
| visualisation_prep  |   | X  | 
| visualisation  |   | X  | 
| SLA_adj_function  |   | X  | 
| SLA_binning_function  |   | X  | 
| SLA_sample_subtract_hmm  |   | X  | 
 <br>


* R-filer med prefixet 'SLA' er udelukkende produceret af SLA. De indgår derfor også kun i dette repository, da de bliver brugt i readthrough_analysis scriptet.
* Følgende 'funktioner' er udarbejdet af os, i samarbejde med SLA:
  - check_upstream_signal
    * I dette kode-stykke udregnes forskellen mellem signalet fra kroppen (elongering) og opstrøms-signalet.
  - hmm_and_ds_fitting
    * Vi forsøger først at fitte en HMM på vores kontrol fratrukket prøven, for at påvise termineringsdefekten. Hvis dette lykkes,                  fortsætter vi med, at fitte en sigmoidal-kurve (vi foretrækker en double-sigmoidal kurve), for at beskrive defekten. 
  - load_GOI
    * Her benytter vi et gen-navn (GOI, i.e. "Gene of Interest") sammen med vores annoteringsfil til at bestemme den region (range) vi er           interesserede i, for det aktuelle gen, med forbehold for strand (+ eller - DNA-streng). Herefter benytter vi denne range til at lave en       dataframe med coverage
      for vores kontrol & sample, på de positioner vi har udledt.
  - log_2_transform
    * Vi laver en log2 transformation af signalet for at gøre det mere symmetrisk, stabilt og lettere fortolkeligt. 
  - normalize_to_gene_body_signal
    * Her benytter vi en algoritme til at finde den 'bedste' TSS og TES ud af alle de transkripter der er beskrevet i annoteringen for det          aktuelle GOI, og derefter normalisere vi signalet ved at trække kroppens median fra.
  - subtract_ctrl_signal
    * Vi trækker kontrolsignalet fra prøvesignalet, hvilket medfører, at det resulterende signal ligger omkring nul, bortset fra                    ved defekten, hvor signalet vil afvige fra nul.
  - wrapper
    * Vores wrapper er designet til at automatisere kørslen af readthrough_analysis på genomeDK's cluster med muligheden for at behandle            flere GOI's i træk, samtidig med at man kan specificere den ønskede beregningskapacitet ved at angive antallet af kerner og                   mængden af RAM.
  
* De følgende funktioner er kun blevet tilpasset af SLA for at forblive kompatible med ændringerne i den tidligere readthrough_analysis-kode.
  - remove_batch_effects
  - visualisation_prep
  - visualisation
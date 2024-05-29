# Readthrough_analysis R-kode oversigt
Under projektforløbet har vi bidraget til at generalisere og automatisere det readthrough_analysis script, som Søren Lykke Andersen (SLA) tidligere har udarbejdet.
Vi har brudt hele den samlede readthrough_analysis kode op i mindre dele (funktioner) for at gøre scriptet mere overskueligt. Dette er ikke af praktiske årsager, og det er ikke hensigten, at scriptet skal køres på denne måde.

### Bidrag
I denne tabel findes en kortfattet oversigt over kode bidrag:

| Funktion      | Fælles indsats| SLA           | 
| ------------- | ------------- | ------------- |
| Content Cell  | Content Cell  | Content Cell  | 
| Content Cell  | Content Cell  | Content Cell  | 



* R-filer med prefixet 'SLA' er udelukkende produceret af SLA. De indgår derfor også kun i dette repository, da de bliver brugt i readthrough_analysis scriptet.
* Følgende 'funktioner' er udarbejdet af os, i samarbejde med SLA:
  - check_upstream_signal
    * I dette kode-stykke udregnes forskellen mellem signalet fra kroppen (elongering) og opstrøms-signalet.
  - hmm_and_ds_fitting
    * Vi forsøger først at fitte en HMM på vores kontrol fratrukket prøven, for at påvise termineringsdefekten. Hvis dette lykkes,                  fortsætter vi med, at fitte en sigmoidal-kurve (vi foretrækker en double-sigmoidal kurve), for at beskrive defekten. 
  - load_GOI
    * Her benytter vi et gen-navn (GOI, i.e. "Gene of Interest") sammen med vores annoteringsfil til at bestemme den region (range) vi er           interesserede i, for det aktuelle gen, med forbehold for strand (+ eller - DNA-streng). Herefter benytter vi denne range til at lave en dataframe med coverage
      for vores kontrol & sample, på de positioner vi har udledt.
  - log_2_transform
    * Vi laver en log2 transformation af signalet for at gøre det mere symmetrisk, stabilt og lettere fortolkeligt. 
  - normalize_to_gene_body_signal
    * Her benytter vi en algoritme til at finde den 'bedste' TSS og TES ud af alle de transkripter der er beskrevet i annoteringen for det          aktuelle GOI, og derefter normalisere vi signalet ved at trække kroppens median fra.
  - subtract_ctrl_signal
    * Vi trækker kontrolsignalet fra prøvesignalet, hvilket medfører, at det resulterende signal ligger omkring nul, bortset fra                    ved defekten, hvor signalet vil afvige fra nul.
  - wrapper
    * Vores wrapper er designet til at automatisere kørslen af readthrough_analysis på genomeDK's cluster med muligheden for at behandle            flere GOI's i træk, samtidig med at man kan specificere den ønskede beregningskapacitet ved at angive antallet af kerner og                   mængden af RAM.
  
* De følgende funktioner er kun blevet let tilpasset af SLA (udelukkende tilpasset med henblik på at fungere sammen med ændret kode fra det tidligere readthrough_analysis script).
  - remove_batch_effects
  - visualisation_prep
  - visualisation

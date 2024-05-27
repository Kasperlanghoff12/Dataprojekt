# Readthrough_analysis R-kode oversigt

Under projektforløbet har vi bidraget til at generalisere og automatisere det readthrough_analysis script, som SLA tidligere har udarbejdet.
Vi har brudt hele den samlede readthrough_analysis kode op i mindre dele (funktioner) for at gøre scriptet mere overskueligt. Dette er ikke af praktiske årsager, og det er ikke hensigten, at scriptet skal køres på denne måde.

### Bidrag

* R-filer med prefixet 'SLA' er udelukkende produceret af SLA.
* Følgende 'funktioner' er udarbejdet i samarbejde med SLA:
  - check_upstream_signal
  - hmmm_and_ds_fitting
  - load_GOI
  - log_2_transform
  - normalize_to_gene_body_signal
  - subtract_ctrl_signal
  - wrapper
  
* De følgende funktioner er kun blevet let tilpasset af SLA (udelukkende tilpasset med henblik på at fungere sammen med ændret kode fra det tidligere readthrough_analysis script).
  - remove_batch_effects
  - visualisation_prep
  - visualisation

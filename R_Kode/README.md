# Readthrough_analysis R-kode oversigt

Under projektforløbet har vi bidraget til at generalisere og automatisere det readthrough_analysis script, som SLA tidligere har udarbejdet.
Vi har brudt hele den samlede readthrough_analysis kode op i mindre dele (funktioner) for at gøre scriptet mere overskueligt. Dette er ikke af praktiske årsager, og det er ikke hensigten, at scriptet skal køres på denne måde.

### Bidrag

R-filer med prefixet 'SLA' er udelukkende produceret af SLA.
Følgende 'funktioner' er udarbejdet i samarbejde med SLA:
- check_upstream_signal
- hmmm_and_ds_fitting
- load_GOI
- log_2_transform
- normalize_to_gene_body_signal
- subtract_ctrl_signal
  
Følgende funktioner er kun rettet minimalt (rettet med henblik på at køre sammen med andet modificeret kode) fra det forrige readthrough_analysis script:
- remove_batch_effects
- visualisation_prep
- visualisation

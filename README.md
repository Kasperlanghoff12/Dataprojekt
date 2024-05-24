# Model for defekter i transskriptionel terminering

## Intro

Molekylærbiologiens centrale dogme beskriver, hvordan det genetiske informationsflow løber fra DNA til RNA til protein. Der findes over 10.000 gener spredt ud i afgrænsede regioner i kromosomerne, som består af DNA. Hvert enkelt gen er forskriften for RNA, som produceres via transskriptionsprocessen. RNA bliver derefter aflæst via translation, hvilket fører til dannelsen af proteiner, som udfører det væld af processer, der foregår i en celle. <br> <br>
<p>
    <img width="450" alt="central_d" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/69349634-7729-42d3-898b-f45b653eb80e">
    <br>
    <em>Figur 1</em>
</p>
Hvert enkelt trin i det centrale dogme er fundamentalt for alt cellulært liv. Transskriptionen – dannelsen af RNA ved aflæsning af DNA’et fra et givet gen – er i sig selv en kompliceret proces, der består af en lang række koordinerede trin, der udføres af et stort antal proteiner. Transskriptionsprocessen kan deles op i tre faser – initiering, elongering og terminering. I forbindelse med initieringen rekrutteres en polymerase til genet vha. en promoter og påbegynder elongeringen, hvorved genet aflæses i en bestemt retning (se figur 2). Når polymerasen møder en terminator, igangsættes termineringen. Mens initieringen typisk sker fra et defineret område ved promoteren, er termineringen en mere stokastisk proces, hvor polymerasen bremses ned efter mødet med terminatoren og stopper mere eller mindre tilfældigt i nedstrømsregionen. På trods af denne stokasticitet er det dog klart, at processen skal være under stram kontrol, da polymerasen skal være fuldstændigt standset, inden det næste gen begynder. <br> <br>
<p>
    <img width="450" alt="transkription" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/343be33a-423c-4277-bba3-be89cdda21c7">
    <br>
    <em>Figur 2</em>
</p>
Således er transskriptions-terminering også en kompliceret proces, der involverer over 50 forskellige proteiner. Defekter i termineringen, hvor polymerasen ikke bremses ordentligt, fører til såkaldt “read-through”, hvor transskriptionen fortsætter længere end normalt i den nedstrøms region. Sådanne defekter kan opstå hos patienter, hvor et eller flere af de involverede proteiner er dysfunktionelle. Det kan også opstå i forbindelse med cellulært stress, som ofte opstår midlertidigt i normale individer, men mere ukontrollerbart i forbindelse med sygdomme som for eksempel kræft. Et åbent spørgsmål i feltet er, hvorvidt read-through blot er en konsekvens af stress, eller om det tjener en funktion og derfor er et forsøg fra cellen på at modvirke stress-tilstanden.
For at kunne forstå termineringen og dens involverede trin er det relevant at kunne udføre en detaljeret beskrivelse af read-through under forskellige betingelser. Vi har i dette projekt været med til at designe en model og algoritme, der kan lave en kvantitativ beskrivelse af read-through på individuelle gener baseret på transskriptionsdata fra celler udsat for forskellige betingelser.

## Data

### Beskrivelse

Hele processen for vores data er vist i figur 3. Efter at DNA- eller RNA-prøver er blevet sekventeret i små dele kaldet "reads", bliver de mange millioner sekvenslæsninger opbevaret i en FASTQ-fil. FASTQ-filen er et tekstbaseret filformat, der indeholder reads som tekststrenge, eksempelvis sekvensen af nukleotider: "GATTTGGGGTTC....". Derudover inkluderer hver FASTQ-fil også information om kvaliteten eller pålideligheden af læsningen for hver enkelt base i sekvensen.
<p>
    <img width="600" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/10cd3282-47dd-49b7-b935-0934898d9121">
    <br>
    <em>Figur 3</em>
</p>

Dette bliver konverteret til en BigWig-fil, som indeholder en værdi for, hvor meget hver position af genets sekvens er afdækket af reads fra FASTQ-filen, dvs. hvor meget data der er læst på hver position. Dette ses i “score”-kolonnen i figur 4, som er selve vores responsvariabel. Vi har også en annoteringsfil, der giver information om, hvor i BigWig-filen, de specifikke gener og transskriptioner er. Vi bruger så BigWig-filen og annoteringsfilen i sammenhold, så vi kan modellere afdækningen af specifikke gener og transskriptioner. Der er tale om en tidsserie, da krav om uafhængighed ikke er opfyldt, hvilket kan ses i autocorrelationsplottet i figur 6.
<p>
    <img width="800" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/5aaccbe1-0c13-4771-90ef-996d895be780">
    <br>
    <em>Figur 4</em>
</p>


I r-koden nedenfor importerer vi den specificerede region i BigWig-filen ud fra hvornår i annotering filen, et bestemt gen optræder.

```{r}
for (fname in ctrl) {
    ctrl_fname = paste0(ctrl_dir, fname, "_", strand_options[strand_sign], ".bw")
    ctrl_bw = import(ctrl_fname, 
                     which = GenomicRanges::GRanges(seqnames = chrom_no_loop, 
                                                    ranges = IRanges::IRanges(start = start_coord_ext, end = end_coord_ext)), 
                     as = "NumericList")[[1]]
    GOI_data[[fname]] = ctrl_bw
  }
```
Et eksempel på vores data ses i figur 5, hvor x-aksen er position på genet, og y- aksen er mængden af data normaliseret. Det ses så hvordan vores kontrol-data, den røde kurve, har minimal udsving efter terminering (nedstrømsregionen), mens test-data, den blå kurve, udsvinger meget efter termineringen. Dermed må der have været en defekt i den blå kurves terminering.

<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/1fe6dbb5-b324-4bba-a9ee-3ef65c098a99">
    <br>
    <em>Figur 5</em>
</p>

Der er muligvis biases i data, da de maskiner, der læser sekvenserne, kan producere støj. Derudover er der altid mulighed for menneskelige fejl. Vi vil dog antage, at biases er minimale, da sekvenserne er blevet læst på et laboratorium i et kontrolleret miljø.
<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/89b5271e-b39c-4069-a133-87aa2563c0a9">
    <br>
    <em>Figur 6</em>
</p>

### Preprocessering

For at generalisere generne så vi kan fitte modeller på dem, har vi brug for at normalisere vores data. Her starter vi med at lægge 1 til alle observationer, så når coverage er 0, forbliver det sådan efter log-transformationen. Siden DNA-strenge både kan læses forlæns og baglæns, var vi også nødt til at tage i betragtning hvilken retning vores datapunkter var.

```{r}
if (strand_sign == "-"){
    GOI_data = GOI_data[order(GOI_data$index, decreasing=TRUE), ]
  }
log2_GOI_data = log2(GOI_data + 1)
```

Udover dette, skulle vi også finde kroppen af generne, så vi kunne undersøge antagelsen om, at kroppen ville være ens for gener med- og uden defekt i termineringen. Udfordringen ved dette var, at vi i annoteringsfilen har flere annotationer af transkripter, hvor der var forskellige koordinater for begyndelsen (TSS) og slutningen (TES) af kroppen. Disse var vi nødt til at sammenligne, så vi var sikre på at finde de rigtige. 

Når vi kender kroppen, kan vi normalisere dem ved at trække medianen af kroppen fra for de enkelte gener. På denne måde sikrer vi en mere robust sammenligning af kontrol- og sample-gener.

```{r}
body.norm.factors = apply(log2_GOI_data[TSS_row:TES_row,], 2, median)
log2_GOI_data.bodynorm = t(t(log2_GOI_data) - body.norm.factors)
```
Til sidst fjernede vi systematisk bias fra kontrol-signalet. Dette gjorde vi ved at trække den normaliserede krops middelværdi for kontrol genet, fra den for sample genet.

```{r}
log2_GOI_data.bodynorm.means[, "ctrl"] = rowMeans(log2_GOI_data.bodynorm[, ctrl])
log2_GOI_data.bodynorm.means[, "sample"] = rowMeans(log2_GOI_data.bodynorm[, sample])

log2_GOI_data.bodynorm.means.ctrlsubtract = data.frame("sample" = log2_GOI_data.bodynorm.means[, "sample"] - log2_GOI_data.bodynorm.means[, "ctrl"])
```


## Algoritme

## Modellering

## Perspektivering

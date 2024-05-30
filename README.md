# Model for defekter i transskriptionel terminering

## Intro

Molekylærbiologiens centrale dogme beskriver, hvordan det genetiske informationsflow løber fra DNA til RNA til protein. Der findes over 10.000 gener spredt ud i afgrænsede regioner i kromosomerne, som består af DNA. Hvert enkelt gen er forskriften for RNA, som produceres via transskriptionsprocessen. RNA bliver derefter aflæst via translation, hvilket fører til dannelsen af proteiner, som udfører det væld af processer, der foregår i en celle. <br> <br>
<p>
    <img width="450" alt="central_d" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/69349634-7729-42d3-898b-f45b653eb80e">
    <br>
    <em>Figur 1 - Det centrale dogme (Transskriptionen foregår i processen fra DNA til
RNA).</em>
</p>
Hvert enkelt trin i det centrale dogme er fundamentalt for alt cellulært liv. Transskriptionen – dannelsen af RNA ved aflæsning af DNA’et fra et givet gen – er i sig selv en kompliceret proces, der består af en lang række koordinerede trin, der udføres af et stort antal proteiner. Transskriptionsprocessen kan deles op i tre faser – initiering, elongering og terminering. I forbindelse med initieringen rekrutteres en polymerase til genet vha. en promoter og påbegynder elongeringen, hvorved genet aflæses i en bestemt retning (se figur 2). Når polymerasen møder en terminator, igangsættes termineringen. Mens initieringen typisk sker fra et defineret område ved promoteren, er termineringen en mere stokastisk proces, hvor polymerasen bremses ned efter mødet med terminatoren og stopper mere eller mindre tilfældigt i nedstrømsregionen. På trods af denne stokasticitet er det dog klart, at processen skal være under stram kontrol, da polymerasen skal være fuldstændigt standset, inden det næste gen begynder. <br> <br>

<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/d2554b07-8fc2-4e61-9a99-672af91fad6b">
    <br>
    <em>Figur 2 - Transskription - Promotoren og terminator er begge DNA-sekvenser, der
markerer henholdsvis starten og afslutningen af transskriptionen.</em>
</p>
Således er transskriptions-terminering også en kompliceret proces, der involverer over 50 forskellige proteiner. Defekter i termineringen, hvor polymerasen ikke bremses ordentligt, fører til såkaldt “read-through”, hvor transskriptionen fortsætter længere end normalt i nedstrømsregionen. Sådanne defekter kan opstå hos patienter, hvor et eller flere af de involverede proteiner er dysfunktionelle. Det kan også opstå i forbindelse med cellulært stress, som ofte opstår midlertidigt i normale individer, men mere ukontrollerbart i forbindelse med sygdomme som for eksempel kræft. Et åbent spørgsmål i feltet er, hvorvidt read-through blot er en konsekvens af stress, eller om det tjener en funktion og derfor er et forsøg fra cellen på at modvirke stress-tilstanden.
For at kunne forstå termineringen og dens involverede trin er det relevant at kunne udføre en detaljeret beskrivelse af read-through under forskellige betingelser. Vi har i dette projekt været med til at designe en model og algoritme, der kan lave en kvantitativ beskrivelse af read-through på individuelle gener baseret på transskriptionsdata fra celler udsat for forskellige betingelser.
<br> <br>
Vores readthrough-analyse proces ses i figur 3. Efter sekventeringen starter vores proces med at klargøre data, som beskrevet i Data-afsnittet. Derefter pre-processerer vi vores data med log-transformering og normalisering for at vi kan sammenligne test- og kontrol-gener. Så identificerer vi transskriptionens start- og slut-sites, som definerer de data-områder som vi modellerer på. Vi fitter så for hvert gen en Hidden Markov Model, og bruger Viterbi decoding til at identificere om der er en defekt, og i hvilket område som defekten sker. Hvis der er det, fitter vi en (double) sigmoid til at beskrive defekten.
<br> <br>
<p>
    <img width="650" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/7af7be05-1366-4a6c-b3cf-beed91c7f237">
    <br>
    <em>Figur 3 - Projekt pipeline</em>
</p>

## Data

### Beskrivelse

Hele processen for vores data er vist i figur 4. Efter at DNA- eller RNA-prøver er blevet sekventeret i små dele kaldet "reads", bliver de mange millioner sekvenslæsninger opbevaret i en FASTQ-fil. FASTQ-filen er et tekstbaseret filformat, der indeholder reads som tekststrenge, eksempelvis sekvensen af nukleotider: "GATTTGGGGTTC....". Derudover inkluderer hver FASTQ-fil også information om kvaliteten eller pålideligheden af læsningen for hver enkelt base i sekvensen. 
<p>
    <img width="600" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/10cd3282-47dd-49b7-b935-0934898d9121">
    <br>
    <em>Figur 4 - Data pipeline </em>
</p>

Dette bliver konverteret til en BigWig-fil, som indeholder en værdi for, hvor meget hver position af genets sekvens er afdækket af reads fra FASTQ-filen, dvs. hvor meget data der er læst på hver position. Dette ses i “score”-kolonnen i figur 5, som er selve vores responsvariabel. Vi har også en annoteringsfil, der giver information om, hvor i BigWig-filen, de specifikke gener og transskriptioner er. BigWig-filen indeholder 10.386 forskellige gener, hvilke spænder over samlet 654.074.921 datapunkter, altså knap 63.000 signalværdier for hvert gen i gennemsnit. Vi bruger så BigWig-filen og annoteringsfilen i sammenhold, så vi kan modellere afdækningen af specifikke gener og transskriptioner. Der er tale om sekventiel data, da krav om uafhængighed ikke er opfyldt, hvilket kan ses i autocorrelationsplottet i figur 7. 
<p>
    <img width="800" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/5aaccbe1-0c13-4771-90ef-996d895be780">
    <br>
    <em>Figur 5 - Eksempel på de filtyper, vi arbejder med</em>
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
Et eksempel på vores data ses i figur 6, hvor x-aksen er position på genet, og y- aksen er mængden af data normaliseret. Det ses så hvordan vores kontrol-data, den røde kurve, har minimale udsving efter terminering (nedstrømsregionen), mens test-data, den blå kurve, udsvinger meget efter termineringen. Dermed må der have været en defekt i den blå kurves terminering.

<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/1fe6dbb5-b324-4bba-a9ee-3ef65c098a99">
    <br>
    <em>Figur 6 - Test-data (blå kurve) og kontrol-data (rød kurve)</em>
</p>

Der er muligvis biases i data, da de maskiner der læser sekvenserne, kan producere støj. Derudover er der altid mulighed for menneskelige fejl. Vi vil dog antage, at biases er minimale, da sekvenserne er blevet læst på et laboratorium i et kontrolleret miljø.
<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/89b5271e-b39c-4069-a133-87aa2563c0a9">
    <br>
    <em>Figur 7 - Autocorrelationsplot for vores observerede variable</em>
</p>

### Preprocessering

For at generalisere generne så vi kan fitte modeller på dem, har vi brug for at normalisere vores data. Her starter vi med at lægge 1 til alle observationer, så vi stadig kan bruge log-transformering, selv når coverage er 0. Siden DNA-strenge både kan læses forlæns og baglæns, var vi også nødt til at tage i betragtning, hvilken retning vores datapunkter var.

```{r}
if (strand_sign == "-"){
    GOI_data = GOI_data[order(GOI_data$index, decreasing=TRUE), ]
  }
log2_GOI_data = log2(GOI_data + 1)
```

Udover dette, skulle vi også finde kroppen af generne, hvilket er beskrevet nærmere i næste afsnit. Når vi kender kroppen, kan vi normalisere dem ved at trække medianen af kroppen fra for de enkelte gener. På denne måde sikrer vi en mere robust sammenligning af kontrol- og sample-gener.

```{r}
body.norm.factors = apply(log2_GOI_data[TSS_row:TES_row,], 2, median)
log2_GOI_data.bodynorm = t(t(log2_GOI_data) - body.norm.factors)
```
Til sidst trækker vi sample-genet fra kontrol-genet, så vi får differencen imellem dem. Således burde vi få en kurve, der viser defekten i sample-genet, og vi er derfor klar til at modellere på terminerings-defekten. I figur 8 vises dette med et gen.

```{r}
log2_GOI_data.bodynorm.ctrlsubtract = log2_GOI_data.bodynorm
  
for (i in 1:length(ctrl)) {
  log2_GOI_data.bodynorm.ctrlsubtract[, i + length(ctrl)] = log2_GOI_data.bodynorm.ctrlsubtract[, i + length(ctrl)] - log2_GOI_data.bodynorm.ctrlsubtract[, i]
  log2_GOI_data.bodynorm.ctrlsubtract[, i] = log2_GOI_data.bodynorm.ctrlsubtract[, i] - log2_GOI_data.bodynorm.ctrlsubtract[, i]
}
```

<p>
    <img width="750" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/87d9bf82-1532-49b2-a4d5-4f1fffe17324">
    <br>
    <em>Figur 8 - (a) Ingen normalisering → (b) Log-transformering → (c) Median fratrukket → (d) Kontrol fraktrukket</em>
</p>


## Identificering af transskriptions-sites

Det er en essentiel del af hele vores pipeline at bestemme hvilken region (dvs. start- & slutpunkt), vi antager transskriptionen foregår i, og derved vil bestemmelsen af denne danne grundlag for det data, vi benytter i den senere modellering. 

For hvert gen har vi typisk flere transskript-varianter, angivet i annoteringen. Derfor er det væsentligt at afgøre hvilket Transscript-Start-Site (TSS) og Transscript-End-Site (TES), vi vil benytte i vores analyse af et gen. Vi har derfor konstrueret en algoritme, der skal bestemme den region (TSS:TES), hvor transskriptions signalet er mest fremtrædende. Med andre ord ønsker vi at identificere de grænser, der med størst sandsynlighed afspejler det fulde transskript.

Dette opnår vi ved at kigge på hhv. stigning og fald i signal i hhv. opstrøms- og nedstrømgsregionen i nærheden af de mulige TSS'er og TES'er. I store træk benytter vi den TSS med umiddelbart størst stigning i signal og den TES med umiddelbart største fald i signal.
Desuden tager algoritmen højde for, hvornår der starter et nyt gen, så 'søgningen' efter en TSS og en TES kun udfolder sig inden for det aktuelle GOI's domæne. 

De to figurer forneden (figur 9. & 10.) demonstrerer, hvilke mulige TSS'er og TES'er algoritmen evaluerer for genet 'GAPDH' og hvilken TSS og TES der vælges. Den røde kurve er vores kontrol-data mens den blå kurve udgør vores prøve.
<table>
  <tr>
    <td>
      <img width="450" alt="TSS & TES kandidater" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/128427973/e37d25c7-b09a-46a4-930c-98d13b8d16e9">
      <br>
      <em>Figur 9 - TSS & TES kandidater</em>
    </td>
    <td>
      <img width="450" alt="Endelig TSS & TES" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/770cb179-25c2-4302-b64e-1f49739899ff">
      <br>
      <em>Figur 10 - Endelig TSS & TES</em>
    </td>
  </tr>
</table>

## Modellering

### HMM

Som beskrevet i vores databeskrivelse, arbejder vi med sekventielle data i form af nukleotidbasepar langs en DNA-streng. I sekventielle data gælder antagelsen om identisk og uafhængig fordeling (i.i.d) ikke. Vi udnytter de sekventielle mønstre som korrelation mellem nærtliggende observationer. Derfor giver det mening for os at anvende Markov-modeller, hvor der er en antagelse om, at fremtidige forudsigelser kun afhænger af de mest nylige observationer. Da vi forsøger at skelne mellem to tilstande, defekt og ikke-defekt i termineringen, introducerer vi en Hidden Markov Model (HMM). Her repræsenterer de to tilstande vores skjulte "states", mens vores observerede data er forskellen mellem vores transformerede kontrol- og sampledata. Modellen ses grafisk i figur 11 samt et gitter-diagram i figur 12, som repræsenterer mulige stier af skjule states.
```Denne model udspringer fra teorien i Cristopher M. Bishop, Pattern Recognition and Machine Learning, 2006, kap 13.```

Vi starter med at definere en HMM med to states og en emissionsmodel, der følger en independent gaussian fordeling med antagelsen om delt kovariansmatrix. Det vil sige at de to states, 2 = defekt og 1 = ingen defekt, hver især følger en normalfordeling, som er uafhængig af den anden men med samme varians. Så vi har:

- Antal tilstande: $`K = \{1, 2\}`$
- Skjulte variable: $`Z = \{z_1, z_2, ..., z_N\}`$
- Observerede variable: $`X = \{x_1, x_2, ..., x_N\}`$

<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/fd64f60d-6cce-44ec-b18e-20a2d79d0139">
    <br>
    <em>Figur 11 - Hidden Markov Model</em>
</p>

Vi initierer også parametrene til modellen: $`\theta = \{\pi, A, \Phi\}`$. Her initieres $\pi$ og A uniformt og parametrene til normalfordelingen vha. K-means:
- Initielle sandsynligheder: $`\pi = \{\pi_1, \pi_2\}, \pi_k = p(z_{1k} = 1) = \frac{1}{\vert K \vert}`$
- Overgangs matrix: $`A = \{a_{kj}\},`$ hvor $` a_{kj} = p(z_{n+1} = j\vert z_n = k) = \frac{1}{\vert K \vert}`$
- Emissions parametre: $`\Phi = \{\phi_k(x_n) \sim \mathcal{N}(\mu_k, \Sigma)\}`$

for $`k \in 1, ..., K`$ tilstande og $`n \in 1, ..., N`$ tidspunkter.

<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/bbeb49e7-80a1-424c-ac9f-f9254dd570a1">
    <br>
    <em>Figur 12 - Gitter-diagram for skjulte "states"</em>
</p>


Dette udførte vi i r vha. funktioner fra STAN pakken. Vi brugte koden:

```{r}
hmm = initHMM(curr.data.list, nStates=2, "IndependentGaussian", sharedCov=TRUE)
```

Nu hvor vi har initieret modellens parametre, kan vi fitte den til vores data ved at bruge maximum likelihood. Vi har den følgende simultane fordeling: <br>
$`p(X, Z\vert \theta) = p(z_1\vert\pi) [ \prod_{n=2} p(z_n \vert z_{n-1}, A) ] \prod_{m=1} p(x_m \vert z_m, \Phi)`$ <br>
Vi får så følgende likelihood-funktion ifgl. loven om total sandsynlighed: <br>
$`L(X,\theta) = p(X\vert\theta)=\sum_Z{p(X,Z\vert\theta)}`$ <br>
Denne maksimeres med EM-algoritmen og den rekursive forward-backward algoritme, ved at opdatere parametrene, således denne sandsynlighed er konvergeret, eller algortimen er kørt igennem 50 iterationer. Vi har brugt følgende kode til dette:
  
```{r}  
hmm_fitted = fitHMM(curr.data.list, hmm, maxIters=50)
```

Nu har vi identificeret de mest sandsynlige stadier individuelt for sekvensen. Vi er dog opmærksomme på, at defekter i transkriptionstermineringen ikke kan genopstå. Derfor kan stadiet praktisk talt kun skifte fra 1 til 2 én gang og tilsvarende kun skifte tilbage fra 2 til 1 én gang. Vi er derfor nødt til at anvende en mere robust metode for at finde den mest sandsynlige sekvens af skjulte stadier. Til dette bruger vi Viterbi-algoritmen til decoding, som finder den sti af skjulte tilstande, som maksimerer $`p(z_1, z_2, ..., z_N \vert x_1, x_2, ..., x_N)`$:

```{r}
viterbi = getViterbi(hmm_fitted, curr.data.list)
states = as.integer(viterbi[['GOI']])
```

På denne måde har vi nu en måde at bestemme hvorvidt, der er defekt i genet ved, at se om sekvensen af gentranskriberingen på noget tidspunkt er i tilstand 2.

### (Double)-Sigmoidal fitting

Efter vi har fundet det mest sandsynlige defekte område vha. HMM, kan vi nu fitte enten en sigmoidal eller double sigmoidal model over dette område. Dette hjælper os med, at være i stand til at kvantificere samt beskrive defekten med flere forskellige mål.

Typisk vil vi gøre brug af sigmoidal, hvis der starter et nyt gen kort efter termineringen af genet. Dvs. signalerne fra de to gener flyder sammen - ellers foretrækker vi at benytte double-sigmoid. 

Grunden til, at vi anvender en sigmoidal eller double-sigmoidal model til at beskrive read-through området, er fordi denne type kurver harmonerer godt med vores antagelse om, at de skjulte tilstande i modellen udelukkende kan bevæge sig fra 1 til 2 og tilbage til 1. Denne egenskab stemmer overens med den gradvise overgang, som en double-sigmoid kurve repræsenterer, fra lav til høj værdi og omvendt. Derudover giver muligheden for at vælge mellem en sigmoidal eller double-sigmoidal kurve os fleksibilitet i vores analyse, da det tillader, at vi kan analysere transskriptioner, hvor signalet flyder sammen med signal fra andre gener. Modellen er vist grafisk på et gen i figur 13.

Vi benytter fitAndCategorize-funktionen fra Sicegar pakken, til at afgøre om vi kan fitte en sigmoidal, double-sigmoidal eller ingen. Måden denne funktion fungerer er at forsøge at fitte både sigmoidal og double-sigmoidal kurver til vores data, som er den normaliserede difference mellem kontrol- og sample genet. Den estimerer altså parametrene A, B som er maksimal signalstyrke, $`h_1, h_1`$ som er vækstrater og $`m_1, m_2`$ som er midtpunkter for den første og anden sigmoid til hhv.:

Sigmoidal: $`\frac{A}{1+e^{-x}}`$,

Double-sigmoidal: $`\frac{A}{1+e^{-h_1(x-m_1)}} + \frac{B}{1+e^{-h_2(x-m_2)}}`$.

<p>
    <img width="450" alt="transkription" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/49984447/1b561a4f-f004-4da0-8905-fdbc1cfa97fd">
    <br>
    <em>Figur 13 - Fitted Double Sigmoid-kurve vist i lilla</em>
</p>


Hvis ingen af modellerne passer til vores data (hvilket undersøges vha. forskellige kriterier), returnerer funktionen "ambiguous". Hvis kun den ene passer, returnerer den navnet på denne, og hvis begge passer, vælger funktionen den model med laveste AIC score, hvilket er et kriterie, som bruges til, at bestemme hvilken machine learning model er bedst for et givent datasæt. Hvis fitAndCategorize returnerer 'ambigious', fitter vi både en sigmoidal og en double-sigmoidal, hvorefter vi vælger den med højest $`R^2`$.

Herefter bruger vi bruger vi enten doublesigmoidalFitFormula eller sigmoidalFitFormula til, at forudsige signalstyrke for området med defekt ud fra de givne (double) sigmoidal parametre.

Disse fits og forudsigelser benyttes til vores endelige resultater, som er vores outputs i vores readthrough-funktion. Nedenfor illustreres disse resultater for nogle et lille udpluk af generne. De forskellige attributter beskrives i dokumentet "RESULTATER.md".

|        | TSS|   TES|  body_diff|    us_diff|us.TSS_diff |rt    | max_rt_length|rt_TES  | rt_start| rt_end|    rt_int|      rt_sum|    rt_max|dsfit |sfit  | rt_end_fitted|extrapolated | best_fit_asymp| best_fit_Rsq| endDeclinePoint_x| rt_int_fitted| rt_sum_fitted|
|:-------|---:|-----:|----------:|----------:|:-----------|:-----|-------------:|:-------|--------:|------:|---------:|-----------:|---------:|:-----|:-----|-------------:|:------------|--------------:|------------:|-----------------:|-------------:|-------------:|
|SAMD11  |   1| 20652| -0.1720951|  0.9637476|NA          |TRUE  |         16009|944574  |    20653|  36661| 0.1720951|  8154.55281| 0.4458650|TRUE  |FALSE |         31491|FALSE        |      0.4346400|    0.3124570|             55075|     0.4974774|     7505.4771|
|KLHL17  |   1|  5136|  0.5820985|  3.0724062|NA          |TRUE  |           762|965719  |     5137|   5587| 0.0113300|   -18.06437| 0.2408704|FALSE |TRUE  |          5898|FALSE        |             NA|   -0.0389367|              5898|    -0.0698842|     -117.7091|
|PLEKHN1 |   1|  9384|  0.1013182| -1.4288504|NA          |TRUE  |         25272|975865  |     9385|  27328| 0.7247202| 12489.52924| 1.7907342|FALSE |TRUE  |         34656|FALSE        |             NA|    0.0050567|             34656|     0.6741845|    16691.0026|
|ISG15   |   1| 13403|  1.2747376|  2.3949412|NA          |TRUE  |          5579|1014540 |    13404|  18743| 0.4583741|  2575.88387| 1.2195022|TRUE  |FALSE |         18839|FALSE        |      0.0000000|    0.1698744|             37455|     0.4359586|     2538.1540|
|AGRN    |   1| 35999| -0.3976053|  0.5741375|NA          |TRUE  |        176118|1056118 |    36000| 212054| 0.3976053| 71589.56775| 1.5671740|TRUE  |FALSE |        105483|FALSE        |      0.3975892|    0.1484812|            148864|     0.3976053|    26887.4265|
|B3GALT6 |   1|  2805|  0.1169609|  4.6357742|NA          |FALSE |             0|0       |        0|      0| 0.0000000|     0.00000|        NA|FALSE |FALSE |             0|FALSE        |             NA|           NA|                NA|     0.0000000|        0.0000|


## Dokumentation
I mappen "R_Kode" ligger dokumentet "KODEOVERSIGT.md", som angiver bidragsyderne til de forskellige kodestykker og indeholder en kort beskrivelse af de enkelte kodestykkers funktionalitet.

I mappen "Demo" ligger dokumentet "RESULTATER.md", som præsenterer resultaterne fra vores analyse, udført med readthrough_analysis scriptet. Dokumentet inkluderer en beskrivelse af vores dataset, relevante bemærkninger samt en visuel demonstration af scriptets analyse af tre forskellige gener.

Dokumenterne giver også en forklaring på resten af indholdet i hver af de respektive mapper.

## Evaluering
Resultatet af vores projekt er et script, der specifikt analyserer og kvantificerer defekter i transkriptionens terminering i et givent gen. Scriptet er generaliseret og kræver derfor ikke meget mere end en kontrol- og en prøveindgang. Overordnet set består scriptet af pre-processing, en algoritme til identifikation af transkriptionens start- og slutpunkter, estimater af HMM-parametre, tilstands-annotering ved hjælp af Viterbi-dekodning og fitting af en sigmoidal eller dobbelt-sigmoidal kurve på defektens signal. Til sidst returnerer scriptet statistik, der beskriver eventuelle defekter. Scriptet, eller værktøjet, er i høj grad automatiseret og strømlinet, hvilket gør det mere tilgængeligt som et generelt værktøj til analyse af transskriptionelle defekter.

Dette stemmer i høj grad overens med det mål, vi beskrev i projektbeskrivelsen. Det er dog værd at nævne, at nøjagtigheden af vores samlede model umiddelbart ikke er tilfredsstillende (årsagen fremgår i RESULTATER.md under afsnittet Bemærkning). 

Undervejs i projektet blev muligheden for at benytte et neuralt netværk til at påvise read-throughs diskuteret. Vi valgte dog ikke at fortsætte ad denne vej, da scriptet oprindeligt kun statisk beskrev de påviste read-throughs, og vi manglede derfor inputdata til prøver, hvor der ikke var read-throughs. For at løse dette problem ville vi være nødt til at re-designe analysen, som scriptet er baseret på. Dette vurderede vi som værende for omfattende med vores tidsramme i betragtning.


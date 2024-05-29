# Data


Vi har, vha. vores wrapper samt en shell script fil, kørt readthrough_analysis scriptet på genomeDK clusteret og analyseret transskription i 10.386 forskellige gener - to udgaver af hver: den ene udgave med termineringsdefekt og den anden udgave uden termineringsdefekt. Dvs. 20.772 kørsler af scriptet.

### Dataset
Dette har resulteret i et samlet dataset på 20.772 observationer. 

I sample_dataset.md filen ligger der 1.000 tilfældigt samplet observationer fra det fulde dataset. 

I dataset.csv filen ligger det fulde dataset i csv-format. 

## Note
Efter at have undersøgt vores dataset, mistænker vi, at readthrough_analysis scriptets analyse er for følsom. Ud fra det vi ved om de 2 gange 10.386 prøver, der er analyset i scriptet, vil vi ca. forvente 50/50 fordeling mellem defekt/ikke-defekt, og i realiteten ser vi en fordeling nærmere 84% defekt mod 16% ikke-defekt.

I kodestykket forneden, har vi udregnet middelværdien af rt_int_fitted (readthrough intensitet (gennemsnitlig eller median signal) i det fitted påviste readthrough) for de observationer, hvor readthrough_analysis påviser en termineringsdefekt for både det dataset, hvor vi forventer en termineringsdefekt (n = 9940) og det hvor vi ikke forvetner en termineringsdefekt (n = 7409). 
```{r}
> mean(unlist(df_positive_rt_true$rt_int_fitted))
[1] 1.10093
> mean(unlist(df_negative_rt_true$rt_int_fitted))
[1] 0.5262786
```
Det fremgår tydeligt, at der er en signifikant forskel i middelværdierne, og det er altså af denne grund, at vi konkluderer, at analysen er for følsom.

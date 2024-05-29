# Resultater af Analyser


Vi har, vha. vores wrapper samt en shell script fil, kørt readthrough_analysis scriptet på genomeDK clusteret og analyseret transskription i 10.386 forskellige gener - to udgaver af hver: den ene udgave med termineringsdefekt og den anden udgave uden termineringsdefekt. Dvs. 20.772 kørsler af scriptet.

### Dataset
Dette har resulteret i et samlet dataset på 20.772 observationer. 

I sample_dataset.md filen ligger der 1.000 tilfældigt samplet observationer fra det fulde dataset. 

I dataset.csv filen ligger det fulde dataset i csv-format. 

## Bemærkning
Efter at have undersøgt vores dataset mistænker vi, at readthrough_analysis scriptets analyse er for følsom. Ud fra vores viden om de 20.772 prøver, der er blevet analyseret i scriptet, ville vi forvente en næsten lige fordeling mellem defekt og ikke-defekt, med omkring 50/50. I virkeligheden observerer vi imidlertid en fordeling på ca. 84% defekt mod 16% ikke-defekt.

I kodestykket nedenfor har vi beregnet middelværdien af rt_int_fitted (median-signalet i det påviste read-through) for observationerne, hvor readthrough_analysis påviser en termineringsdefekt både i datasættet, hvor vi forventer en termineringsdefekt (n = 9940), og i det, hvor vi ikke forventer en termineringsdefekt (n = 7409).
```{r}
> mean(unlist(df_positive_rt_true$rt_int_fitted))
[1] 1.10093
> mean(unlist(df_negative_rt_true$rt_int_fitted))
[1] 0.5262786
```
Det fremgår tydeligt, at der er en signifikant forskel i middelværdierne, og det er derfor, at vi konkluderer, at analysen er for følsom.

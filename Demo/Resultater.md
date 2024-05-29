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

## Demonstration
I tabellen nedenfor præsenteres tre figurer, der repræsenterer tre gener, som vi har analyseret ved hjælp af vores readthrough_analysis scriptet. Hver figur indeholder to plots: det øverste plot viser det normaliserede signal fra prøven sammen med det normaliserede signal fra kontrollen, mens det nederste plot viser forskellen mellem det normaliserede prøvesignal og det normaliserede kontrolsignal.

TSS & TES er markeret med de sorte prikkede vertikale linjer, mens de lyserøde stiplede vertikale linjer angiver området for vores (double)-sigmoidal fitting. Den fede lyserøde horizontale streg under plottet angiver det område, hvor vores HMM befinder sig i tilstand 2. Den lyserøde stiplede horizontale streg under nogle af plotsne angiver, at vores HMM er tilbage i tilstand 1.
<table>
  <tr>
    <td>
      <img width="450" alt="TSS & TES kandidater" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/128427973/a028a564-4fca-4cce-90a1-6e3244d0ce1c">
      <br>
      <em>Figur 1 - HELLS</em>
    </td>
    <td>
      <img width="450" alt="Endelig TSS & TES" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/128427973/961ed17d-e911-4966-bb6e-e6b3313dffa0">
      <br>
      <em>Figur 2 - NOP56</em>
    </td>
    <td>
      <img width="450" alt="Endelig TSS & TES" src="https://github.com/Kasperlanghoff12/Dataprojekt/assets/128427973/0e794841-fa82-4983-9f33-41d851277b25">
      <br>
      <em>Figur 3 - RPL10</em>
    </td>
  </tr>
</table>

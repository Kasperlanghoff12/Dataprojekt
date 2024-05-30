# Resultater af Analyser


Vi har, vha. vores wrapper samt en shell script fil, kørt readthrough_analysis scriptet på genomeDK clusteret og analyseret transskription i 10.386 forskellige gener - to udgaver af hver: den ene udgave med termineringsdefekt og den anden udgave uden termineringsdefekt. Dvs. 20.772 kørsler af scriptet.

### Dataset
Dette har resulteret i et samlet dataset på 20.772 observationer. 

I sample_dataset.md filen ligger der 1.000 tilfældigt samplet observationer fra det fulde dataset. 

I dataset.csv filen ligger det fulde dataset i csv-format. 

Tabellen for neden viser en oversigt over de 22 variable i datasettet samt en forklaring af hver variables betydning.
| Variabel | Forklaring |
|----------|----------|
| TSS    | Transscription Start Site   |
| TES    | Transscription End Site   |
| body_diff    | Gns. forskel mellem norm. krop for prøve og kontrol  |
| us_diff    | Gns. forskel mellem opstrøms-region og krop for prøve  |
| us.TSS_diff	    | Alternativ til us_diff, bruges sjældent    |
| rt    | Read-through (transskriptionel defekt) - TRUE/FALSE   |
| max_rt_length    | Længde af rt-region   |
| rt_TES    | TES hvorfra rt starter - bestemt af HMM   |
| rt_start    | Read-through start (i forhold til uTSS) bestemt af HMM - bruges til sigmoidal eller dobbelt-sigmoidal fitting   |
| rt_end   | Read-through slut (i forhold til uTSS) bestemt af HMM  |
| rt_int   | Read-through intensitet (middel- eller medianværdi af signalet) i det HMM bestemte read-through  |
| rt_sum   | Read-through intensitet integreret i det HMM bestemte read-through |
| rt_max   | Maximale read-through intensitet (middel- eller medianværdi af signalet) i det HMM bestemte read-through  |
| dsfit   | Double-Sigmoidal fitted - TRUE/FALSE  |
| sfit   | Sigmoidal fitted - TRUE/FALSE  |
| rt_end_fitted   | Afslutning på read-through (relativ til uTSS) bestemt gennem Sigmoidal eller double-Sigmoidal model  |
| extrapolated   | Er rt_end_fitted bestemt ved ekstrapolering ud over dataområdet for det nuværende GOI - TRUE/FALSE|
| best_fit_asymp   | Asymptotisk værdi for den fittede dobbelt-Sigmoid kurve (bruges til beregning af rt_end_fitted) (kun relevant for dobbeltsigmoidal tilpasning)  |
| best_fit_Rsq   | $R^2$ værdi for den bedst fittede kurve (kun relevant for Sigmoidal eller double-Sigmoidal fitting)  |
| endDeclinePoint_x   | Endepunkt for den tilpassede dobbelt-Sigmoid kurve (bruges til at reducere de gemte data, hvis det er muligt) (kun relevant for dobbelt-Sigmoid fitting)  |
| rt_int_fitted   | Read-through intensitet (mean eller median signal) i det fittede påviste read-through  |
| rt_sum_fitted   | Read-through intensitet integreret i det fittede påviste read-through  |


## Bemærkning
Efter at have undersøgt vores dataset mistænker vi, at readthrough_analysis scriptets analyse er for følsom. Ud fra det, vi på forhånd ved om de 20.772 prøver, der er blevet analyseret i scriptet, ville vi forvente en næsten lige fordeling mellem defekt og ikke-defekt - altså omkring 50% defekt og 50% ikke-defekt. I virkeligheden observerer vi imidlertid en fordeling på ca. 84% defekt mod 16% ikke-defekt.

I kodestykket nedenfor har vi beregnet middelværdien af rt_int_fitted (median-signalet i det påviste read-through) for observationerne, hvor readthrough_analysis påviser en termineringsdefekt både i datasættet, hvor vi forventer en termineringsdefekt (n = 9940), og i det, hvor vi ikke forventer en termineringsdefekt (n = 7409).
```{r}
> mean(unlist(df_positive_rt_true$rt_int_fitted))
[1] 1.10093
> mean(unlist(df_negative_rt_true$rt_int_fitted))
[1] 0.5262786
```
Det fremgår tydeligt, at der er en signifikant forskel i middelværdierne, og det er altså derfor, at vi mistænker, at analysen er for følsom.

## Demonstration
I tabellen nedenfor præsenteres tre figurer, der repræsenterer tre gener (gennavne: HELLS, NOP56 & RPL10), som vi har analyseret ved hjælp af vores readthrough_analysis scriptet. Hver figur indeholder to plots: det øverste plot viser det normaliserede signal fra prøven sammen med det normaliserede signal fra kontrollen, mens det nederste plot viser forskellen mellem det normaliserede prøvesignal og det normaliserede kontrolsignal.
### Visuelle Indikatorer i Plots
I plottet er TSS og TES markeret med tynde, sorte, prikkede vertikale linjer, mens de lyserøde, stiplede vertikale linjer indikerer området for vores (dobbelt)-sigmoidale fitting. En fed, lyserød, horisontal streg under plottet indikerer det område, hvor vores Hidden Markov Model (HMM) befinder sig i tilstand 2, mens de lyserøde, stiplede, horisontale streger under nogle af plottene angiver, at vores HMM er vendt tilbage til tilstand 1.
### Fortolkning
I alle tre analyser er resultatet det samme: en påvist termineringsdefekt. Dette er tydeligt demonstreret i de nederste plots, hvor den dobbelte sigmoid-kurve viser, at prøven fortsat udsender et betydeligt signal, efter at transskriptionen skulle være termineret.
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

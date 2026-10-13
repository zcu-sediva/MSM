# Zadání úlohy 2 — KMA/MSM 2026/2027

**Téma:** Testy o poloze — jednovýběrové, dvouvýběrové a vícevýběrové testy, ANOVA, vícerozměrné testy

**Zadáno:** 3. 11. 2026 (přednáška 7) &nbsp;|&nbsp; **Odevzdání:** 24. 11. 2026

## Zadání

Cílem práce je demonstrovat znalost statistických metod pro analýzu polohy u jednorozměrných a vícerozměrných dat. Využijte parametrické i neparametrické přístupy a výsledky interpretujte.

1. **Volba dat**
   - Zvolte vícerozměrnou datovou sadu, která obsahuje alespoň dvě spojité proměnné a pozorování rozčleněná do alespoň tří skupin.
   - Věnujte pozornost dostatečnému počtu pozorování ve skupinách.
   - *Příklad: [Iris flower data set](https://en.wikipedia.org/wiki/Iris_flower_data_set) — délka a šířka kališního a okvětního lístku (4 proměnné) pro tři druhy kosatců, v každé skupině 50 pozorování. Další datové sady lze najít např. na [UCI Machine Learning Repository](https://archive.ics.uci.edu/) nebo [OpenML](https://www.openml.org/).*

2. **Jednovýběrové testy o poloze**
   - Pro každou spojitou proměnnou stanovte nulovou a alternativní hypotézu (např. test, zda se střední hodnota rovná zvolené referenční hodnotě) a zdůvodněte volbu referenční hodnoty (např. průměr všech dat, teoretická hodnota).
   - Proveďte test parametrickým i neparametrickým přístupem.
   - Uveďte p-hodnotu, rozhodnutí o hypotéze a interpretaci v kontextu dat; diskutujte rozdíly mezi parametrickým a neparametrickým přístupem.
   - Doplňte vhodnou vizualizaci, včetně kontroly normality.
   - *Příklad: ověřujeme, zda se průměrná délka kališního lístku (sepal length) u každého druhu rovná hodnotě μ₀ = 5,84 cm (průměr všech kosatců), a zda se medián délky okvětního lístku (petal length) u každého druhu rovná 4,35 cm (medián všech kosatců).*

3. **Dvouvýběrové testy o poloze**
   - Pro každou spojitou proměnnou porovnejte polohu vždy dvou skupin mezi sebou.
   - Postupujte obdobně jako v bodě 2 (hypotézy, parametrický i neparametrický test, interpretace).
   - Výsledky vizualizujte, např. krabicovými grafy po skupinách, grafem intervalů spolehlivosti pro rozdíly a maticí (heatmapou) p-hodnot atd.
   - *Příklad: porovnáváme průměry délek a šířek kališních i okvětních lístků mezi jednotlivými druhy (vždy dva druhy mezi sebou) dvouvýběrovým t-testem. Výsledky zobrazíme ve formě matice p-hodnot a grafem s intervaly spolehlivosti.*

4. **Vícevýběrové testy o poloze**
   - Pro každou spojitou proměnnou otestujte shodu polohy ve všech skupinách současně (např. ANOVA, neparametricky Kruskalův–Wallisův test).
   - Proveďte post-hoc analýzu (mnohonásobné porovnávání) a výsledky porovnejte s testováním „po dvojicích“ z bodu 3.
   - Výsledky vizualizujte, např. krabicovými nebo houslovými grafy všech skupin, grafem simultánních intervalů spolehlivosti z post-hoc analýzy a diagnostickými grafy reziduí (normalita, shoda rozptylů) atd

5. **Vícerozměrný přístup**
   - Analogicky bodům 2–4 porovnejte polohu vícerozměrných dat.
   - Výsledky vizualizujte, např. maticí bodových grafů obarvenou podle skupin, konfidenčními elipsami podle Mahalanobisovy vzdálenosti, grafem simultánních intervalů spolehlivosti pro vektory středních hodnot a diagnostickými grafy reziduí (normalita, shoda kovariančních matic) atd
   - *Příklad: párové porovnání vektorů středních hodnot (Hotellingův T² test), porovnání kovariančních matic, porovnání vektorů středních hodnot všech skupin (MANOVA) atd.*

## Software a odevzdání

- Úlohu lze zpracovat v **Pythonu** nebo **MATLABu**; jiný software je možný jen po předchozí domluvě s vyučujícími.
- Rozsah a forma odevzdávané práce (poster + kód) se řídí [obecnými požadavky na semestrální práce](../README.md).
- Do STAGu odevzdává každý student/ka **poster i kód** v adresáři `Prijmeni_jmeno_02`; u dvojic oba stejnou verzi.

---
*Blanka Šedivá, KMA/MSM 2026/2027*

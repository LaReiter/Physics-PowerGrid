### Lars Reiter Nielsen
### Bachelorprojekt 2016 i Fysik
### Matlab scripts til projektet: Stability of power grid topologies
### Lars Reiter Nielsen

* BEMÆRK: RK4_N_Node scriptsne er MEGET tunge. Alle SimScript_BraessParadox scriptsne arbejder over en hel dag, 
afhængigt af parameter valg.

***
AuxFunctions indeholder hjælpefunktionerne:

NnodeSim:   Generelt script der bruger RK4 metoden til at finde samtlige numeriske løsninger
til kuramatos model under udvalgte begyndelsesbetingelser. Finder desuden kraft vektor samt faseordens
parameter.

NnodeSim_CritVal:  Givet topologi specifikationer, returnerer langstids gennemsnitlige faseordens parameter.

CreateAdj:  Givet antal knuder, forbindelser og et "udgangspunkts koblingsmatrix", returnerer
koblingsmatrix med valgte specifikationer.

CreateAdj2:  Giver mulighed for at have forskellig koblingsstyrke selvom udgangspunkt koblingsmatrix
er givet. "Udgangspunkt" koblingsmatrix kan også være forskellig dimension fra antallet af knuder.

findKc: Givet topologi specifikationer, returnerer den kritiske koblingsstyrke samt det langstids
gennemsnitlige faseordens parameter.


***
RK4_two_node indeholder scriptsne:

Two_node_analysis: Her blev undersøgt robustheden af to-knude netværket under en pertubation af effekten. Svarende til
et øget forbrug.

Two_node_Analysis_dissipation: Her finder vi den kritske koblingsstyrke K_c der udviser stabil system opførsel, alt imens
systemet perturberes.

Two_node_Analysis_Kc: Her bestemmer vi den kritiske koblingsstyrke.

Two_node_Analysis_Kc_StableFixPoint: Her redegøres for at der findes et stabilt og et ikke
stabilt fixpunkt (som løber løbsk!) i to knude topologien.

***

RK4_five_node indeholder scriptsne:

Find_Opt_5Node:   Dette script blev brugt til at konstruere randomiserede 5 knude netværk der udviser Braess Paradox,
såvel på aktuelle som på nye forbindelser.

five_node_topologi_BraessParadox og five_node_topologi_BraessParadox:  Disse script undersøger braess paradox på
to udvalgte 5 knude netværk.

five_node_topologi_CriticalValue:   Dette script bestemmer den kritiske koblingsstyrke på 5 knude netværket.

oscil:   Hjælpe funktion benyttet til at udregne sin(theta_diff) delen i modellen

Braess_Paradox_replicate_8Nodes:    En replikation af en topologi fra "Braess paradox in oscillator networks, desynchronization and power outage",
der udviser Braess Paradox. s. 4 i artiklen.

***

RK4_N_node indeholder scriptsne:

SimScript_BraessParadox_LinkStr:  Finder Braess Paradox på aktuelle forbindelser af randomiserede 30-node topologier
ved at forstærke koblingsstyrken.

SimScript_BraessParadox_NewLinks: Finder Braess Paradox på potentielle forbindelser af randomiserede 30-node topologier
ved at tilføje nye forbindelser.

SimScript_PowervsTopology:  Undersøger hvorvidt der en gennemgående kritiske koblingsstyrke for randomiserede 30 node topologier

SimScript_BraessParadox_IncComplexity:  Tilføjer n = 2,4 ... 20 tilfældigt placerede generatorer (med tilfældige forbindelser)
til netværket og undersøger hyppigheden af Braess Paradox på aktuelle forbindelser

SimScript_TopologyStabilityVSGenLinks:    Undersøger hvordan stabilitetn afhænger af antal forbindelser der udgår fra generatorne.


***

Kontakt: lenny9876@gmail.com 
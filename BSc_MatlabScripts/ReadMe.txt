### Lars Reiter Nielsen
### Bachelorprojekt 2016 i Fysik
### Matlab scripts til projektet: Stability of power grid topologies
### Lars Reiter Nielsen

* BEM�RK: RK4_N_Node scriptsne er MEGET tunge. Alle SimScript_BraessParadox scriptsne arbejder over en hel dag, 
afh�ngigt af parameter valg.

***
AuxFunctions indeholder hj�lpefunktionerne:

NnodeSim:   Generelt script der bruger RK4 metoden til at finde samtlige numeriske l�sninger
til kuramatos model under udvalgte begyndelsesbetingelser. Finder desuden kraft vektor samt faseordens
parameter.

NnodeSim_CritVal:  Givet topologi specifikationer, returnerer langstids gennemsnitlige faseordens parameter.

CreateAdj:  Givet antal knuder, forbindelser og et "udgangspunkts koblingsmatrix", returnerer
koblingsmatrix med valgte specifikationer.

CreateAdj2:  Giver mulighed for at have forskellig koblingsstyrke selvom udgangspunkt koblingsmatrix
er givet. "Udgangspunkt" koblingsmatrix kan ogs� v�re forskellig dimension fra antallet af knuder.

findKc: Givet topologi specifikationer, returnerer den kritiske koblingsstyrke samt det langstids
gennemsnitlige faseordens parameter.


***
RK4_two_node indeholder scriptsne:

Two_node_analysis: Her blev unders�gt robustheden af to-knude netv�rket under en pertubation af effekten. Svarende til
et �get forbrug.

Two_node_Analysis_dissipation: Her finder vi den kritske koblingsstyrke K_c der udviser stabil system opf�rsel, alt imens
systemet perturberes.

Two_node_Analysis_Kc: Her bestemmer vi den kritiske koblingsstyrke.

Two_node_Analysis_Kc_StableFixPoint: Her redeg�res for at der findes et stabilt og et ikke
stabilt fixpunkt (som l�ber l�bsk!) i to knude topologien.

***

RK4_five_node indeholder scriptsne:

Find_Opt_5Node:   Dette script blev brugt til at konstruere randomiserede 5 knude netv�rk der udviser Braess Paradox,
s�vel p� aktuelle som p� nye forbindelser.

five_node_topologi_BraessParadox og five_node_topologi_BraessParadox:  Disse script unders�ger braess paradox p�
to udvalgte 5 knude netv�rk.

five_node_topologi_CriticalValue:   Dette script bestemmer den kritiske koblingsstyrke p� 5 knude netv�rket.

oscil:   Hj�lpe funktion benyttet til at udregne sin(theta_diff) delen i modellen

Braess_Paradox_replicate_8Nodes:    En replikation af en topologi fra "Braess paradox in oscillator networks, desynchronization and power outage",
der udviser Braess Paradox. s. 4 i artiklen.

***

RK4_N_node indeholder scriptsne:

SimScript_BraessParadox_LinkStr:  Finder Braess Paradox p� aktuelle forbindelser af randomiserede 30-node topologier
ved at forst�rke koblingsstyrken.

SimScript_BraessParadox_NewLinks: Finder Braess Paradox p� potentielle forbindelser af randomiserede 30-node topologier
ved at tilf�je nye forbindelser.

SimScript_PowervsTopology:  Unders�ger hvorvidt der en gennemg�ende kritiske koblingsstyrke for randomiserede 30 node topologier

SimScript_BraessParadox_IncComplexity:  Tilf�jer n = 2,4 ... 20 tilf�ldigt placerede generatorer (med tilf�ldige forbindelser)
til netv�rket og unders�ger hyppigheden af Braess Paradox p� aktuelle forbindelser

SimScript_TopologyStabilityVSGenLinks:    Unders�ger hvordan stabilitetn afh�nger af antal forbindelser der udg�r fra generatorne.


***

Kontakt: lenny9876@gmail.com 
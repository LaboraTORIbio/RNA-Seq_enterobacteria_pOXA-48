library(ggplot2)

### BP related to iron in the main RNA-Seq dataset

BP_main <- read.csv('table_count_BPs_iron.tsv', header = TRUE, sep = "\t")

strain_order_1 <- c("C325", "EC10", "MG1655", "KPN04", "KPN07", "KPN10")
func_order_1 <- c("response to iron(III) ion ", "heme oxidation ", "heme metabolic process ", "siroheme biosynthetic process ", "heme O biosynthetic process ", "heme biosynthetic process ", "intracellular iron ion homeostasis ", "intracellular sequestering of iron ion ", "protein maturation by iron-sulfur cluster transfer ", "iron incorporation into metallo-sulfur cluster ", "[2Fe-2S] cluster assembly ", "iron-sulfur cluster assembly ", "heme transmembrane transport ", "iron ion transmembrane transport ", "iron ion transport ", "iron import into cell ", "siderophore transmembrane transport ", "siderophore transport ", "siderophore-dependent iron import into cell ", "siderophore biosynthetic process ", "enterobactin biosynthetic process ")
custom_colors_1 <- c("indianred2", "lightsalmon", "gray95",
                     "plum2", "#BA68C8", "mediumorchid4")

ggplot(BP_main, aes(fill=factor(Strain, level=strain_order_1), x=Count, y=factor(GO_BP, level=func_order_1))) +
  geom_bar(position="stack", stat="identity", color="black", size=0.3) +
  theme_bw() + scale_x_continuous(limits = c(-5,5), breaks = seq(-5, 5, by = 1)) +
  geom_vline(xintercept = 0, size = 1, size = 0.9) +
  scale_fill_manual(values = custom_colors_1)


### Gene Products related to iron the main RNA-Seq dataset

Prod_main <- read.csv('table_count_Products_iron.tsv', header = TRUE, sep = "\t")

strain_order_2 <- c("C325", "EC10", "MG1655", "CF13", "J57", "K147", "KPN15", "KPN18", "KPN04", "KPN07", "KPN10")
prod_order_1 <- c("superoxide dismutase [Fe] (sodB)", "Ni/Fe-hydrogenase b-type cytochrome subunit (hyaC)", "Ni/Fe-hydrogenase cytochrome b subunit (hybB)", "heme lyase NrfEFG subunit NrfG (nrfG)", "heme lyase NrfEFG subunit NrfF (nrfF)", "cytochrome c nitrite reductase Fe-S protein (nrfC)", "cytochrome c nitrite reductase pentaheme subunit (nrfB)", "heme lyase CcmF/NrfE family subunit", "heme exporter protein CcmB (ccmB)", "cytochrome c biogenesis heme-transporting ATPase CcmA (ccmA)", "protoheme IX farnesyltransferase (cyoE)", "protoheme IX biogenesis protein HemY (hemY)", "protein-methionine-sulfoxide reductase heme-binding subunit MsrQ (msrQ)", "iron-containing redox enzyme family protein", "iron-containing alcohol dehydrogenase", "BolA family iron metabolism protein IbaG (ibaG)", "succinate dehydrogenase iron-sulfur subunit SdhB (sdhB)", "fumarate reductase iron-sulfur protein (frdB)", "aldehyde dehydrogenase iron-sulfur subunit (paoA)", "Rieske 2Fe-2S domain-containing protein", "2Fe-2S ferredoxin-like protein", "YfhL family 4Fe-4S dicluster ferredoxin", "4Fe-4S dicluster domain-containing protein", "4Fe-4S binding protein", "(4Fe-4S)-binding protein", "Fe-S cluster assembly protein SufB (sufB)", "Fe-S cluster assembly scaffold SufA (sufA)", "Fe-S cluster assembly scaffold IscU (iscU)", "Fe-S cluster assembly transcriptional regulator IscR (iscR)", "iron-sulfur cluster assembly protein IscA (iscA)", "iron-sulfur cluster repair protein YtfE (ytfE)", "Fe-S biogenesis protein NfuA (nfuA)", "iron-sulfur cluster carrier protein ApbC (apbC)", "non-heme ferritin-like protein", "non-heme ferritin (ftnA)", "ferritin-like domain-containing protein", "bacterioferritin-associated ferredoxin (bfd)", "bacterioferritin (bfr)", "siderophore-interacting protein", "catecholate siderophore receptor Fiu", "bifunctional siderophore receptor/adhesin Iha (iha)", "TonB-dependent siderophore receptor", "iron uptake system protein EfeO", "iron ABC transporter substrate-binding protein", "iron export ABC transporter permease subunit FetB (fetB)", "iron ABC transporter ATP-binding protein FetA (fetA)", "Fe(3+)-hydroxamate ABC transporter substrate-binding protein FhuD (fhuD)", "Fe3+-hydroxamate ABC transporter ATP-binding protein FhuC (fhuC)", "Fe(3+)-hydroxamate ABC transporter permease FhuB (fhuB)", "heme ABC transporter permease", "iron-dicitrate ABC transporter permease FecC (fecC)", "fec operon regulator FecR (fecR)", "Fe(2+) transporter permease subunit FeoB (feoB)", "ferrous iron transporter A (feoA)", "iron-enterobactin ABC transporter ATP-binding protein (fepC)", "siderophore enterobactin receptor FepA (fepA)", "enterobactin transporter EntS (entS)", "enterobactin non-ribosomal peptide synthetase EntF (entF)", "enterobactin biosynthesis bifunctional isochorismatase/aryl carrier protein EntB (entB)")
custom_colors_3 <-c("indianred2", "lightsalmon", "gray95", "darkgoldenrod1", "#A5D6A7",
                    "lightblue1", "lightskyblue", "skyblue3", "plum2", "#BA68C8", "mediumorchid4")

ggplot(Prod_main, aes(fill=factor(Strain, level=strain_order_2), x=Count, y=factor(Product, levels=prod_order_1))) +
  geom_bar(position="stack", stat="identity", color="black", size=0.3) +
  theme_bw() + scale_x_continuous(limits = c(-5,5), breaks = seq(-5, 5, by = 1)) +
  geom_vline(xintercept = 0, size = 0.9) +
  scale_fill_manual(values = custom_colors_3)


### Gene Products related to iron the new RNA-Seq (DeltaLysR dataset)

Prod_new <- read.csv('table_count_Products_iron_new.tsv', header = TRUE, sep = "\t")
Prod_new$Count <- Prod_new$Count * -1  # to plot the inverse comparison

prod_order_2 <- c("heme-binding protein", "non-heme ferritin (ftnA)", "bacterioferritin (bfr)", "Ni/Fe-hydrogenase cytochrome b subunit (hybB)", "non-heme iron oxygenase ferredoxin subunit", "BolA family iron metabolism protein IbaG (ibaG)", "fumarate reductase iron-sulfur protein (frdB)", "YfhL family 4Fe-4S dicluster ferredoxin", "iron-sulfur cluster insertion protein ErpA (erpA)", "Fe-S cluster assembly protein IscX (iscX)", "Fe-S cluster assembly scaffold IscU (iscU)", "Fe-S cluster assembly transcriptional regulator IscR (iscR)", "iron-sulfur cluster assembly protein IscA (iscA)", "Fe-S biogenesis protein NfuA (nfuA)", "TonB-dependent siderophore receptor", "iron uptake system protein EfeO", "iron ABC transporter ATP-binding protein FetA (fetA)", "Fe(3+)-hydroxamate ABC transporter substrate-binding protein FhuD (fhuD)", "Fe3+-hydroxamate ABC transporter ATP-binding protein FhuC (fhuC)", "Fe(3+)-hydroxamate ABC transporter permease FhuB (fhuB)", "Fe(2+) transporter permease subunit FeoB (feoB)", "ferrous iron transporter A (feoA)", "Fe(3+)-siderophore ABC transporter permease (fepD)", "siderophore enterobactin receptor FepA (fepA)", "iron-enterobactin ABC transporter permease (fepG)", "iron-enterobactin ABC transporter ATP-binding protein (fepC)", "enterobactin transporter EntS (entS)", "enterobactin non-ribosomal peptide synthetase EntF (entF)", "enterobactin biosynthesis bifunctional isochorismatase/aryl carrier protein EntB (entB)", "ferric iron uptake transcriptional regulator (fur)")
custom_colors_2 <- c("darkgoldenrod1", "lightskyblue")

ggplot(Prod_new, aes(fill=Strain, x=Count, y=factor(Product, levels=prod_order_2))) +
  geom_bar(position="stack", stat="identity", color="black", size=0.3) +
  theme_bw() + scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, by = 1)) +
  geom_vline(xintercept = 0, size = 0.9) +
  scale_fill_manual(values = custom_colors_2)


# Tip Top Table

<!--
## A tool to evaluate the topological preference of a TMH based on a population of TMHs with known topology

This runs the TMH sequence through a positionally dependent matrix of residue scores and checks the total score between forwards/backwards runs of the TMH. A greater difference indicates a greater topological preference. The advantage of this method is that the individual contribution of each residue are calculated, and whilst the accuracy of the predictor may not always be the highest overall, it allows for sensitive evaluation of the topology of a TMP without the need for hidden layers in neural networks or HMMs.

Enter your input below. Note that this is not the full protein sequence, nor a fasta formatted sequence. The sequence should be the predicted, or experimentally derived TMH with ±5 flanking residues.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/JamesABaker/TMH-open-topology/blob/master/)
-->

## Input

UniProt ID

## Output

TMH positions, topology, evidence, Info on TMH ...

# Goals

## Minimum

-   [ ] Verify if in MPTOPO
      -  [x] TMH positions
      -  [ ] I/O topology
-   [ ] Verify if in TOPDB
      -  [x] TMH positions
      -  [ ] I/O topology
-   [ ] Verify if in UniProt
      -  [x] TMH positions
      -  [ ] I/O topology
-   [x] Record source evidence

## Useful

-   [ ] Pore residue score
-   [ ] Delta H for each TMHs
-   [ ] TipTop Score
-   [ ] Hydrophobicity

## Above and beyond

-   [ ] Include Better Predict TMHs Than UniProt's Transmem
-   [ ] Remove signal peptides
-   [ ] Implement in Django
-   [ ] Beta barrels

<!--

TMs In Protein TOPology = TiPTop
TYpical Protein TOPology = TypTop
TypIcal Protein TOPology = TipTop
Tip Top Protein Topology

 -->

# Back-end sources and references

 `topdb_all.xml` from <http://topdb.enzim.hu> Dobson, L., Langó, T., Reményi, I. & Tusnády, G. E. Expediting topology data gathering for the TOPDB database. Nucleic Acids Res. 43, D283–D289 (2015).

 `mptopoTblXml.xml` from <http://blanco.biomol.uci.edu/mptopo/> Jayasinghe, S. MPtopo: A database of membrane protein topology. Protein Sci. 10, 455–458 (2001).

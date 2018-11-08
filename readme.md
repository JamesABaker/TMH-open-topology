# TipTop Table

<!--
## A tool to evaluate the topological preference of a TMH based on a population of TMHs with known topology

This runs the TMH sequence through a positionally dependent matrix of residue scores and checks the total score between forwards/backwards runs of the TMH. A greater difference indicates a greater topological preference. The advantage of this method is that the individual contribution of each residue are calculated, and whilst the accuracy of the predictor may not always be the highest overall, it allows for sensitive evaluation of the topology of a TMP without the need for hidden layers in neural networks or HMMs.

Enter your input below. Note that this is not the full protein sequence, nor a fasta formatted sequence. The sequence should be the predicted, or experimentally derived TMH with ±5 flanking residues.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/JamesABaker/TMH-open-topology/blob/master/)
-->

## Input

UniProt ID

## Output

TMH positions, topology, evidence ...

# To Do

## Minimum

-   [ ] Verify if in MPTOPO and get TMH positions and I/O topology
-   [ ] Verify if in TOPDB and get TMH positions and I/O topology
-   [ ] Verify if in UniProt and get TMH positions and I/O topology
-   [ ] Record source evidence

## Useful

-   [ ] Pore residue?
-   [ ] Delta H for each TMHs
-   [ ] TipTop Score

## Above and beyond

-   [ ] Include Better Predict TMHs Than Transmem
  -   [ ] Inside outside topology
  -   [ ] Remove signal peptides
-   [ ] Brew/Apt/Bash installer script?
-   [ ] Implement in Django



<!--
TMs In Protein TOPology = TiPTop
 -->

 # Back-end sources and references

 `topdb_all.xml` from http://topdb.enzim.hu Dobson, L., Langó, T., Reményi, I. & Tusnády, G. E. Expediting topology data gathering for the TOPDB database. Nucleic Acids Res. 43, D283–D289 (2015).

 `mptopoTblXml.xml` from http://blanco.biomol.uci.edu/mptopo/ Jayasinghe, S. MPtopo: A database of membrane protein topology. Protein Sci. 10, 455–458 (2001).

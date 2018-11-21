# Tip Top Table

<!--
## A tool to evaluate the topological preference of a TMH based on a population of TMHs with known topology

This runs the TMH sequence through a positionally dependent matrix of residue scores and checks the total score between forwards/backwards runs of the TMH. A greater difference indicates a greater topological preference. The advantage of this method is that the individual contribution of each residue are calculated, and whilst the accuracy of the predictor may not always be the highest overall, it allows for sensitive evaluation of the topology of a TMP without the need for hidden layers in neural networks or HMMs.

Enter your input below. Note that this is not the full protein sequence, nor a fasta formatted sequence. The sequence should be the predicted, or experimentally derived TMH with ±5 flanking residues.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/JamesABaker/TMH-open-topology/blob/master/)
-->

## Input

UniProt IDs in a list separated by a new line.

```
V4JVY4
V4MIU5
V4LN14
A4S3Y1
V4LFC3
A4S6G3
V4LV30
A4SB49
V4L9D6
A4SA02
```

## Output

A comma separated table of transmembrane helix positions, topology, evidence source, info on TMH ...

# Goals

## Minimum

-   [x] Verify if in MPTOPO
      -  [x] TMH positions
      -  [ ] I/O topology
-   [x] Verify if in TOPDB
      -  [x] TMH positions
      -  [ ] I/O topology
-   [x] Verify if in UniProt
      -  [x] TMH positions
      -  [ ] I/O topology
-   [x] Record source evidence

## Useful
-   [ ] Membrane location
-   [ ] Pore residue score
-   [ ] Delta H for each TMHs
-   [ ] TipTop Score
-   [ ] Hydrophobicity
-   [ ] Include Better Predict TMHs Than UniProt's Transmem

## Above and beyond


-   [ ] Remove signal peptides
-   [ ] Implement in Django
-   [ ] Beta barrels

<!--

TMs In Protein TOPology = TiPTop
TYpical Protein TOPology = TypTop
TypIcal Protein TOPology = TipTop
Tip Top Protein Topology Table

 -->

# Back-end sources and references

UniProt `uniprot_bin/*.txt` <https://www.uniprot.org/> Bateman, A. et al. UniProt: the universal protein knowledgebase. Nucleic Acids Res. 45, D158–D169 (2017).

 TopDB `topdb_all.xml` from <http://topdb.enzim.hu> Dobson, L., Langó, T., Reményi, I. & Tusnády, G. E. Expediting topology data gathering for the TOPDB database. Nucleic Acids Res. 43, D283–D289 (2015).

 MPTOPO `mptopoTblXml.xml` from <http://blanco.biomol.uci.edu/mptopo/> Jayasinghe, S. MPtopo: A database of membrane protein topology. Protein Sci. 10, 455–458 (2001).

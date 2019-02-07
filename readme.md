# Tip Top Table

<!--
## A tool to evaluate the topological preference of a TMH based on a population of TMHs with known topology

This runs the TMH sequence through a positionally dependent matrix of residue scores and checks the total score between forwards/backwards runs of the TMH. A greater difference indicates a greater topological preference. The advantage of this method is that the individual contribution of each residue are calculated, and whilst the accuracy of the predictor may not always be the highest overall, it allows for sensitive evaluation of the topology of a TMP without the need for hidden layers in neural networks or HMMs.

Enter your input below. Note that this is not the full protein sequence, nor a fasta formatted sequence. The sequence should be the predicted, or experimentally derived TMH with ±5 flanking residues.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/JamesABaker/TMH-open-topology/blob/master/)
-->

## Input



UniProt IDs in a list separated by a new line.

Tricky TMHs for testing:
P32897

 This should be fetched from:
 Human SwissProt query https://www.uniprot.org/uniprot/?query=reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22&format=fasta&force=true&sort=score

## Output

A comma database of transmembrane helix positions, topology, evidence source, and computed information about a TMH.

# Goals

## Minimum

-   [x] Verify if in MPTOPO
      -  [x] TMH positions
      -  [x] I/O topology
-   [x] Verify if in TOPDB
      -  [x] TMH positions
      -  [x] I/O topology
-   [x] Verify if in UniProt
      -  [x] TMH positions
      -  [x] I/O topology
-   [x] Record source evidence
-   [ ] Database format
-   [ ] Automate database building

## Useful

-   [x] Membrane location
-   [ ] Pore residue score
-   [ ] Delta H for each TMHs
-   [ ] TipTop Score
-   [ ] Hydrophobicity
-   [ ] Include Better Predict TMHs Than UniProt's Transmem

## Above and beyond

-   [x] Remove signal peptides
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

 OPM `opm/*.txt` from <https://opm.phar.umich.edu/> Lomize, M. A., Pogozheva, I. D., Joo, H., Mosberg, H. I. & Lomize, A. L. OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res. 40, D370–D376 (2012).

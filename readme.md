[![CodeFactor](https://www.codefactor.io/repository/github/jamesabaker/tip-top-table/badge)](https://www.codefactor.io/repository/github/jamesabaker/tip-top-table)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Build Status](https://travis-ci.com/JamesABaker/Tip-Top-Table.svg?token=Ge7HGfcCxBUTsSpxLzMX&branch=master)](https://travis-ci.com/JamesABaker/Tip-Top-Table.svg?token=Ge7HGfcCxBUTsSpxLzMX&branch=master)

# Database of Variants in TMPs.

This is a tool to help to evaluate disease variants in the context of transmembrane proteins.
It curates TMHs boundaries and topology from several sources and cross references those against several variant databases.

<!--
TMs In Protein TOPology = TiPTop
TYpical Protein TOPology = TypTop
TypIcal Protein TOPology = TipTop
Tip Top Protein Topology Table
VarTMPs
 -->

# Back-end sources and references

UniProt `uniprot_bin/*.txt` <https://www.uniprot.org/> Bateman, A. et al. UniProt: the universal protein knowledgebase. Nucleic Acids Res. 45, D158–D169 (2017).

TopDB `topdb_all.xml` from <http://topdb.enzim.hu> Dobson, L., Langó, T., Reményi, I. & Tusnády, G. E. Expediting topology data gathering for the TOPDB database. Nucleic Acids Res. 43, D283–D289 (2015).

MPTOPO `mptopoTblXml.xml` from <http://blanco.biomol.uci.edu/mptopo/> Jayasinghe, S. MPtopo: A database of membrane protein topology. Protein Sci. 10, 455–458 (2001).

OPM `opm/*.txt` from <https://opm.phar.umich.edu/> Lomize, M. A., Pogozheva, I. D., Joo, H., Mosberg, H. I. & Lomize, A. L. OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res. 40, D370–D376 (2012).

TMSOC from <http://tmsoc.bii.a-star.edu.sg/> Wong, W.-C., Maurer-Stroh, S., Schneider, G. & Eisenhaber, F. Transmembrane helix: simple or complex. Nucleic Acids Res. 40, W370–W375 (2012).

CD-HIT from <https://github.com/weizhongli/cdhit> Huang, Y., Niu, B., Gao, Y., Fu, L. & Li, W. CD-HIT Suite: A web server for clustering and comparing biological sequences. Bioinformatics 26, 680–682 (2010).

FunFam

ΔG TM insertion from <https://github.com/ElofssonLab/dgpred> Hessa, T. et al. Molecular code for transmembrane-helix recognition by the Sec61 translocon. Nature 450, 1026–1030 (2007).

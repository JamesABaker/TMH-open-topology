import Bio
from Bio import SeqIO
from Bio import SwissProt

filename = "scripts/external_datasets/uniprot_bin/Q92956.txt"

for record in SeqIO.parse(filename, "swiss"):
      for i, f in enumerate(record.features):
          if f.type == "TOPO_DOM":
              print(f.qualifiers["description"])

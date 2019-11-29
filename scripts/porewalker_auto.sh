#! /bin/sh
# Example Usage:
# sh porewalker_auto.sh <4 letter id code>
input="vartmh.txt"
#id=$1 # This should be a 4 letter PDB id
while IFS= read id
do
  dir=/nfs/nobackup/thornton/pdbsum/PoreWalker/vartmh$id
  mkdir $dir
  chmod 777 $dir
  pathfilename=$dir/vartmh$id.pdb #The pdb file must match the directory pathfilename

  # Fetches the pdb structure from the PDBe in its biological unit
  wget -O $pathfilename.gz http://www.ebi.ac.uk/pdbe/static/entry/download/$id-assembly-1.cif.gz
  # Fetches rcsb structural unit
  # wget -O $pathfilename.gz http://www.rcsb.org/pdb/files/$id.pdb.gz
  gunzip -d $pathfilename.gz
  python cif2pdb.py $id-assembly-1.cif $id.pdb
  perl -w /nfs/research1/thornton/www/software/cgi-bin/data/PoreWalker/pdb_format.pl $pathfilename > $dir/temp.pdb # Cleans the PDB file for porewalker
  mv $dir/temp.pdb $pathfilename



  # gubbins required for notification of completion
  # Don't include .pdb in the other file names
  echo vartmh$id > $dir/vartmh$id-id-list.el
  echo "bakerjames@ebi.ac.uk" > $dir/vartmh$id-mail.el
  messfile=$dir/vartmh$id-message.el
  echo "Calculation is ready, please visit the following link:" > $messfile
  echo "<P>http://www.ebi.ac.uk/thornton-srv/software/cgi-bin/data/PoreWalker/start.pl?vartmh$id" >> $messfile

  #ensure that it is chmod 777 again, otherwise the scripts can delete it, or do anything with it (as they are running under romans user id).
  chmod -R 777 $dir/*
  echo vartmh$id >> /nfs/nobackup/thornton/pdbsum/PoreWalker/latest.lst
done <"$input"

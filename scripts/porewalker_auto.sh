#! /bin/sh
# Example Usage:
# sh porewalker_auto.sh <4 letter id code>

id=$1 # This should be a 4 letter PDB id
dir=batch$id
mkdir $dir
chmod 777 $dir
pathfilename=$dir/batch$id.pdb #The pdb file must match the directory pathfilename

# Fetches the pdb structure from the PDB
wget -O $pathfilename.gz http://www.rcsb.org/pdb/files/$id.pdb.gz
gzip -d $pathfilename.gz
perl -w /nfs/research1/thornton/www/software/cgi-bin/data/PoreWalker/pdb_format.pl $pathfilename > $dir/temp.pdb # Cleans the PDB file for porewalker
mv $dir/temp.pdb $pathfilename



# gubbins required for notification of completion
# Don't include .pdb in the other file names
echo batch$id > $dir/batch$id-id-list.el
echo "bakerjames@ebi.ac.uk" > $dir/batch$id-mail.el
messfile=$dir/batch$id-message.el
echo "Calculation is ready, please visit the following link:" > $messfile
echo "<P>http://www.ebi.ac.uk/thornton-srv/software/cgi-bin/data/PoreWalker/start.pl?batch$id" >> $messfile

#ensure that it is chmod 777 again, otherwise the scripts can delete it, or do anything with it (as they are running under romans user id).
chmod -R 777 $dir/*
echo batch$id >> latest.lst

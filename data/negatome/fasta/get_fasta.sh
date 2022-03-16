WEBADDR=https://www.uniprot.org/uniprot
for i in {1..1733}
do
  ID=$(sed -n $i'p' ./ids.txt)
  echo $ID
  wget $WEBADDR/$ID'.fasta'
done

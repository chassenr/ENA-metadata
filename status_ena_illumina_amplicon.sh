# run on server:
mkdir ENA_xref
cd ENA_xref

# check with primary study accession
mkdir xref_primary_study_accnos
cd xref_primary_study_accnos
while read line
do
  curl -s -o ${line}"_xref.txt" "https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession="${line}"&source=EuropePMC"
done < ../ena_primary_study_accnos.txt

# filter studies with no hits (file size is 17)
ls -l *xref.txt | sed 's/ \+/\t/g' | cut -f5 | sort | uniq -c
find . -maxdepth 1 -size -20c | grep "xref" > ../no_hits_primary.txt
# 5065 no hits
find . -maxdepth 1 -size -20c | grep "xref" | xargs rm
cd ..

# check with secondary study accession
mkdir xref_secondary_study_accnos
cd xref_secondary_study_accnos
while read line
do
  curl -s -o ${line}"_xref.txt" "https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession="${line}"&source=EuropePMC"
done < ../ena_secondary_study_accnos.txt

# filter studies with no hits (file size is 17)
find . -maxdepth 1 -size -20c | grep "xref" > ../no_hits_secondary.txt
# 4729 no hits
find . -maxdepth 1 -size -20c | grep "xref" | xargs rm
cd ..



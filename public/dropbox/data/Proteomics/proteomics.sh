#!/bin/sh

rm -f proteins_unsorted.txt

# Get full list of proteins
for f in RPPA/*; do
	tail -n +3 $f | cut -f1 >> proteins_unsorted.txt
done

# Sort full list of proteins
sort proteins_unsorted.txt | uniq > proteins.txt

echo "MappedFileName" > proteomics.txt

for f in RPPA/*; do
	# Append patient
	head -n 1 $f | cut -f 2 >> proteomics.txt
	
	# Sort protein file
	tail -n +3 $f | sort > $f.sorted
	
	# Left join with the list of all proteins
	join -a 1 -t $'\t' -e '0' -o auto proteins.txt $f.sorted > tmp.txt
	cat tmp.txt > proteins.txt
	
	# Cleanup
	rm -f $f.sorted
done

# Transpose the proteins
while read protein; do
	echo $protein | tr " " "\n" > protein.txt
	paste proteomics.txt protein.txt > tmp.txt
	cat tmp.txt > proteomics.txt
done < proteins.txt

# Get Clinical patient to MappedFileName mapping
cut -f 1,3 ProteomicsMapping.txt | awk '{print $2 "\t" $1}' > tmp.txt

# Get header
head -n 1 tmp.txt > clinical.txt

# Sort Clinical patients
sort tmp.txt >> clinical.txt

# Join patients with expression values, and remove MappedFileName
join --header -t $'\t' clinical.txt proteomics.txt | cut -f 2- | sort > tmp.txt
cat tmp.txt > proteomics.txt

# Cleanup
rm -f proteins_unsorted.txt tmp.txt protein.txt proteins.txt clinical.txt

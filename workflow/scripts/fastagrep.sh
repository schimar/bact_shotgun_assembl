for file in contigs.fasta;do
	grep -F “>” $file | sed -e ‘s/_/ /g’ |sort -nrk 6 |awk ‘$6>=4.0 && $4>=500 {print $0}’| sed -e ‘s/ /_/g’|sed -e ‘s/>//g’>$file.txt
	echo sequences to keep
	wc -l $file.txt
	echo running fastagrep.pl
	../../../scripts/fastagrep.pl -f $file.txt $file > HCov.$file
	echo sequences kept
	grep -c “>” HCov.$file
done

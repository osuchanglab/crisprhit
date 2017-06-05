#!/bin/bash
crisprhit="../crisprhit.py"
bout="spacer_v_protospacer.tab"
chout="benchmark.crisprhit"
opt=$1
if [[ ! -s "protospacers_rev.fas.nin" ]]; then
	makeblastdb -dbtype nucl -in protospacers_rev.fas
fi
if [[ ! -s "$bout" || $opt == "blast" ]]; then
	blast="blastn -db protospacers_rev.fas -query spacer.fna -outfmt 7 -max_target_seqs 50000 -evalue 1 -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -word_size 5 -out $bout"
	echo Running BLASTN command: $blast
	$blast
	echo Finished running BLASTN
fi
if [[ ! -s "$chout" || $opt == "crisprhit" ]]; then
	ch="$crisprhit --verbose --PAM CTT protospacers_rev.fas spacer.fna $bout"
	echo "Running crisprhit command: $ch > $chout"
	$ch > $chout
	echo Finished running crisprhit
fi

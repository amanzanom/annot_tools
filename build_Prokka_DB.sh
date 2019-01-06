#!/bin/bash

INGBK=$1;
OUTPREFIX=$2;

perl ~/software/ANNOTATION/OTHER/prokka/bin/prokka-genbank_to_fasta_db \
	--verbose \
	--format genbank \
	--idtag 'locus_tag' \
	--pseudo \
	--hypo \
	--gcode 11 \
	$INGBK 1> $OUTPREFIX.faa 2> log_prokkaDB;

~/software/ANNOTATION/OTHER/prokka/binaries/linux/makeblastdb \
	-dbtype prot \
	-in $OUTPREFIX.faa \
	-out $OUTPREFIX;

mv $OUTPREFIX.p* ~/software/ANNOTATION/OTHER/prokka/db/genus;


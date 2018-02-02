# crisprhit


# Installation

No installation required. After cloning the git repository, move the script to a directory in your PATH or add the cloned directory to your PATH. Install BioPython and its dependencies and you will be good to go.

## Dependencies

Required

* [python v. 2.7.8+](https://www.python.org)
* [Biopython v. 1.65+](http://biopython.org/DIST/docs/install/Installation.html#sec12)

## Rules

Rules for spacer-protospacer interactions were adapted from [Fineran et al. 2014](http://www.pnas.org/content/111/16/E1629.long).

# Usage

The purpose of this script is to identify the outcomes of spacer-protospacer interactions based on a BLAST search of the spacers to a target.

The basic syntax is:

`crisprhit.py target.fasta spacers.fasta spacer_blast_to_target.tab`

The basic criteria for a BLASTN search that is compatible with the short sequence of a CRISPR spacer is:

`blastn -evalue 1 -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -word_size 5`

## Example workflow

1. Extract spacer sequences from CRISPR locus in bacterium of interest.
    * [CRISPR-Finder](http://crispr.i2bc.paris-saclay.fr/Server/)
    * [CRISPRDetect](http://brownlabtools.otago.ac.nz/CRISPRDetect/predict_crispr_array.html)
2. Generate blastDB for target sequence using makeblastdb.
3. Run blastn search with predicted spacers against target.
    * Use settings `blastn -evalue 1 -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -word_size 5`
    * Use blast output type 6 or 7 (tab-delimited).
4. Use spacer file, target sequence file, and blast output as input to crisprhit.py.

# Output types

**FASTA** output type results in protospacer sequences printed in FASTA format. Headers include the target accession and position, spacer name, as well as metadata about the hits, mismatches, and predicted quality of the spacer-protospacer interaction.

For example:

```
 Accesion    Positions  Spacer_name     Metadata about hits and mismatches                                  Quality
>KX911187.1:c9266..9298_FH100_spacer139 hits:21 mismatch:14 PAM_hits:1 seed_hits:3 priming_mm:4 stable_mm:4 quality:other
CCCACCTCGTCGGAGGCGGCATCTTTGCTGGCAtg
Protospacer Sequence             ^^ PAM sequence
```

**Table** output is the default output type. Table output includes tab-delimited information about the Spacer-Protospacer matches.

For example:

```
#name	id	proto_seq	spacer_seq	spacer_revcom	PAM	start	end	strand	spacer	hits	misses	PAM_hits	seed_hits	priming_mm	stable_mm	hq_mm	quality	guide
ps00022	KX911187.1	CATCTACTCGGATTTTGATGAAGCTCAATGGCA	TCGTCCTTTTCTTCATCAAAAAGCGCCTCGATG	CATCGAGGCGCTTTTTGATGAAGAAAAGGACGA	tg	2043	2011	minus	FH100_spacer17	21	14	1	2	5	1	1	priming	**STTTTQST P I     I    XI PP PI   
ps00023	KX911187.1	GTACGAGAAAATCCTTCAGGAAGAAAAGGACAT	TCGTCCTTTTCTTCATCAAAAAGCGCCTCGATG	CATCGAGGCGCTTTTTGATGAAGAAAAGGACGA	ca	8889	8857	minus	FH100_spacer17	20	15	0	5	5	4	1	other	**TTSSSISS   I  X PI PP XQXP   I  X
ps00024	KX911187.1	CGCCTGGGGGTTCTCTCGTTGTGTGAGGATATTT	AAATATCCTCACACAACGAGAGAACCCCCAGGCG	CGCCTGGGGGTTCTCTCGTTGTGTGAGGATATTT	tt	14121	14088	minus	FH100_spacer18	36	0	2	7	0	0	0	perfect	**SSSSSISS   I     I     I     I    
ps00025	KX911187.1	TCATGAGATGGAGAGGCGACCGCGGCAAGCGAA	TTCGTTTGCCGCGGTCGCCTCTCCATCTCATGA	TCATGAGATGGAGAGGCGACCGCGGCAAACGAA	tt	5835	5803	minus	FH100_spacer20	34	1	2	6	0	0	0	hq	**SSSSTISS   I     I     I     I   

```

The majority of the columns are self-explanatory.

| Column       | Meaning                    |
| ------------ | -------------------------- |
|name          | Arbitrary protospacer name |
|id            | Accession of target sequence |
|proto_seq     | Protospacer DNA sequence |
|spacer_seq    | Spacer (query) DNA sequence |
|spacer_revcom | Spacer reverse complement sequence. Used for direct comparisons to protospacer sequence. |
|PAM           | Protospacer adjacent motif (PAM) sequence. |
|start         | Position of start of match in target sequence. |
|end           | Position of end of match in target sequence. |
|strand        | Strand for protospacer sequence. |
|spacer        | Spacer (query) name. |
|hits          | Number of matched bases between spacer and protospacer (i.e. perfect basepairing). |
|misses        | Number of mis-matched bases between spacer and protospacer. |
|PAM_hits      | Number of positions of PAM sequence that are correct. |
|seed_hits     | Number of matched bases in the seed positions. |
|priming_mm    | Number of mis-matches that are likely to produce priming events. |
|stable_mm     | Number of mis-matches that are likely to produce stable events. |
|hq_mm         | Number of mis-matches that are unlikely to produce priming/stable events. |
|quality       | Predicted quality of spacer-protospacer interaction. |
|guide         | Guide showing where the PAM, Seed, and mismatches are in the Spacer-Protospacer interaction. |

Valid qualities are perfect, hq, hq-priming, priming, stable, and other.

Perfect and hq are predicted to result in direct interference between the spacer targeted protospacer sequence.

Priming are predicted to result in acquisition of new spacer sequences by the CRISPR locus.

Spacers with stable quality are predicted to have no interaction with the protospacer sequence, leaving the target in-tact.

Spacers with other quality do not have enough signal in the spacer-protospacer interaction to accurately predict the type of outcome. These typically have more mismatches with even amounts of priming and stable mismatches.

The guide can be used to visually inspect the spacer-protospacer interactions. For example, for protospacer ps00022 above:

```
ps00022	
spacer	  TCGTCCTTTTCTTCATCAAAAAGCGCCTCGATG
guide	**STTTTQST P I     I    XI PP PI   

Meaning of Symbols:
* - PAM
S - Seed positions

Mismatches:
T - Mismatches in seed positions
I - Positions that allow for interference, even with mismatches. Interference positions are shown in the guide even when no mismatches are present.
P - Mismatch that encourages priming events. (mismatches of Cytosine)
X - Mismatch that encourages stable events.  (mismatches of Guanine)

Combination Mismatches:
Y - Mismatch that encourages stable events in the interference positions.
Q - Mismatch that encourages priming events in the interference positions.
```

**Basepair** output is similar to 'normal' BLAST output in that it shows a visual alignment of the spacer-protospacers with the guide in-line with the text. The alignment shows the spacer sequence aligned with the reverse (not reverse complement) of the protospacer sequence, to simulate the basepairing that occurs during CRISPR-dependent immunity.

Here is the same alignment as above (ps00022) in basepair output format:

```
name     Accession  Positions Strand Spacer_Name   Quality
#ps00022 KX911187.1 2043:2011 minus FH100_spacer17 priming
          **STTTTQST P I     I    XI PP PI   
            TCGTCCTTTTCTTCATCAAAAAGCGCCTCGATG
            |     |   |||||||||||  ||  | ||||
2043      gtACGGTAACTCGAAGTAGTTTTAGGCTCATCTAC           2011
```

# Options

| Option       | Values [default]         | Function |
| ------------ | ------------------------ | -------- |
|-h, --help    |                          | show this help message and exit |
|--filetype    | FASTA, GENBANK, [auto]   | Input file type provided for target sequence. |
|--PAM         | string [TT]              | Putative PAM sequence. |
|--outtype     | [all], PAM               | Output sequence type.  Valid for FASTA output only |
|--outfmt      | fasta, [table], basepair | Output file format. |
|--length      | full, seed, [fill]       | Matching protospacer lengths. |
|--plength     | integer [2]              | PAM length to include. |
|--width       | integer [60]             | Output column width. Valid for basepair output file format. |
|--match       | True, [False], Partial   | Force PAM match. |
|--spacers     | [all], perfect, hq, priming, partial, other, stable | Limit output to a particular spacer category. |
|--filters     | [all], pre, first, second, third, post              | Change spacer filters to control quality. Early filters are more confident (pre, first, second) than later filters. |
|--mmlimit     | integer [15]             | Set mismatch limit. |
|--prefix      | string [ps]              | Prefix for arbitrary spacer names. |
|--hits        | [all], top               | Amount of hits to report for each spacer. top restricts to the single best hit against the target. |
|-v, --verbose |                          | Print progress messages. |
|-q, --quiet   |                          | Hide warning messages. |
|--debug       |                          | Print debugging messages. |
|-V, --version |                          | Print version message and quit. |

# History

v1.0.0 - 2018-02-01 - Developer version merged to stable version. Many improvements to the quality of the algorithm.

v0.9.0 - 2017-06-01 - First revision released to GitHub.

# Credits

* [Edward Davis](mailto:davised.dev@gmail.com)
* [Dr. Jeff Chang](mailto:changj@science.oregonstate.edu)

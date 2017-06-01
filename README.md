# crisprhit


# Installation

No installation required. After cloning the git repository, move the script to a directory in your PATH or add the cloned directory to your PATH. Install BioPython and its dependencies and you will be good to go.

## Dependencies

Required

* [python v. 2.7.8+](https://www.python.org)
* [Biopython v. 1.65+](http://biopython.org/DIST/docs/install/Installation.html#sec12)

# Usage

The purpose of this script is to identify the outcomes of spacer-protospacer interactions based on a BLAST search of the spacers to a target.

The basic syntax is:

`crisprhit.py target.fasta spacers.fasta spacer_blast_to_target.tab`

The basic criteria for a BLASTN search that is compatible with the short sequence of a CRISPR spacer is:

`blastn -evalue 1 -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -word_size 5`

Tab-delimited output is by default.

# History

v0.9.0 - 2017-06-01 - First revision released to GitHub.

# Credits

* [Edward Davis](mailto:davised.dev@gmail.com)
* [Dr. Jeff Chang](mailto:changj@science.oregonstate.edu)

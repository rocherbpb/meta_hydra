## Species discovery using COI barcoding

## UNDER CONSTRUCTION!

### Specimen identification versus species discovery
Specimen identification involves matching the specimen or query sequence to sequences in a reference library.
Species discovery involves characterising barcode diversity, delimiting this diversity into species and, effectively, the creation of a reference library. The former, ultimately, on the latter.  

### COI Barcode data
Raw data will typically consist of 658bp sequences from the mtDNA COI locus acquired via [Folmer et al (1994)](https://pubmed.ncbi.nlm.nih.gov/7881515) primers. 

Much of the analyses described below relies on the contruction of a sequence alignment. This can be done using the [Muscle](https://academic.oup.com/nar/article/32/5/1792/2380623) algorithm implemented in [SeaView](https://academic.oup.com/mbe/article/27/2/221/970247) (available [here](http://doua.prabi.fr/software/seaview)). Amino acid and codon positions can be determined by [TranslatorX](https://academic.oup.com/nar/article/38/suppl_2/W7/1094709) (available [here](http://translatorx.co.uk)) by running the alignment with "Guess most likely reading frame", "Invertebrate mitochondrial" genetic code, and the Muscle algorithm setting. Alignments are also performed natively in several analytical tools found in the [BOLD workbench](http://www.boldsystems.org).

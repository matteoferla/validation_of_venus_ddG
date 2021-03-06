# Benchmarks of the ∆∆G calculations in Venus

[Venus](https://venus.sgc.ox.ac.uk/) does not do a full minimisation but cheats by only minimising the neighbourhood 
before and after introducing the variant.
Therefore, it is important to assess the accuracy of Venus ddG calculations.

* For experimental details see [details](details.md)
* For the refactored code for scoring see [scoring notes](scoring.md)
* For the refactored code for analysis see [analyses notes](analyses.md)
* For the dataset of results see XXX.

For previous iterations and other files, see:

* For the code used for the O2567 dataset see [scoring](notes/O2567_scoring.md).
* For the code used for the ProTherm-star dataset see [scoring](notes/ProTherm-star_scoring.md).
* For past summary see [preliminary results](notes/preliminary_results.md)
* Other files within [notes folder](notes).

The dataset run were:

* O2567
* ProTherm-star
* S1342

These are subsets of the ProTherm database and are aimed at being more balanced in terms of various factors.
These factors differ, so refer to the original papers for details.
For [example](notes/residue_distro.md), O2567 is balanced in the number of mutants of a given residues mutating to another, 
whereas S1342 is dominated by alanine scanning mutants. Etc. etc. etc.


Do note that the code presented in the markdown was originally on Jupyter notebooks.
So there may be leftover lines that should not worry the reader (see [note on extra fluff](notes/extra_fluff.md))
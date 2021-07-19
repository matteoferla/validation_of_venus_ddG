## Scorefunction choice

In addition to the choice of cartesian or not, there is also the weights themselves.

Whereas most of the modern litterature uses `ref2015`, this does not mean it is the most up to date.

Using my helper module (`pip3 install pyrosetta-help`) to inspect the scorefunctions:

```python
from pyrosetta_help import WeightWatcher
ww = WeightWatcher()
for sfn in sorted(ww.possible_scorefxn_names):
    comments = ww.get_scorefxn_comments(sfn).replace('#', '').replace('\n', ' - ')
    print(f"* `{sfn}` {comments}")
```
I get the following:

* `HBNet` 
* `PPI_discrimination`  These weights are for discriminating native complexes from decoys -  for small-molecule inhibitors of protein-protein interactions. -  These weights were developed as described in the publication: -  "Enhancements to the Rosetta energy function enable improved identfication -     of small-molecules that inhibit protein-protein interactions" -  Kelow, Bazzoli, and Karanicolas (submitted) -  The reference weights here are not meaningful. -  These weights are for scoring only, NOT anything where the structure moves.
* `ProQ2` 
* `ProQM` 
* `RS_centroid`  score 6 of rosetta++
* `abinitio_remodel_cen` 
* `allfold_soft_rep` 
* `antibody_design_talaris2013` Weights for RosettaAntibodyDesign - Jared Adolf-Bryfogle June 2014 - Dunbrak Lab - Weight set based on talaris2013 with dihedral, coordinate, and atom_pair constraints - The atom_pair constraints are there for SiteConstraints between the epitope and paratope. - They have been optimized (by eye) to not dominate the scorefunction as there are many. This is still not optimum, and may have to be done dynamically based on the number of constraints.    - Reference energies will be reweighted to capitulate native CDR sequence recoveries - I also have chainbreak and linear chainbreak at 100 as I hate chainbreaks.
* `beta`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 - METHOD_WEIGHTS ref 1.8394 3.6196 -2.3716 -2.7348 1.0402 0.83697 -0.45461 0.73287 -1.5107 0.18072 0.60916 -0.93687 -2.4119 -0.18838 -1.2888 -0.77834 -1.0874 1.9342 1.6906 1.2797
* `beta_cart`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 - METHOD_WEIGHTS ref 1.8394 3.6196 -2.3716 -2.7348 1.0402 0.83697 -0.45461 0.73287 -1.5107 0.18072 0.60916 -0.93687 -2.4119 -0.18838 -1.2888 -0.77834 -1.0874 1.9342 1.6906 1.2797
* `beta_cst`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 - METHOD_WEIGHTS ref 1.8394 3.6196 -2.3716 -2.7348 1.0402 0.83697 -0.45461 0.73287 -1.5107 0.18072 0.60916 -0.93687 -2.4119 -0.18838 -1.2888 -0.77834 -1.0874 1.9342 1.6906 1.2797
* `beta_genpot`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 - METHOD_WEIGHTS ref 1.8394 3.6196 -2.3716 -2.7348 1.0402 0.83697 -0.45461 0.73287 -1.5107 0.18072 0.60916 -0.93687 -2.4119 -0.18838 -1.2888 -0.77834 -1.0874 1.9342 1.6906 1.2797
* `beta_genpot_cart`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 - METHOD_WEIGHTS ref 1.8394 3.6196 -2.3716 -2.7348 1.0402 0.83697 -0.45461 0.73287 -1.5107 0.18072 0.60916 -0.93687 -2.4119 -0.18838 -1.2888 -0.77834 -1.0874 1.9342 1.6906 1.2797
* `beta_genpot_cst`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 - METHOD_WEIGHTS ref 1.8394 3.6196 -2.3716 -2.7348 1.0402 0.83697 -0.45461 0.73287 -1.5107 0.18072 0.60916 -0.93687 -2.4119 -0.18838 -1.2888 -0.77834 -1.0874 1.9342 1.6906 1.2797
* `beta_genpot_soft`  beta_nov16_soft -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 -     - this is "soft" version of beta_nov16, equivalent to "soft_rep" (for score12)
* `beta_july15`  hpark Aug 28 2015 -  There are two reference weights below; first one is optimized in the minimized context, second is in packer.  -  by default the second will be used, but you can switch to the first one if you work with minimization. - METHOD_WEIGHTS ref 1.66209 3.70349 -2.00354 -2.38137 0.65595 0.81210 -0.42905 2.27849 -0.47142 1.68043 2.59866 -1.09912 -0.97642 -1.60738 -0.76238 0.23890 0.72083 2.28760 1.86226 0.41729
* `beta_july15_cart`  hpark Aug 28 2015 -  There are two reference weights below; first one is optimized in the minimized context, second is in packer.  -  by default the second will be used, but you can switch to the first one if you work with minimization. - METHOD_WEIGHTS ref 1.66209 3.70349 -2.00354 -2.38137 0.65595 0.81210 -0.42905 2.27849 -0.47142 1.68043 2.59866 -1.09912 -0.97642 -1.60738 -0.76238 0.23890 0.72083 2.28760 1.86226 0.41729
* `beta_july15_cst`  hpark Aug 28 2015 -  There are two reference weights below; first one is optimized in the minimized context, second is in packer.  -  by default the second will be used, but you can switch to the first one if you work with minimization. - METHOD_WEIGHTS ref 1.66209 3.70349 -2.00354 -2.38137 0.65595 0.81210 -0.42905 2.27849 -0.47142 1.68043 2.59866 -1.09912 -0.97642 -1.60738 -0.76238 0.23890 0.72083 2.28760 1.86226 0.41729
* `beta_nov15`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) -  - METHOD_WEIGHTS ref 1.82468 3.75479 -2.14574 -2.72453 1.21829 0.79816 -0.30065 2.30374 -0.71458 1.66147 2.15735 -1.34026 -1.94321 -1.45095 -0.59474 -0.28969 1.15175 2.64269 2.26099 0.58223
* `beta_nov15_cart`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) -  - METHOD_WEIGHTS ref 1.82468 3.75479 -2.14574 -2.72453 1.21829 0.79816 -0.30065 2.30374 -0.71458 1.66147 2.15735 -1.34026 -1.94321 -1.45095 -0.59474 -0.28969 1.15175 2.64269 2.26099 0.58223
* `beta_nov15_cst`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) -  - METHOD_WEIGHTS ref 1.82468 3.75479 -2.14574 -2.72453 1.21829 0.79816 -0.30065 2.30374 -0.71458 1.66147 2.15735 -1.34026 -1.94321 -1.45095 -0.59474 -0.28969 1.15175 2.64269 2.26099 0.58223
* `beta_nov15_soft`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    This is "soft" version of beta_nov15, equivalent to "soft_rep" to score12, which will be useful when designing at less accurate conditions. - 
* `beta_nov16`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016
* `beta_nov16_cart`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 
* `beta_nov16_cst`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016
* `beta_nov16_soft`  beta_nov16_soft -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 -     - this is "soft" version of beta_nov16, equivalent to "soft_rep" (for score12)
* `beta_peptide` 
* `beta_peptide_mm` 
* `beta_peptide_soft_rep_design` 
* `beta_peptide_soft_rep_mm`
* `beta_soft`  beta_nov16_soft -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016 -     - this is "soft" version of beta_nov16, equivalent to "soft_rep" (for score12)
* `beta_soft_rep` 
* `cen_std` 
* `cen_std_smooth` Smooth version of Cen Std using Frank Dimaio's updated centroid statistics.  Highly recommended.
* `corrections_conway2016`  The scorefile (and corresponding rotamer library correction) is described in Conway and DiMaio, Protein Science 2016. -  NOTE! This scorefile should be used with the flag -dun10_dir rotamer/corrections_conway2016 - 
* `covalent_labeling_fa`  Ref15 weights with the addition of covalent_labeling_fa to be used when scoring with sparse covalent labeling exp data
* `ddg` 
* `ddg_monomer` 
* `design_hpatch` 
* `dna` 
* `dna_gb`  dna interface weight set with the gb_elec term and without the pair term -  trained (single pass) in mini by J. Ashworth 2008 -  ps. warning: gb_elec is very slow (expect 5-10x slowdown for packing)
* `dna_no_gb` 
* `docking` 
* `empty` 
* `enzdes` 
* `enzdes_polyA_min` 
* `fldsgn_cen` 
* `franklin2018`  franklin2018: A membrane environemt energy function with implicit lipids -  author: Rebecca F. Alford (raalford3@jhu.edu) -  cite: Alford et al. (2019) "Protein structure prediction and design in a biologically-realistic implicit membrane" - METHOD_WEIGHTS ref 1.82468 3.75479 -2.14574 -2.72453 1.21829 0.79816 -0.30065 2.30374 -0.71458 1.66147 2.15735 -1.34026 -1.94321 -1.45095 -0.59474 -0.28969 1.15175 2.64269 2.26099 0.58223
* `franklin2019`  franklin2019: A membrane environemt energy function with implicit lipids -  author: Rebecca F. Alford (raalford3@jhu.edu) -  cite: Alford et al. (2019) "Protein structure prediction and design in a biologically-realistic implicit membrane" - METHOD_WEIGHTS ref 1.82468 3.75479 -2.14574 -2.72453 1.21829 0.79816 -0.30065 2.30374 -0.71458 1.66147 2.15735 -1.34026 -1.94321 -1.45095 -0.59474 -0.28969 1.15175 2.64269 2.26099 0.58223
* `gauss`  weights trained by optE -  fixed params: -  free params: fa_atr gauss rama fa_dun p_aa_pp sa
* `gauss1` 
* `gen_born` 
* `genpot_full`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016
* `genpot_full_cart`  beta_nov16 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2016
* `hpol_only_relaxed_set_5` 
* `hpol_wcpd_relaxed_set_4` 
* `hrf_dynamics` 
* `hrf_ms_labeling` 
* `hydrate_score12` To use this energy function with the defaults for which it was optimized, the option -restore_pre_talaris_2013_behavior must be passed. -  - Reoptimized for the Hydrate/SPaDES protocol which introduces new scoring terms
* `interchain_cen` 
* `interface` 
* `latest_greatest` This is actually a Talaris13 derived scorefunction... so defo not latest
* `ligand` 
* `ligand_gen_born` 
* `ligand_soft_rep` 
* `ligand_soft_rep_gen_born` 
* `ligandprime` Original ligand weights with the reference energy optimized for sequence recovery of the binding pocket. OptE used and trained on 51 proteins. - Dataset used for optimization was the CSAR dataset http://www.csardock.org/ 2012 - Benchmark of weights on 234 proteins resulted in a 7% sequence recovery
* `loop_fragsample`  loop_fragsample  used in SlidingWindowLoopClosure - 
* `mainchain_potential_generation`  A scorefunction used by the make_mainchain_potential application. -  This is essentially mm_std_fa_elec_dslf_fa13.wts with the intra-residue -  electrostatics and intra_residue hydrogen bonding turned on, and no -  unfolded or constriant terms. -  Added by Vikram K. Mulligan (vmulligan@flatironinstitute.org) on -  30 November 2018. 
* `make_rot_lib` 
* `make_rot_lib_orig`
* `make_rot_lib_orig_bl_ba` 
* `membrane_highres`  @file: membrane_highres -  @desc: All atom energy function for membrane proteins: combines Rosetta's -  all atom energy function with membrane environment and solvation from the -  lazaridis-karplus implicit membrane model. -  @note: deprecated - use mpframework instead -  @cite: Barth, 2007, PNAS
* `membrane_highres_Menv_smooth`  @file: membrane_highres_Menv_smooth -  @desc: All atom energy function for membrane proteins: combines Rosetta's -  all atom energy function with membrane environment and solvation from the -  lazaridis-karplus implicit membrane model. -  @note: continuous for minimization, deprecated - use mpframework instead -  @cite: Barth, 2007, PNAS
* `mm_std` 
* `mm_std_cart` 
* `mm_std_cart_no_hbond` 
* `mm_std_fa_elec_dslf_fa13` Test addition of fa_elec and dslf_fa13 to mm_std - fa_elec weight approximately proportional to the talaris2013 weight relative to the hbond term weights.
* `mm_std_fa_elec_dslf_fa13_cart` Test addition of fa_elec and dslf_fa13 to mm_std - fa_elec weight approximately proportional to the talaris2013 weight relative to the hbond term weights.
* `mm_std_fa_elec_dslf_fa13_split_unfolded` Test addition of fa_elec and dslf_fa13 and split unfolded terms to mm_std - fa_elec weight approximately proportional to the talaris2013 weight relative to the hbond term weights.
* `mm_std_fa_elec_dslf_fa13_split_unfolded_cart` Test addition of fa_elec and dslf_fa13 and split unfolded terms to mm_std - fa_elec weight approximately proportional to the talaris2013 weight relative to the hbond term weights.
* `mm_std_no_hbond` 
* `motif_dock_score` 
* `mpdocking_cen_14-7-23_no-penalties` for mpdocking for MPframework paper - membrane - docking
* `mpdocking_cen_14-7-25_low-penalties` test for mpdocking - membrane - docking
* `mpdocking_fa_14-7-23_no-penalties` for mpdocking for MPframework paper - membrane - docking
* `mpdocking_fa_14-7-25_low-penalties` test for mpdocking - membrane - docking
* `mpframework_cen_2006`  file: mpframework_cen_2014.wts -  @desc: Scores low resolution pairwise, packing density, and environment energies -  in a 5-layer membrane environment. Also includes penalties for structures unfavorably -  embedded in the membrane. -  @note: supported by the membrane framework -  @cite: [Yarov-Yaravoy, 2006, Proteins] - Mlipo 1.0
* `mpframework_docking_cen_2015` for mpdocking for MPframework paper - membrane - docking
* `mpframework_docking_fa_2015` for mpdocking for MPframework paper - membrane - docking
* `mpframework_fa_2007`  file: mpframework_fa_2014.wts -  @desc: All atom energy function for membrane proteins: combines Rosetta's -  all atom energy function with membrane environment and solvation from the -  lazaridis-karplus implicit membrane model. -  @note: supported by the membrane framework - but use mpframework_smooth_fa instead -  @cite: Barth, 2007, PNAS
* `mpframework_pHmode_fa_2015`  file: mpframework_pHmode_fa_2014.wts -  @desc: Combines the smoothed membrane all atom energy function with -  pH mode - enables the packer to sample protonated rotamers where -  applicable at a given pH (specified via command line) -  @note: supported by the membrane framework -  @cite: Terms from [Barth, 2007, PNAS] and [Kilambi, 2013, Proteins]
* `mpframework_smooth_fa_2012`  file: mpframework_smooth_fa_2014.wts -  @desc: All atom energy function for membrane proteins: combines Rosetta's -  all atom energy function with membrane environment and solvation from the -  lazaridis-karplus implicit membrane model. -  @note: supported by the membrane framework - continuous for minimization -  @cite: Barth, 2007, PNAS
* `mpframework_smooth_fa_2012_cart`  file: mpframework_smooth_fa_2014.wts -  @desc: All atom energy function for membrane proteins: combines Rosetta's -  all atom energy function with membrane environment and solvation from the -  lazaridis-karplus implicit membrane model. -  @note: supported by the membrane framework - continuous for minimization -  @cite: Barth, 2007, PNAS
* `mpframework_symdock_cen_2015`  file: mpframework_symdock_cen_2014.wts -  @desc: Low resolution membrane scoring function with low resolution -  protein-protein docking score terms (used in regular docking protocol) -  @note: Supported by the membrane framework - 2014 -  @cite: Terms from [Yarov-Yaravoy, 2006, Proteins] and [Gray, 2003, Journal -  of Molecular Biology] - membrane - docking
* `mpframework_symdock_fa_2015`  Rosetta Membrane Framework Symmetric Protein-Protein Docking Highres weights -  Description: High resolution weights for symmetric docking for membranes - combines -  smoothed high resolution energy function with fa_elec term from highres asymmetric docking - docking
* `none`  blank weights file, for convenient scorefxn setup from command line via -set_weights.
* `nv_env_design` 
* `opte` 
* `orbitals` 
* `orbitals_dun10` 
* `orbitals_talaris2013`  The orbitals scorefunction was optimized using the defaults determined at the Talaris2013 score function conference. -  Changes include addition of hack_elec term, optimization of reference energies with optE on Jane Richardson's HiQ54 benchmark, -  with the following commands: -corrections:chemical:icoor_05_2009 -smooth_hack_elec -hackelec_min_dis 1.3 -hackelec_r_option false  -  -expand_st_chi2sampling -set_atom_properties fa_standard:ONH2:LK_DGFREE:-5.85 fa_standard:NH2O:LK_DGFREE:-7.8 fa_standard:Narg:LK_DGFREE:-10.0 fa_standard:OH:LK_DGFREE:-6.70 -mistakes:restore_pre_talaris_2013_behavior -  - 
* `orbitals_talaris2013_softrep`  The orbitals scorefunction was optimized using the defaults determined at the Talaris2013 score function conference. -  Changes include addition of hack_elec term, optimization of reference energies with optE on Jane Richardson's HiQ54 benchmark, -  with the following commands: -corrections:chemical:icoor_05_2009 -smooth_hack_elec -hackelec_min_dis 1.3 -hackelec_r_option false  -  -expand_st_chi2sampling -set_atom_properties fa_standard:ONH2:LK_DGFREE:-5.85 fa_standard:NH2O:LK_DGFREE:-7.8 fa_standard:Narg:LK_DGFREE:-10.0 fa_standard:OH:LK_DGFREE:-6.70 -mistakes:restore_pre_talaris_2013_behavior -  - 
* `pb_elec` 
* `pre_talaris_2013_standard` 
* `ref2015`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) 
* `ref2015_cart`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) 
* `ref2015_cart_cst`  ref2015_cart_cst -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) 
* `ref2015_cst`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. -       The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing)
* `ref2015_memb`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    Two sets of reference weight are provided. - The first is for use in "minimization context" (e.g., RTmin, min_pack, or sidechain relax). -       The second, and default set, is for use in "packing context" (e.g. Rotamer trials or packing) -  added adjustments to a membrane protein energy function (Jonathan Weinstein, Assaf Elazar and Sarel Fleishman, 2017)
* `ref2015_soft`  beta_nov15 -    beta energy function following parameter refitting (Frank DiMaio and Hahnbeom Park), November 2015 -  -    This is "soft" version of beta_nov15, equivalent to "soft_rep" to score12, which will be useful when designing at less accurate conditions. - 
* `remodel_cen` 
* `rna_minimize`  Rosetta Lennard-Jones - ENLARGE_H_LJ_WDEPTH  turns on stronger repulsion between hydrogens -  Hydrogen bonds and solvation - fa_sol          0.25  nonpolar only (see NO_LK_POLAR_DESOLVATION below) - NO_HB_ENV_DEP  no hydrogen-bond/geom-sol dependence on burial -- nucleic acid bonds get too weak. -  electrostatics (not in H-bonds) -  RNA torsion terms - linear_chainbreak 5.0  strong enough to force chainbreak closure - RNA_TORSION_POTENTIAL      RNA11_based_new  Fang's latest, derived from Richardson RNA11 set. - RNA_SYN_G_POTENTIAL_BONUS  -1.5             RNA11 penalized syn-G too much. - RNA_SUITENESS_BONUS        test/1z_6n_2[_bonus  Example of torsion corrections. -  intra-residue matching inter-residue - PUT_INTRA_INTO_TOTAL  applies to fa_atr, fa_rep, fa_sol, geom_sol_fast, and hbond terms. - fa_atr 1 - fa_rep 0.55 - fa_sol 1.0 - lk_ball_wtd 1.0 - fa_elec 1.0 - hbond_sr_bb 1.0 - hbond_lr_bb 1.0 - hbond_bb_sc 1.0 - hbond_sc 1.0
* `rna_res_level_energy4`  Rosetta Lennard-Jones -  H-bonds & solvation -  RNA torsion terms -  electrostatics (not in H-bonds ) -  bonuses/costs for free/instantiated moieties - intermol        1.0  should be 1.0 to maintain kT scale - loop_close      1.0  should be 1.0 to maintain kT scale - free_suite      2.0  should be less than each ref - free_2HOprime   1.0  should be 1.0 to maintain kT scale - ref             1.0  should be 1.0, apply METHOD_WEIGHTS without scaling. - other_pose      1.0  should be 1.0: contribution of 'sister' poses - linear_chainbreak 5.0  strong enough to force chainbreaks -  first 20 reference weights are protein; then four for DNA; then four for RNA [G,A,C,U] - NO_HB_ENV_DEP              no hydrogen-bond dependence on burial -- nucleic acid bonds get too weak. - RNA_TORSION_POTENTIAL      RNA11_based_new  Fang's latest, derived from Richardson RNA11 set. - RNA_SYN_G_POTENTIAL_BONUS  -1.5             RNA11 penalized syn-G too much. - RNA_SUITENESS_BONUS        test/1z_6n_2[_bonus  helps favor correct conformation for UUCG. Example of torsion corrections. - ENLARGE_H_LJ_WDEPTH   turns on stronger repulsion between hydrogens
* `rnp_ddg` 28 entries, 0's are DNA, last 4 RNA: G A C U - Terms in both protein and RNA sfxns - Terms in protein sfxn only - Terms in RNA sfxn only
* `rnp_ddg_stepwise`  rnp_ddg with extra stepwise-specific terms -  use this score function for RNA-protein ddG calculations with stepwise  28 entries, 0's are DNA, last 4 RNA: G A C U - Terms in both protein and RNA sfxns - Terms in protein sfxn only - Terms in RNA sfxn only -  stepwise stuff -  bonuses/costs for free/instantiated moieties - intermol        1.1  should be 1.0 to maintain kT scale - loop_close      1.1  should be 1.0 to maintain kT scale - free_suite      2.2  should be less than each ref - free_2HOprime   1.1  should be 1.0 to maintain kT scale - other_pose      1.1  should be 1.0: contribution of 'sister' poses
* `rnp_lores` 
* `sasa_only`  The Talaris2014 score function represents a small modification to the -  Talaris2013 score function: the weights are all scaled upwards so that -  fa_atr has a weight of 1.0; then the hydrogen bond strengths are decreased -  by 20%.  This has the effect of keeping the hbond strengths fixed while -  everything else gets stronger. -  -  The benchmarking performed for the O'Meara et al. (2014) hbond paper -  showed that weakening the hbond weights by 20% improved sequence recovery, -  rotamer recovery, and decoy discrimination.  This weight set is not -  (currently) the official gold standard weight set for Rosetta, though its -  use is encouraged.  It may be activated by using the -talaris2014 flag. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.8%) was -  chosen. - 
* `score0`  score0 from rosetta++, used in stage 1 of the -  ClassicAbinitio protocol. -  Score 0 has a vdw weight of 1, in R++, but then it divides -  the vdw score by 10 and rounds down to nearest integer. -  Mini does not round down to the nearest integer.
* `score0_memb`  score0 from rosetta++, used in stage 1 of the -  ClassicAbinitio protocol. -  Score 0 has a vdw weight of 1, in R++, but then it divides -  the vdw score by 10 and rounds down to nearest integer. -  Mini does not round down to the nearest integer. -  added adjustments to a membrane protein energy function (Jonathan Weinstein, Assaf Elazar and Sarel Fleishman, 2017)
* `score0_membrane`  score0 from rosetta++, used in stage 1 of the -  ClassicAbinitio protocol. -  Score 0 has a vdw weight of 1, in R++, but then it divides -  the vdw score by 10 and rounds down to nearest integer. -  Mini does not round down to the nearest integer.
* `score1`  score1 from rosetta++, used in stage 2 of ClassicAbinitio -  protocol from Rosetta++.
* `score12` To use this energy function with the defaults for which it was optimized, the option -restore_pre_talaris_2013_behavior must be passed.
* `score12_cart` 
* `score12_full` To use this energy function with the defaults for which it was optimized, the option -restore_pre_talaris_2013_behavior must be passed.
* `score12_justdisulfides` 
* `score12_no_hb_env_dep` 
* `score12_unfolded` 
* `score12_w_bicubic_interpolation_corrections`  For use with -    -use_bicubic_interpolation -  -  See Leaver-Fay et al. Methods in Enzymology 2013 for details - 
* `score12_w_corrections` 
* `score12minpack`  OptE 132. Minpack.  JSR54 training set. fa_max_dis 6.0. no_his_his_pairE.
* `score12prime`  Standard.wts + score12 patch weights, refit reference energies -  with optE for sequence profile recovery. -  Reference energies normalized so that they sum to one.
* `score13` 
* `score13_env_hb` 
* `score1_memb`  score1 from rosetta++, used in stage 2 of ClassicAbinitio -  protocol from Rosetta++. -  added adjustments to a membrane protein energy function (Jonathan Weinstein, Assaf Elazar and Sarel Fleishman, 2017)
* `score1_smooth`  smooth variant of score1 from rosetta++
* `score2`  score2 from rosetta++, used in stage 3 of ClassicAbinitio protocol.
* `score2_memb`  score2 from rosetta++, used in stage 3 of ClassicAbinitio protocol. -  added adjustments to a membrane protein energy function (Jonathan Weinstein, Assaf Elazar and Sarel Fleishman, 2017)
* `score2_smooth`  smooth variant of score2 from rosetta++
* `score3`  score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
* `score3_cenrot` vdw               1.6 opt weights 3.2
* `score3_memb`  score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol. -  added adjustments to a membrane protein energy function (Jonathan Weinstein, Assaf Elazar and Sarel Fleishman, 2017)
* `score3_rob`  score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
* `score3_smooth`  smooth variant of score3 from rosetta++
* `score4_cenrot_cartmin` 
* `score4_cenrot_desgin` 
* `score4_cenrot_relax` 
* `score4_cenrot_relax_cart` 
* `score4_cenrot_repack` 
* `score4_smooth` 
* `score4_smooth_cart` 
* `score4_smooth_cen_relax` 
* `score5`  score5.wts, used in stage 3 of ClassicAbinitio protocol
* `score5_memb`  score5.wts, used in stage 3 of ClassicAbinitio protocol -  added adjustments to a membrane protein energy function (Jonathan Weinstein, Assaf Elazar and Sarel Fleishman, 2017)
* `score5_smooth`  smooth variant of score 5
* `score_membrane`  file: score_membrane.wts -  @desc: Scores low resolution pairwise, packing density, and environment energies -  in a 5-layer membrane environment. Also includes penalties for structures unfavorably -  embedded in the membrane. -  @note: Deprecated: Consider using mpframework_cen_2014 (2014) -  @cite: [Yarov-Yaravoy, 2006, Proteins]
* `score_membrane_cen` 
* `score_membrane_center_normal` 
* `score_membrane_env` 
* `scorefacts` 
* `scorefacts_cart` 
* `small_radii` 
* `smooth_etable` 
* `soft_rep` 
* `soft_rep_design` 
* `soft_rep_gen_born` 
* `sp2_correction` Weights file to be used with the hbond sp2 correction version 11c. - This is score12prime with Sp2 H-Bond correction version 11c -  * short range bb hbonds are not down weighted by a half, as they are in score12, but rather have the same strength as the long-range bb/bb hbonds. -  * reference energies refit using OptE (Andrew's OptE run 209) -  * the solvation term, increased by 15% -  * should be used in tandem with these flag -    -corrections:hbond_sp2_correction -  optE209 jsr55_pnataa_sc12dun10_sp2v11_olf11c_hesm16_fmx6_nofapair_hbenv_refonly_final
* `sp2_correction_dun02` Weights file to be used with the hbond sp2 correction version 11c. - This is score12prime with Sp2 H-Bond correction version 11c -  * short range bb hbonds are not down weighted by a half, as they are in score12, but rather have the same strength as the long-range bb/bb hbonds. -  * reference energies refit using OptE (Andrew's OptE run 208) -  * the solvation term, increased by 15% -  * should be used in tandem with these flags -    -corrections:hbond_sp2_correction -    -dun10 false -  optE208 jsr55_pnataa_sc12bicubic_sp2v11_olf11c_hesm16_fmx6_nofapair_hbenv_refonly_final
* `sp2_correction_dun02_no_bicubic_spline` Weights file to be used with the hbond sp2 correction version 11c. - This is score12prime with Sp2 H-Bond correction version 11c -  * short range bb hbonds are not down weighted by a half, as they are in score12, but rather have the same strength as the long-range bb/bb hbonds. -  * reference energies refit using OptE (Andrew's OptE run 213) -  * the solvation term, increased by 15% -  * should be used in tandem with these flags -    -corrections:hbond_sp2_correction -    -use_bicubic_interpolation false -    -dun10 false - optE213 jsr55_pnataa_sc12_sp2v11_olf11c_hesm16_fmx6_nofapair_hbenv_refonly_final
* `talaris2013`  The Talaris2013 score function combines several improvements to the previous -  default score function, score12: the 2010 Dunbrack Rotamer Library, -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  This set of defaults was elected at the energy function meeting at the -  Talaris conference center in Seattle in May 2013. They were demoted -  from default at the March 2016 WinterRosettaCon, in favor of talaris2014. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.4%) was -  chosen. - 
* `talaris2013_calibrated`  The Talaris2013 score function combines several improvements to the previous -  default score function, score12: the 2010 Dunbrack Rotamer Library, -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.4%) was -  chosen. -  - The following weights were modified bya scaling factor of 1/0.83 based on the testing Matt omera did - Talaris weights - fa_atr 0.8 - fa_rep 0.44 - fa_sol 0.75 - fa_intra_rep 0.004 - fa_elec 0.7 - pro_close 1 - hbond_sr_bb 1.17 - hbond_lr_bb 1.17 - hbond_bb_sc 1.17 - hbond_sc 1.1 - dslf_fa13 1.0 - rama 0.2 - omega 0.5 - fa_dun 0.56 - p_aa_pp 0.32 - ref 1
* `talaris2013_cart`  The Talaris2013 score function combines several improvements to the previous -  default score function, score12: the 2010 Dunbrack Rotamer Library, -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.4%) was -  chosen. -  -  -- talaris2013_cart -- -  cart_bonded 0.5 and pro_close 0 -  Reference energies adjusted on a subset of Ding & Dokholyan 38 set using -  bonded params checked in with this file (HIS_D added, PRO backbone same as others)
* `talaris2013_cst`  The Talaris2013 score function combines several improvements to the previous -  default score function, score12: the 2010 Dunbrack Rotamer Library, -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.4%) was -  chosen. -  -  This _cst version of Talaris assigns a weight of 1.0 to constraints, and is otherwise -  identical to the talaris set. To be used with protocols that use constraints (e.g. enzdes). - 
* `talaris2013_downAla`  The Talaris2013 score function combines several improvements to the previous -  default score function, score12: the 2010 Dunbrack Rotamer Library, -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.4%) was -  chosen. -  - DASM changed  weights: 0.592942 0.354993 -1.28682 -1.55374 0.43057 0.140526 0.357498 0.831803 -0.287374 0.602328 0.158677 -0.94198 -0.219285 -1.17797 -0.14916 0.176583 0.16454 0.744844 0.92933 0.131696 - DASM to downweight ala (they are in alphabetycal order based on one letter code):  
* `talaris2013_dun02`  This weight set was generated with all of the talaris2013 improvements -  except for the 2010 rotamer library.  The reference energies here -  were fit with the 2002 rotamer library.  That is, the score function -  improvements (over score12) that it contains are: -  -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  If you need to use this weight set, then use these flags on your command line: -  -dun10 false -  -score:weights talaris2013_dun02.wts -  -  Where the second flag (-score:weights) may not be neccessary if you're using -  RosettaScripts and have specified your score function from an XML file. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.6%) was -  chosen. - 
* `talaris2013_occ_sol`  The Talaris2013 score function combines several improvements to the previous -  default score function, score12: the 2010 Dunbrack Rotamer Library, -  the sp2 hydrogen bond potential, an explicit electrostatics term with a -  distance dependent dielectric (and a removal of the previous knowledge- -  based electrostatic potential, fa_pair), an adjustment to the LK_DGFREE -  parameters for four atom types, the 05.2009 ideal coordinates for -  the amino acids, an expansion of hydroxyl sampling for serine and -  threonine, the use of bicubic-spline interpolation in our knowledge- -  based potentials, an improved disulfide potential, and an analytic -  evaluation of our Lennard-Jones and EEF1 potentials. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.4%) was -  chosen. -  -  JAB - Talaris2013 with occ_Hbond_sol_fitted patch from Karanicolas Lab - Temporary file for antibody design analysis.
* `talaris2013_split_unfolded` Alternate version of the Talaris2013 weights to test the split unfolded energy term.
* `talaris2014`  The Talaris2014 score function represents a small modification to the -  Talaris2013 score function: the weights are all scaled upwards so that -  fa_atr has a weight of 1.0; then the hydrogen bond strengths are decreased -  by 20%.  This has the effect of keeping the hbond strengths fixed while -  everything else gets stronger. -  -  The benchmarking performed for the O'Meara et al. (2014) hbond paper -  showed that weakening the hbond weights by 20% improved sequence recovery, -  rotamer recovery, and decoy discrimination.  This weight set is -  currently the official gold standard weight set for Rosetta. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.8%) was -  chosen. - 
* `talaris2014_cart`  The Talaris2014 score function represents a small modification to the -  Talaris2013 score function: the weights are all scaled upwards so that -  fa_atr has a weight of 1.0; then the hydrogen bond strengths are decreased -  by 20%.  This has the effect of keeping the hbond strengths fixed while -  everything else gets stronger. The Talaris2014_cart weight set replaces -  the pro_close term with cart_bonded. -  -  The benchmarking performed for the O'Meara et al. (2014) hbond paper -  showed that weakening the hbond weights by 20% improved sequence recovery, -  rotamer recovery, and decoy discrimination. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.1%) was -  chosen. - 
* `talaris2014_cst`  The Talaris2014 score function represents a small modification to the -  Talaris2013 score function: the weights are all scaled upwards so that -  fa_atr has a weight of 1.0; then the hydrogen bond strengths are decreased -  by 20%.  This has the effect of keeping the hbond strengths fixed while -  everything else gets stronger. -  -  The benchmarking performed for the O'Meara et al. (2014) hbond paper -  showed that weakening the hbond weights by 20% improved sequence recovery, -  rotamer recovery, and decoy discrimination. -  -  Reference energies were fit with optE on Jane Richardson's HiQ54 benchmark -  set in triplicate, and tested on the Ding & Dokholyan 38 set.  The -  set of reference energies with the highest sequence recovery (39.8%) was -  chosen. -  -  This _cst version of talaris2014 assigns a weight of 1.0 to constraints, and is otherwise -  identical to the talaris set. To be used with protocols that use constraints (e.g. enzdes).

Specifically:

* `beta` &larr;  old beta
* `beta_genpot` &larr; now `genpot_full`
* `beta_july15` &larr; prequel to ref2015
* `beta_nov15` &larr; synonym of ref2015
* `beta_nov16` &larr; ???
* `corrections_conway2016` &larr; ???
* `franklin2019` &larr; membrane
* `ref2015`  &larr; main
* `score0` &larr; ancient (rosetta++)
* `score1` &larr; ancient (rosetta++)
* `score12` &larr; old
* `score13` &larr; old
* `score2` &larr; ancient (rosetta++)
* `score3` &larr; ancient (rosetta++)
* `score5` &larr; ancient (rosetta++)
* `talaris2013` &larr; old
* `talaris2014` &larr; old

Two scorefunctions that stand out are `corrections_conway2016` and `beta_nov16`.

Unfortunately, the suggested flag for `corrections_conway2016` does not work.

```python
pyrosetta.rosetta.basic.options.set_boolean_option('dun10_dir rotamer/corrections_conway2016', True)  #RuntimeError
pyrosetta.rosetta.basic.options.set_boolean_option('dun10_dir rotamer', True)  #RuntimeError
pyrosetta.rosetta.basic.options.set_boolean_option('corrections_conway2016', True)  #RuntimeError
```

However, looking at the paper it seems not relevant and an offshoot of Talaris2014 scores.

```python
print(ww.compare(['score12', 'talaris2013', 'talaris2014', 'ref2015', 'corrections_conway2016', 'beta_nov16']).to_markdown())
```

| name                   |   fa_atr |   fa_rep |   fa_sol |   fa_intra_rep |   pro_close |   fa_pair |   hbond_sr_bb |   hbond_lr_bb |   hbond_bb_sc |   hbond_sc |   dslf_ss_dst |   dslf_cs_ang |   dslf_ss_dih |   dslf_ca_dih |   rama |   omega |   fa_dun |   p_aa_pp |   ref |   ref_ALA |   ref_ARG |   ref_ASN |   ref_ASP |   ref_CYS |   ref_GLN |   ref_GLU |   ref_GLY |   ref_HIS |   ref_ILE |   ref_LEU |   ref_LYS |   ref_MET |   ref_PHE |   ref_PRO |   ref_SER |   ref_THR |   ref_TRP |   ref_TYR |   ref_VAL |   fa_elec |   dslf_fa13 |   yhh_planarity |   fa_intra_sol_xover4 |   lk_ball_wtd |   rama_prepro |   fa_intra_atr_xover4 |   fa_intra_rep_xover4 |   fa_dun_dev |   fa_dun_rot |   fa_dun_semi |   lk_ball |   lk_ball_iso |   lk_ball_bridge |   lk_ball_bridge_uncpl |   fa_intra_elec |   hxl_tors |
|:-----------------------|---------:|---------:|---------:|---------------:|------------:|----------:|--------------:|--------------:|--------------:|-----------:|--------------:|--------------:|--------------:|--------------:|-------:|--------:|---------:|----------:|------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|------------:|----------------:|----------------------:|--------------:|--------------:|----------------------:|----------------------:|-------------:|-------------:|--------------:|----------:|--------------:|-----------------:|-----------------------:|----------------:|-----------:|
| score12                |      0.8 |     0.44 |   0.65   |          0.004 |        1    |      0.49 |         0.585 |          1.17 |          1.17 |        1.1 |             1 |             1 |             1 |             1 |   0.2  |   0.5   |     0.56 |      0.32 |     1 |  0.16     |  -0.98    | -0.89     |  -0.67    |  1.7      |  -0.97    |  -0.81    | -0.17     |  0.56     |  0.24     | -0.1      | -0.65     | -0.34     |   0.63    |  0.02     | -0.37     |  -0.27    |   0.91    |  0.51     |  0.29     |     0     |        0    |           0     |                0      |             0 |          0    |                     0 |                  0    |        0     |         0    |          0    |      0    |          0    |             0    |                   0    |               0 |          0 |
| talaris2013            |      0.8 |     0.44 |   0.75   |          0.004 |        1    |      0    |         1.17  |          1.17 |          1.17 |        1.1 |             0 |             0 |             0 |             0 |   0.2  |   0.5   |     0.56 |      0.32 |     1 |  0.592942 |  -0.14916 | -0.94198  |  -1.28682 |  0.354993 |  -1.17797 |  -1.55374 |  0.140526 |  0.357498 |  0.831803 |  0.602328 | -0.287374 |  0.158677 |   0.43057 | -0.219285 |  0.176583 |   0.16454 |   0.92933 |  0.131696 |  0.744844 |     0.7   |        1    |           0     |                0      |             0 |          0    |                     0 |                  0    |        0     |         0    |          0    |      0    |          0    |             0    |                   0    |               0 |          0 |
| talaris2014            |      1   |     0.55 |   0.9375 |          0.005 |        1.25 |      0    |         1.17  |          1.17 |          1.17 |        1.1 |             0 |             0 |             0 |             0 |   0.25 |   0.625 |     0.7  |      0.4  |     1 |  0.773742 |  -0.32436 | -1.19118  |  -1.63002 |  0.443793 |  -1.51717 |  -1.96094 |  0.173326 |  0.388298 |  1.0806   |  0.761128 | -0.358574 |  0.249477 |   0.61937 | -0.250485 |  0.165383 |   0.20134 |   1.23413 |  0.162496 |  0.979644 |     0.875 |        1.25 |           0.625 |                0      |             0 |          0    |                     0 |                  0    |        0     |         0    |          0    |      0    |          0    |             0    |                   0    |               0 |          0 |
| ref2015                |      1   |     0.55 |   1      |          0.005 |        1.25 |      0    |         1     |          1    |          1    |        1   |             0 |             0 |             0 |             0 |   0    |   0.4   |     0.7  |      0.6  |     1 |  1.32468  |  -0.09474 | -1.34026  |  -2.14574 |  3.25479  |  -1.45095 |  -2.72453 |  0.79816  | -0.30065  |  2.30374  |  1.66147  | -0.71458  |  1.65735  |   1.21829 | -1.64321  | -0.28969  |   1.15175 |   2.26099 |  0.58223  |  2.64269  |     1     |        1.25 |           0.625 |                1      |             1 |          0.45 |                     0 |                  0    |        0     |         0    |          0    |      0    |          0    |             0    |                   0    |               0 |          0 |
| corrections_conway2016 |      1   |     0.55 |   0.9375 |          0.005 |        1.25 |      0    |         1.17  |          1.17 |          1.17 |        1.1 |             0 |             0 |             0 |             0 |   0.25 |   0.625 |     0    |      0.4  |     1 |  0.773742 |  -0.32436 | -1.19118  |  -1.63002 |  0.443793 |  -1.51717 |  -1.96094 |  0.173326 |  0.388298 |  1.0806   |  0.761128 | -0.358574 |  0.249477 |   0.61937 | -0.250485 |  0.165383 |   0.20134 |   1.23413 |  0.162496 |  0.979644 |     0.875 |        1.25 |           0.625 |                0.9375 |             0 |          0    |                     1 |                  0.55 |        0.375 |         0.7  |          0.7  |      0    |          0    |             0    |                   0    |               0 |          0 |
| beta_nov16             |      1   |     0.55 |   1      |          0     |        1.25 |      0    |         1     |          1    |          1    |        1   |             0 |             0 |             0 |             0 |   0    |   0.48  |     0    |      0.61 |     1 |  2.3386   |  -1.281   | -0.873554 |  -2.2837  |  3.2718   |  -1.0644  |  -2.5358  |  1.2108   |  0.134426 |  1.0317   |  0.729516 | -1.6738   |  1.2334   |   1.4028  | -5.1227   | -1.1772   |  -1.425   |   3.035   |  0.964136 |  2.085    |     1     |        1.25 |           0     |                1      |             0 |          0.5  |                     1 |                  0.55 |        0.69  |         0.76 |          0.78 |      0.92 |         -0.38 |            -0.33 |                  -0.33 |               1 |          1 |

Regarding `beta_nov16`, it seems to be a progression of the `beta` series. So is intriguing.

For the full list of weights in tabular format see [weights.csv](data/weights.csv), which was generated via:

```python
import re
from pyrosetta_help import WeightWatcher
ww = WeightWatcher()

data = []
residue_order = ['ref_ALA',  # single letter alphabetical
                 'ref_CYS',
                 'ref_ASP',
                 'ref_GLU',
                 'ref_PHE',
                 'ref_GLY',
                 'ref_HIS',
                 'ref_ILE',
                 'ref_LYS',
                 'ref_LEU',
                 'ref_MET',
                 'ref_ASN',
                 'ref_PRO',
                 'ref_GLN',
                 'ref_ARG',
                 'ref_SER',
                 'ref_THR',
                 'ref_VAL',
                 'ref_TRP',
                 'ref_TYR']
for name in ww.possible_scorefxn_names:
    block = ww.get_scorefxn_block(name)
    years_mentioned = sorted(set(re.findall(r'(?<=\D)20\d\d(?=\D)', block)))
    datum = dict(name=name,
                 comments = '',
                 max_year_mentioned = max(years_mentioned) if len(years_mentioned) else float('nan')
                )
    for line in block.split('\n'):
        line = line.strip()
        rex1 = re.match(r'(\w+)\s+([\.\d]+)', line)
        rex2 = re.match(r'(\w+)\s+(.*)', line)
        if not line:
            continue
        elif line[0] == '#':
            datum['comments'] += line.replace('#', '').strip()+'\n'
        elif 'METHOD_WEIGHTS ref ' in line:
            datum = {**datum, **dict( zip(residue_order, line.split()[2:]) )}
        elif rex1:
            datum[rex1.group(1)] = float(rex1.group(2))
        elif rex2:
            datum[rex2.group(1)] = rex2.group(2)
        elif ' ' in line:
            print(line)
            raise Exception
        else:
            datum[line] = True
    # end loop per line
    data.append(datum)
# end loop per name

import pandas as pd
weights = pd.DataFrame(data)

def get_derivation(row):
    if str(row['fa_sol']) == 'nan':
        return 'rosetta++'
    elif row['fa_sol'] == 0.6500:
        return 'score12'
    elif row['fa_sol'] == 0.7500:
        return 'talaris2013'
    elif row['fa_sol'] == 0.9375:
        return 'talaris2014'
    elif row['fa_sol'] != 1.0000:
        print(row['name'], row['fa_sol'], row['fa_atr'], row['fa_rep'])
        return 'unknown'
    elif row['fa_sol'] == 1.16497:
        return 'mm_std'
    elif row['fa_sol'] == 0.583958:
        return 'orbitals'
    else:
        return 'ref2015'
    
weights['suspected_derivation'] = weights.apply(get_derivation, 1)

def w2(name):
    if 'fa_' in name:
        return 0
    elif 'ref_' in name:
        return 1
    elif '_constraint' in name:
        return 0.5
    elif name.upper() == name:
        return 20
    else:
        return 10

nice = sorted(weights.columns.values, key=w2)
    
weights[nice].to_csv('/Users/matteo/Coding/weights.csv')
```

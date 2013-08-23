MUSIC Observing Run: August 2013
================================

Observing Stuff
---------------

Just files needed for observing:

 * clean_params_science_music.txt
 * mapping_params_science_music.txt
 * login_info.txt.gpg
 * targets.txt [formatted for UIP]
 * targets.reg

Data Examination / Analysis
---------------------------

There are scripts in the scripts and individual observing nights' directories.  

seds.py:
Create SEDs for all files following the YYMMDD_ob# naming scheme in subdirectories, assuming
they're MUSIC auto-reduced map files in IDLsave format.

convolve_match_makefits.py:
Create beam-matched maps of each file sampled onto the highest-resolution grid
with the lowest-resolution beam size.  Optionally, save each of these and their
unsharp-masked companions as .fits files.

viewer.py:
Some trivial but handy quick-look data viewing scripts.

#!/usr/bin/env python3

import os
import subprocess
import glob
import sys
import shutil # For checking executable existence and file operations
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# Removed 'import pandas as pd' - no longer needed

# --- Configuration ---

# 1. Define the sequences to extract per group
#    Key: Name for the analysis group (used for output filenames)
#    Value: List of exact sequence IDs to extract.
#           IMPORTANT: Assume the FIRST ID in each list is the 'query' or reference
#                      for that gene group, used if the FASTA header lacks its symbol.
sequences_to_extract = {
    'Cavin3': [
    # Query
    'ENSMUSP00000044979.3', # Mouse Cavin3

    # Mammalian Cavin3 Orthologs (Best hits)
    'ENSFCAP00000035976.3', # Cat
    'ENSP00000307292.3',   # Human
    'ENSLAFP00000009191.3', # Elephant
    'ENSEEUP00000011563.1', # Hedgehog
    'ENSBTAP00000026323.5', # Cattle
    'ENSMODP00000023188.3', # Opossum

    # Mammalian Cavin2 Paralogs (Best hits per species)
    'ENSFCAP00000040014.1', # Cat
    'ENSP00000305675.4',   # Human
    'ENSBTAP00000024617.3', # Cattle
    'ENSMODP00000014872.1', # Opossum
    'ENSLAFP00000025451.1', # Elephant
    'ENSEEUP00000008600.1', # Hedgehog

    # Non-mammalian Best Hits (All are Cavin2)
    'ENSPMRP00000005185.1', # Lizard (Cavin2)
    'ENSXETP00000107494.1', # Frog (Cavin2)
    'ENSGALP00010014540.1', # Chicken (Cavin2)
    'ENSDARP00000121255.1', # Zebrafish (cavin2b - Best Zebrafish Hit)
    'ENSDARP00000153266.1', # Zebrafish (cavin2a - Representative)


    # Representative Cavin1 Paralogs (Mammals & Non-mammals)
    'ENSMODP00000018587.2', # Opossum
    'ENSGALP00010030535.1', # Chicken (PTRF/Cavin1)
    'ENSLAFP00000004735.3', # Elephant
    'ENSFCAP00000041960.1', # Cat
    'ENSBTAP00000102795.1', # Cattle
    'ENSPMRP00000025014.1', # Lizard
    'ENSP00000349541.4',   # Human
    'ENSXETP00000101996.1', # Frog
    'ENSDARP00000151250.1', # Zebrafish (cavin1b - Representative)
    'ENSDARP00000099909.3', # Zebrafish (cavin1a - Representative)

    # Representative Cavin4 Paralogs (Mammals & Non-mammals)
    'ENSFCAP00000041458.1', # Cat
    'ENSEEUP00000014249.1', # Hedgehog
    'ENSMODP00000059422.1', # Opossum
    'ENSBTAP00000029327.4', # Cattle
    'ENSP00000418668.1',   # Human
    'ENSGALP00010011040.1', # Chicken
    'ENSPMRP00000016921.1', # Lizard
    'ENSDARP00000121336.1', # Zebrafish (cavin4a - Representative)
    'ENSDARP00000073748.4', # Zebrafish (cavin4b - Representative)
],
     'Dsg3_Analysis': [    # Query
         # Query
         'ENSMUSP00000064718.7',  # Mouse DSG3

         # Mammalian DSG3 Orthologs
         'ENSFCAP00000022219.4',  # Cat DSG3
         'ENSP00000257189.4',  # Human DSG3
         'ENSBTAP00000002538.4',  # Cattle DSG3
         'ENSMODP00000013762.4',  # Opossum DSG3
         'ENSEEUP00000003624.1',  # Hedgehog DSG3
         # Elephant DSG3 not found

         # Mammalian DSG4 Paralogs
         'ENSMODP00000046333.1',  # Opossum DSG4
         'ENSP00000311859.4',  # Human DSG4
         'ENSFCAP00000017616.2',  # Cat DSG4
         'ENSLAFP00000014865.3',  # Elephant DSG4
         'ENSBTAP00000045458.2',  # Cattle DSG4
         'ENSEEUP00000011646.1',  # Hedgehog DSG4

         # Mammalian DSG1 Paralogs
         'ENSMODP00000013830.4',  # Opossum DSG1
         'ENSP00000257192.4',  # Human DSG1
         'ENSLAFP00000014857.3',  # Elephant DSG1
         'ENSFCAP00000049125.2',  # Cat DSG1
         'ENSBTAP00000018382.5',  # Cattle DSG1
         'ENSEEUP00000007320.1',  # Hedgehog DSG1

         # Mammalian DSG2 Paralogs
         'ENSP00000261590.8',  # Human DSG2
         'ENSBTAP00000069779.2',  # Cattle DSG2
         'ENSLAFP00000024941.1',  # Elephant DSG2
         'ENSFCAP00000024430.3',  # Cat DSG2
         'ENSMODP00000013697.4',  # Opossum DSG2
         # Hedgehog DSG2 not found

         # Non-Mammalian Best Hits / Key Paralogs
         'ENSGALP00010005873.1',  # Chicken DSG2 (Best Chicken Hit)
         'ENSGALP00010005862.1',  # Chicken DSG1/related
         'ENSPMRP00000007925.1',  # Lizard DSG2 (Best Lizard Hit)
         'ENSPMRP00000007949.1',  # Lizard DSG1/related
         'ENSXETP00000074721.2',  # Frog DSG1 (Best Frog Hit)
         'ENSXETP00000104811.1',  # Frog DSG2
         # Zebrafish hits missing
         ],

    'Lalba_Analysis': [    # Query
    'ENSMUSP00000023726.4', # Mouse LALBA

    # Mammalian LALBA Orthologs
    'ENSP00000301046.2',   # Human LALBA
    'ENSFCAP00000056431.1', # Cat LALBA
    'ENSBTAP00000007701.2', # Cattle LALBA
    'ENSLAFP00000011887.2', # Elephant LALBA
    'ENSMODP00000026437.3', # Opossum LALBA

    # Key Lysozyme Paralogs (Mammalian)
    'ENSBTAP00000031167.5', # Cattle LYZ/related
    'ENSBTAP00000101980.1', # Cattle LYZ/related (another locus?)
    'ENSBTAP00000027401.3', # Cattle LYSB
    'ENSBTAP00000015846.4', # Cattle LYZF1
    'ENSBTAP00000013640.5', # Cattle LYZL1
    'ENSBTAP00000007924.4', # Cattle LYZ2
    'ENSBTAP00000027399.5', # Cattle LYZ3
    'ENSBTAP00000038081.2', # Cattle LYZ
    'ENSFCAP00000027271.1', # Cat LYZ/related
    'ENSFCAP00000005029.2', # Cat LYZ
    'ENSP00000261267.2',   # Human LYZ
    'ENSP00000364650.3',   # Human LYZL1
    'ENSP00000364467.2',   # Human LYZL2
    'ENSLAFP00000005203.2', # Elephant LYZ/related
    'ENSMODP00000037594.2', # Opossum LYZ/related

    # Non-Mammalian Outgroups (Best hits - Lysozymes)
    'ENSPMRP00000019718.1', # Lizard LYZ/related (Best hit)
    'ENSGALP00010016495.1', # Chicken LYZ (Best hit)
    'ENSXETP00000109344.1'], # Frog LYZ/related (Best hit)],
    'Pip_Analysis': ['ENSMUSP00000074818.8', # Mouse Pip (assuming this is query)
                   'ENSBTAP00000008743.4', 'ENSFCAP00000044419.2', 'ENSLAFP00000005939.2',
                   'ENSLAFP00000023013.1', 'ENSP00000291009.3'],
    'Plet1_Analysis': ['ENSMUSP00000110118.2', # Mouse Plet1 (assuming this is query)
                     'ENSBTAP00000020999.6', 'ENSFCAP00000060825.1', 'ENSMODP00000038360.2',
                     'ENSP00000341412.2'],
    'TOP2A_Analysis': [
    # Query
    'ENSMUSP00000068896.8', # Mouse TOP2A

    # Mammalian TOP2A Orthologs
    'ENSLAFP00000027558.1', # Elephant TOP2A
    'ENSP00000411532.1',   # Human TOP2A
    'ENSFCAP00000008944.4', # Cat TOP2A
    'ENSBTAP00000053858.3', # Cattle TOP2A
    'ENSMODP00000016461.3', # Opossum TOP2A
    # Hedgehog TOP2A not found

    # Non-Mammalian TOP2A Orthologs (!!! Found matching symbols !!!)
    'ENSGALP00010033066.1', # Chicken TOP2A
    'ENSPMRP00000025907.1', # Lizard TOP2A
    'ENSXETP00000086934.2', # Frog top2a
    'ENSDARP00000139353.2', # Zebrafish top2a

    # Mammalian TOP2B Paralogs
    'ENSFCAP00000020385.3', # Cat TOP2B
    'ENSBTAP00000089599.1', # Cattle TOP2B
    'ENSP00000396704.2',   # Human TOP2B
    'ENSLAFP00000013292.3', # Elephant TOP2B
    'ENSMODP00000050300.1', # Opossum TOP2B
    'ENSEEUP00000005775.1', # Hedgehog TOP2B

    # Non-Mammalian TOP2B Paralogs
    'ENSGALP00010003508.1', # Chicken TOP2B
    'ENSXETP00000064526.2', # Frog top2b
    'ENSPMRP00000032571.1', # Lizard TOP2B
    'ENSDARP00000108088.2', # Zebrafish top2b
],
'Art3':[    # ── Mouse ──  (query 基准)
    "ENSMUSP00000113493.2",    # Mouse ART3

    # ── 哺乳动物：ART3 正交 ──
    "ENSP00000514990.1",       # Human ART3
    "ENSFCAP00000061968.1",    # Cat ART3
    "ENSBTAP00000087902.1",    # Cattle ART3
    "ENSLAFP00000013864.3",    # Elephant ART3
    "ENSEEUP00000012185.1",    # Hedgehog ART3
    "ENSMODP00000024079.3",    # Opossum ART3

    # ── 家族内 paralog（ART5，取 Human 为代表）──
    "ENSP00000352992.4",       # Human ART5   ← 检测 ART3/5 分化节点

    # ── 外群：爬-鸟-两栖 ──
    "ENSGALP00010004346.1",    # Chicken ART3  ← 直接高分命中
    "ENSXETP00000080186.2",    # Frog art5     ← 最近家族 paralog
    "ENSGALP00010003273.1",    # Chicken ART1A ← 深外群 (ART1/PTRF 类)
    "ENSXETP00000000014.5"],    # Frog art1     ← 深外群 (ART1)]
    'Afm':[
    # Query
    'ENSMUSP00000108804.3', # Mouse AFM

    # Mammalian AFM Orthologs
    'ENSP00000226355.3',   # Human AFM
    'ENSFCAP00000026478.2', # Cat AFM
    'ENSBTAP00000094759.1', # Cattle AFM
    'ENSLAFP00000009646.2', # Elephant AFM (provisional)
    'ENSMODP00000023933.3', # Opossum AFM
    # Hedgehog AFM not found

    # Mammalian AFP Paralogs
    'ENSBTAP00000022772.3', # Cattle AFP
    'ENSP00000226359.2',   # Human AFP
    'ENSFCAP00000021758.1', # Cat AFP
    'ENSMODP00000038307.2', # Opossum AFP
    'ENSLAFP00000027657.1', # Elephant AFP
    'ENSEEUP00000002410.1', # Hedgehog AFP

    # Mammalian ALB Paralogs
    'ENSP00000295897.4',   # Human ALB
    'ENSFCAP00000028536.3', # Cat ALB
    'ENSBTAP00000022763.5', # Cattle ALB
    'ENSLAFP00000000908.4', # Elephant ALB
    'ENSMODP00000023917.2', # Opossum ALB
    'ENSEEUP00000006959.1', # Hedgehog ALB

    # Mammalian GC Paralogs
    'ENSFCAP00000004482.6', # Cat GC
    'ENSBTAP00000033386.4', # Cattle GC
    'ENSMODP00000055418.1', # Opossum GC
    'ENSP00000426683.1',   # Human GC
    'ENSLAFP00000013695.2', # Elephant GC
    'ENSEEUP00000004340.1', # Hedgehog GC

    # Non-Mammalian Outgroups (Albuminoid family members)
    'ENSGALP00010006705.1', # Chicken AFP (Best Chicken Hit)
    'ENSGALP00010006681.1', # Chicken ALB
    'ENSGALP00010005944.1', # Chicken GC
    'ENSXETP00000050547.4', # Frog ALB (Best Frog Hit)
    'ENSXETP00000045132.4', # Frog GC
    'ENSDARP00000106674.3', # Zebrafish GC (Best Zebrafish Hit)
    # Lizard representatives missing from data
],
    'Icam2':[
    # Query
    'ENSMUSP00000001055.9', # Mouse ICAM2

    # Mammalian ICAM2 Orthologs
    'ENSP00000462579.1',   # Human ICAM2
    'ENSFCAP00000047272.2', # Cat ICAM2
    'ENSLAFP00000004183.3', # Elephant ICAM2
    'ENSMODP00000011683.3', # Opossum ICAM2
    'ENSBTAP00000073014.2', # Cattle ICAM2 (candidate)
    'ENSEEUP00000004748.1', # Hedgehog ICAM2 (candidate)

    # Mammalian Paralogs (ICAM1, ICAM3, ICAM4, ICAM5, VCAM1, TCAM1)
    'ENSBTAP00000025884.7', # Cattle ICAM related
    'ENSBTAP00000034654.3', # Cattle TCAM1
    'ENSBTAP00000076658.1', # Cattle ICAM1
    'ENSBTAP00000081446.1', # Cattle ICAM5
    'ENSBTAP00000020903.2', # Cattle ICAM3
    'ENSFCAP00000025193.2', # Cat ICAM related
    'ENSFCAP00000019585.3', # Cat ICAM3
    'ENSFCAP00000008749.5', # Cat ICAM5
    'ENSFCAP00000008747.5', # Cat ICAM1
    'ENSP00000160262.3',   # Human ICAM3
    'ENSP00000221980.3',   # Human ICAM5
    'ENSP00000264832.2',   # Human ICAM1
    'ENSP00000370147.2',   # Human ICAM4
    'ENSP00000359133.1',   # Human VCAM1
    'ENSMODP00000047483.1', # Opossum ICAM5
    'ENSMODP00000040597.1', # Opossum ICAM1-like
    'ENSMODP00000053020.1', # Opossum ICAM3
    'ENSLAFP00000019338.1', # Elephant ICAM related
    'ENSLAFP00000018446.1', # Elephant ICAM5
    'ENSEEUP00000012508.1', # Hedgehog ICAM5
    'ENSEEUP00000010953.1', # Hedgehog VCAM1

    # Non-Mammalian Outgroups
    'ENSPMRP00000036814.1', # Lizard ICAM related (Best Lizard hit)
    'ENSXETP00000108710.1', # Frog ICAM related (Best Frog hit)
],
    'S100a3':[
    # Query
    'ENSMUSP00000142747.2', # Mouse S100A3

    # Mammalian S100A3 Orthologs
    'ENSP00000357701.1',   # Human S100A3
    'ENSLAFP00000025329.1', # Elephant S100A3
    'ENSMODP00000052348.1', # Opossum S100A3
    # Cattle, Cat, Hedgehog S100A3 not found

    # Critical Non-Mammalian S100A3 Hit
    'ENSPMRP00000022660.1', # Lizard S100A3 !!! Found matching symbol !!!

    # Mammalian Paralogs (Representative selection)
    'ENSLAFP00000012560.2', # Elephant S100A4
    'ENSLAFP00000018112.1', # Elephant S100A2
    'ENSLAFP00000025087.1', # Elephant S100A6
    'ENSLAFP00000002220.2', # Elephant S100A1
    'ENSLAFP00000012558.3', # Elephant S100A5
    'ENSLAFP00000012505.3', # Elephant S100B
    'ENSLAFP00000010903.1', # Elephant S100Z
    'ENSBTAP00000000589.3', # Cattle S100A2
    'ENSBTAP00000097885.1', # Cattle S100A4
    'ENSBTAP00000006806.4', # Cattle S100A1
    'ENSBTAP00000081557.1', # Cattle S100A5
    'ENSBTAP00000063447.1', # Cattle S100B
    'ENSBTAP00000026904.5', # Cattle S100Z
    'ENSFCAP00000018025.3', # Cat S100A2
    'ENSFCAP00000053589.1', # Cat S100A4
    'ENSFCAP00000023487.1', # Cat S100A1
    'ENSFCAP00000013040.4', # Cat S100A5
    'ENSFCAP00000041539.1', # Cat S100A6 (No Symbol, but likely)
    'ENSFCAP00000036720.1', # Cat S100B
    'ENSFCAP00000048112.1', # Cat S100P
    'ENSFCAP00000041362.2', # Cat S100Z
    'ENSP00000357703.1',   # Human S100A4
    'ENSP00000357698.2',   # Human S100A2
    'ENSP00000357709.1',   # Human S100A6
    'ENSP00000292169.2',   # Human S100A1
    'ENSP00000357706.2',   # Human S100A5
    'ENSP00000380769.1',   # Human S100B
    'ENSP00000296370.3',   # Human S100P
    'ENSP00000483535.1',   # Human S100Z
    'ENSP00000271638.2',   # Human S100A11
    'ENSMODP00000050064.1', # Opossum S100A4
    'ENSMODP00000021678.3', # Opossum S100A1
    'ENSMODP00000051818.1', # Opossum S100A5/related
    'ENSMODP00000024649.2', # Opossum S100Z
    'ENSEEUP00000013800.1', # Hedgehog S100A1

    # Non-Mammalian Outgroups (Representatives of S100 family)
    'ENSGALP00010042260.1', # Chicken S100A4 (Best Chicken Hit)
    'ENSGALP00010015659.1', # Chicken S100Z
    'ENSGALP00010042258.1', # Chicken S100A6
    'ENSGALP00010042314.1', # Chicken S100A1
    'ENSGALP00010019174.1', # Chicken S100B
    'ENSGALP00010041048.1', # Chicken S100A11
    'ENSPMRP00000022675.1', # Lizard S100A5-like
    'ENSPMRP00000028423.1', # Lizard S100A1
    'ENSPMRP00000022682.1', # Lizard S100A6-like
    'ENSPMRP00000026190.1', # Lizard S100Z
    'ENSXETP00000037225.2', # Frog S100Z (Best Frog Hit)
    'ENSXETP00000086171.1', # Frog S100A1
    'ENSXETP00000063874.2', # Frog S100P
    'ENSDARP00000056560.5', # Zebrafish S100Z (Best Zebrafish Hit)
    'ENSDARP00000123790.1', # Zebrafish S100A1
    'ENSDARP00000072444.4', # Zebrafish S100T
    'ENSDARP00000074756.4', # Zebrafish S100B
    'ENSDARP00000033294.6', # Zebrafish S100A10b
],
    'Oas2':[
    # Query
    'ENSMUSP00000080209.7', # Mouse OAS2

    # Mammalian OAS2 Orthologs
    'ENSP00000376362.3',   # Human OAS2
    'ENSFCAP00000000294.4', # Cat OAS2
    'ENSBTAP00000068569.2', # Cattle OAS2
    'ENSLAFP00000011781.3', # Elephant OAS2
    'ENSMODP00000042792.1', # Opossum OAS2 (candidate)
    # Hedgehog OAS2 not clearly found

    # Mammalian OAS3 Paralogs
    'ENSLAFP00000011417.4', # Elephant OAS3
    'ENSMODP00000004072.3', # Opossum OAS3
    'ENSP00000228928.7',   # Human OAS3
    'ENSFCAP00000000293.4', # Cat OAS3
    'ENSEEUP00000004289.1', # Hedgehog OAS3

    # Mammalian OAS1 Paralogs (incl. variants Y/Z)
    'ENSMODP00000003962.3', # Opossum OAS1
    'ENSP00000449053.2',   # Human OAS1
    'ENSLAFP00000004658.4', # Elephant OAS1
    'ENSBTAP00000042418.3', # Cattle OAS1 (candidate)
    'ENSFCAP00000000292.5', # Cat OAS1
    'ENSBTAP00000016856.3', # Cattle OAS1Y
    'ENSBTAP00000082376.1', # Cattle OAS1Z
    'ENSEEUP00000010728.1', # Hedgehog OAS1/OASL (candidate)

    # Mammalian OASL Paralogs
    'ENSP00000257570.4',   # Human OASL
    'ENSMODP00000031339.3', # Opossum OASL
    'ENSFCAP00000008943.4', # Cat OASL
    'ENSLAFP00000015000.2', # Elephant OASL
    'ENSBTAP00000004270.5', # Cattle OASL

    # Non-Mammalian Outgroups
    'ENSPMRP00000024205.1', # Lizard OAS1A-like (Best Lizard Hit)
    'ENSGALP00010040957.1', # Chicken OASL (Best Chicken Hit)
    # Frog, Zebrafish hits not significant
],
    'Csn2':[    # Query
    'ENSMUSP00000143341.2', # Mouse CSN2

    # Mammalian CSN2 Orthologs
    'ENSLAFP00000012348.3', # Elephant CSN2
    'ENSFCAP00000000593.5', # Cat CSN2
    'ENSBTAP00000094835.1', # Cattle CSN2
    'ENSP00000341030.3',], # Human CSN2
'Ltf':[    "ENSFCAP00000057170.1",  # Cat LTF (-) best e-value
    "ENSFCAP00000039270.3",  # Cat LTF (2nd-best)
    "ENSFCAP00000020707.4",  # Cat TF   (non-LTF paralogue)

    "ENSBTAP00000078872.1",  # Cattle LTF (-) best
    "ENSBTAP00000001704.4",  # Cattle LTF (2nd)
    "ENSBTAP00000009564.6",  # Cattle TF

    "ENSGALP00010037107.1",  # Chicken TF

    "ENSLAFP00000006597.3",  # Elephant LTF (-)
    "ENSLAFP00000010419.3",  # Elephant TF-like (annotation lacks gene_symbol)

    "ENSXETP00000082148.1",  # Frog tf

    "ENSEEUP00000010966.1",  # Hedgehog TF-like

    "ENSP00000405719.2",     # Human LTF (-) best
    "ENSP00000405546.1",     # Human LTF (2nd)
    "ENSP00000385834.3",     # Human TF

    "ENSPMRP00000022817.1",  # Lizard LTF (-) best
    "ENSPMRP00000022826.1",  # Lizard LTF (2nd)
    "ENSPMRP00000010623.1",  # Lizard MELTF

    "ENSMODP00000044603.1",  # Opossum LTF (-) best
    "ENSMODP00000051639.1",  # Opossum LTF (2nd)
    "ENSMODP00000027349.3",  # Opossum MELTF

    "ENSDARP00000135247.2",  # Zebrafish tfa
]
}

# 2. List of Paths to the source peptide FASTA files
source_fasta_files = [
    # --- MODIFY THIS PATH TO YOUR FASTA FILE ---
    "/data02/sunxuebo/project/scrnaseq/mammaltree/data/blacklist/all_species.pep.fa"
]

# 3. Output Directories (Renumbered)
base_output_dir = "./phylo_pipeline_results_multi_source"
mafft_output_subdir = "mafft_alignments"
raxml_output_subdir = "raxml_trees"

# 4. Tool Executables and Parameters (Renumbered)
mafft_executable = "mafft"
raxml_executable = "raxmlHPC-PTHREADS"
mafft_threads = 20
mafft_options = ["--auto"]
raxml_threads = 20
raxml_bootstrap_reps = 500
raxml_model = "PROTGAMMALGF"
raxml_seed_bootstrap = 12345
raxml_seed_parsimony = 54321

# 5. Species Mapping (Prefix -> Species Name) (Renumbered)
prefix_to_species_map = {
    'ENSGALP': 'Chicken', 'ENSMODP': 'Opossum', 'ENSPMRP': 'Lizard',
    'ENSP': 'Human', 'ENSMUSP': 'Mouse', 'ENSXETP': 'Frog',
    'ENSBTAP': 'Cattle', 'ENSEEUP': 'Hedgehog', 'ENSFCAP': 'Cat',
    'ENSLAFP': 'Elephant',
    'ENSDARP':'Zebrafish'
}

# 6. Minimum sequences for RAxML (Renumbered)
MIN_SEQS_FOR_RAXML = 4

# --------------------- End of Configuration ---------------------

# --- Helper Functions ---
def run_command(command_list, cwd=None, check=True, stdout_file=None, stderr_capture=True):
    """Runs a shell command using subprocess, handles errors, captures stderr."""
    print(f"Running command: {' '.join(command_list)}")
    stderr_pipe = subprocess.PIPE if stderr_capture else None
    stdout_dest = stdout_file if stdout_file else None
    try:
        process = subprocess.run(command_list, cwd=cwd, check=check,
                                 stdout=stdout_dest, stderr=stderr_pipe,
                                 text=True, encoding='utf-8')
        if stderr_capture and process.stderr:
             stderr_lower = process.stderr.lower()
             if "error" in stderr_lower or "warning" in stderr_lower or process.returncode != 0 or len(process.stderr) < 500 :
                 print("--- STDERR ---")
                 print(process.stderr.strip())
                 print("--------------")
        print(f"Command finished {'successfully' if process.returncode == 0 else 'with exit code ' + str(process.returncode)}.")
        return process
    except FileNotFoundError as e:
        print(f"\nError: Command '{command_list[0]}' not found - {e}.")
        print(f"Please ensure '{command_list[0]}' is installed and in your system's PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nError executing command: {' '.join(command_list)}")
        print(f"Return code: {e.returncode}")
        if e.stderr: print(f"--- STDERR ---\n{e.stderr.strip()}\n--------------")
        if check: raise
        return e
    except Exception as e:
        print(f"\nAn unexpected error occurred running command: {e}")
        raise

def get_species_from_id_regex(sequence_id, compiled_regex, prefix_map):
    """Uses a pre-compiled regex to determine species name (prioritizes longer prefixes)."""
    match = compiled_regex.match(sequence_id)
    if match:
        matched_prefix = match.group(1)
        return prefix_map.get(matched_prefix)
    if sequence_id.startswith("ENSMUSP"): return "Mouse" # Specific fallback
    if sequence_id.startswith("ENSP"): return "Human" # Specific fallback
    return "UnknownSpecies"

def parse_gene_symbol(description_string):
    """Extracts gene symbol using regex, looking for 'gene_symbol:SYMBOL' pattern in FASTA description."""
    if not isinstance(description_string, str):
        return None
    # Regex to find "gene_symbol:" followed by non-whitespace characters
    match = re.search(r'gene_symbol:(\S+)', description_string, re.IGNORECASE)
    if match:
        symbol = match.group(1)
        # Remove potential trailing brackets/source info if attached directly
        symbol = symbol.split('[')[0]
        # Handle cases like 'gene_symbol:ATP6 ' (trailing space)
        symbol = symbol.strip()
        return symbol
    # --- Removed the less reliable description: fal
    return None # Return None if no specific 'gene_symbol:' tag found

# --- Main Workflow ---

print("--- Starting Phylogenetic Pipeline (Using Multiple Source FASTA Files) ---")

# 1. Validate Inputs and Tools
print("\n--- 1. Validating Inputs & Tools ---")
# Check executables, source FASTA, sequence dictionary (as before) ...
if shutil.which(mafft_executable) is None: print(f"Error: MAFFT executable '{mafft_executable}' not found."); sys.exit(1)
print(f"Found MAFFT: {shutil.which(mafft_executable)}")
if shutil.which(raxml_executable) is None: print(f"Error: RAxML executable '{raxml_executable}' not found."); sys.exit(1)
print(f"Found RAxML: {shutil.which(raxml_executable)}")
valid_source_files = []
for f_path in source_fasta_files:
    if not os.path.exists(f_path): print(f"Error: Source FASTA file not found: {f_path}."); sys.exit(1)
    valid_source_files.append(f_path)
print(f"Using source FASTA files: {', '.join(valid_source_files)}")
if not sequences_to_extract: print("Error: 'sequences_to_extract' is empty."); sys.exit(1)

# Create output directories (as before) ...
mafft_output_dir = os.path.join(base_output_dir, mafft_output_subdir)
raxml_output_dir = os.path.join(base_output_dir, raxml_output_subdir)
os.makedirs(mafft_output_dir, exist_ok=True)
os.makedirs(raxml_output_dir, exist_ok=True)
print(f"MAFFT output will be in: {os.path.abspath(mafft_output_dir)}")
print(f"RAxML output will be in: {os.path.abspath(raxml_output_dir)}")

# 2. Prepare Species Regex and Index Source Files (as before) ...
print("\n--- 2. Indexing Source Files & Preparing Species Regex ---")
try:
    sorted_prefixes = sorted(prefix_to_species_map.keys(), key=len, reverse=True)
    if not sorted_prefixes: raise ValueError("Prefix map is empty.")
    regex_pattern_str = r"^(" + "|".join(re.escape(p) for p in sorted_prefixes) + r")"
    species_regex = re.compile(regex_pattern_str)
    print(f"Compiled species regex pattern: {species_regex.pattern}")
except Exception as e: print(f"Error preparing species regex: {e}"); sys.exit(1)

source_indices = []
total_indexed = 0
for f_path in valid_source_files:
    print(f"Indexing: {os.path.basename(f_path)} ...")
    try:
        index = SeqIO.index(f_path, "fasta")
        count = len(index)
        if count == 0: print(f"  Warning: Source FASTA file '{f_path}' appears empty.")
        else: print(f"  Indexed {count} sequences."); source_indices.append(index); total_indexed += count
    except Exception as e: print(f"  Error indexing file {f_path}: {e}."); sys.exit(1)
if not source_indices: print("Error: Failed to index any source FASTA files."); sys.exit(1)
print(f"Total sequences indexed across {len(source_indices)} file(s): {total_indexed}")

# --- Removed BLAST Loading Step ---

# 3. Process Each Group: Extract, Format, Align, Build Tree (Renumbered)
print("\n--- 3. Processing Sequence Groups ---") # Renumbered Step
failed_groups = []
processed_group_count = 0

for group_name, id_list in sequences_to_extract.items():
    print(f"\n>>> Processing Group: {group_name} ({len(id_list)} IDs requested)")
    processed_group_count += 1
    sequences_for_group = []
    found_ids_this_group = set()
    missing_ids_in_group = []
    query_id_for_group = id_list[0] if id_list else None # Assumes first ID is query

    safe_group_name = re.sub(r'[\\/*?:"<>| ]', '_', group_name)
    temp_mafft_input_file = os.path.join(mafft_output_dir, f"{safe_group_name}_temp_mafft_input.fa")
    mafft_alignment_file = os.path.join(mafft_output_dir, f"{safe_group_name}.aln.fa")

    # a. Extract sequences and format headers
    print("  a. Extracting sequences and formatting headers...")
    for seq_id in id_list:
        record_found_in_any_file = None
        for index in source_indices:
            if seq_id in index:
                try: record_found_in_any_file = index[seq_id]; break
                except Exception as e: print(f"    Warning: Error retrieving '{seq_id}': {e}")

        if record_found_in_any_file:
            if seq_id not in found_ids_this_group:
                # --- MODIFIED HEADER FORMATTING (using FASTA description) ---
                original_id = record_found_in_any_file.id # Keep original ID base
                species = get_species_from_id_regex(original_id, species_regex, prefix_to_species_map)

                # Try to get gene symbol directly from FASTA description
                gene_symbol = parse_gene_symbol(record_found_in_any_file.description)

                # Fallback 1: If symbol not found AND this is the query ID, use group name
                if gene_symbol is None and original_id == query_id_for_group:
                     base_group_name = group_name.split('_')[0]
                     if base_group_name:
                          gene_symbol = base_group_name
                          # print(f"  DEBUG: Using group name '{gene_symbol}' as fallback symbol for query {original_id}")

                # Fallback 2: If still no symbol, use "UnknownGene"
                if gene_symbol is None:
                    gene_symbol = "UnknownGene"
                    print(f"    Warning: Could not parse gene symbol for '{original_id}' from its FASTA header. Using '{gene_symbol}'. Header: >{original_id} {record_found_in_any_file.description[:100]}...") # Print warning

                # Construct the new header
                clean_id = original_id.split('|')[0].split(' ')[0]
                new_header_id = f"{clean_id}|{species}|{gene_symbol}"
                # --- END MODIFICATION ---

                formatted_record = SeqRecord(record_found_in_any_file.seq, id=new_header_id, description="", name="")
                sequences_for_group.append(formatted_record)
                found_ids_this_group.add(seq_id)
        else:
            missing_ids_in_group.append(seq_id)

    print(f"    Found {len(found_ids_this_group)} out of {len(id_list)} requested unique IDs.")
    if missing_ids_in_group:
        print(f"    Error: Could not find required IDs: {', '.join(missing_ids_in_group)}")
        print(f"  Skipping group '{group_name}' due to missing sequences.")
        failed_groups.append(group_name + " (Missing IDs)")
        continue

    if len(sequences_for_group) < 2:
        print(f"  Skipping group '{group_name}': Fewer than 2 sequences found/retrieved.")
        failed_groups.append(group_name + " (<2 Seqs)")
        continue

    # b. Run MAFFT (as before) ...
    print(f"  b. Running MAFFT (Output: {os.path.basename(mafft_alignment_file)})...")
    mafft_success = False
    try:
        # Write sequences with the new headers to the temporary file
        with open(temp_mafft_input_file, "w") as temp_out: SeqIO.write(sequences_for_group, temp_out, "fasta")
        mafft_command = [mafft_executable, f"--thread", str(mafft_threads)] + mafft_options + [temp_mafft_input_file]
        with open(mafft_alignment_file, "w") as outfile: run_command(mafft_command, stdout_file=outfile, check=True)
        mafft_success = True
        print(f"    MAFFT alignment saved successfully.")
    except Exception as e:
        print(f"  Error during MAFFT alignment for group '{group_name}': {e}")
        if os.path.exists(mafft_alignment_file):
            try: os.remove(mafft_alignment_file); print(f"    Cleaned up incomplete MAFFT output.")
            except OSError as remove_err: print(f"    Warning: Could not remove MAFFT output file {mafft_alignment_file}: {remove_err}")
        failed_groups.append(group_name + " (MAFFT Error)")
    finally:
        # Ensure temporary file is removed even if MAFFT fails
        if os.path.exists(temp_mafft_input_file):
            try: os.remove(temp_mafft_input_file)
            except OSError as remove_err: print(f"    Warning: Could not remove temporary MAFFT input file {temp_mafft_input_file}: {remove_err}")

    # c. Run RAxML (as before) ...
    if not mafft_success: continue
    print(f"  c. Running RAxML v8 (Results in: {os.path.abspath(raxml_output_dir)})...")
    if len(sequences_for_group) < MIN_SEQS_FOR_RAXML:
        print(f"  Skipping RAxML for group '{group_name}': Only {len(sequences_for_group)} sequences (requires >= {MIN_SEQS_FOR_RAXML}).")
        continue
    raxml_output_label = safe_group_name
    final_tree_file = os.path.join(raxml_output_dir, f"RAxML_bipartitions.{raxml_output_label}")
    if os.path.exists(final_tree_file):
        print(f"    Output file {os.path.basename(final_tree_file)} already exists. Skipping RAxML run.")
        continue
    try:
        # RAxML uses the alignment file which now has the correct headers
        raxml_cmd = [
            raxml_executable, "-f", "a", "-m", raxml_model,
            "-s", os.path.abspath(mafft_alignment_file), # Input alignment file
            "-w", os.path.abspath(raxml_output_dir), "-n", raxml_output_label,
            "-T", str(raxml_threads), "-x", str(raxml_seed_bootstrap),
            "-p", str(raxml_seed_parsimony), "-#", str(raxml_bootstrap_reps)
        ]
        run_command(raxml_cmd, check=True)
        print(f"    RAxML analysis completed successfully for {group_name}.")
        print(f"    Key output tree: {os.path.basename(final_tree_file)}")
    except Exception as e:
        print(f"  Error during RAxML analysis for group '{group_name}': {e}")
        failed_groups.append(group_name + " (RAxML Error)")

# 4. Final Summary (Renumbered, logic remains the same)
print("\n--- 4. Workflow Summary ---") # Renumbered Step
total_groups = len(sequences_to_extract)
success_groups = total_groups - len(failed_groups)
print(f"Attempted processing for {total_groups} groups.")
print(f"Successfully completed the pipeline (MAFFT and RAxML where applicable) for {success_groups} groups.")
if failed_groups:
    failed_reasons = defaultdict(list)
    for group_info in failed_groups:
        parts = group_info.split(" (", 1)
        group = parts[0]
        reason = parts[1][:-1] if len(parts) > 1 else "Unknown Reason"
        failed_reasons[reason].append(group)
    print("Failed or skipped groups:")
    for reason, groups in failed_reasons.items(): print(f"  - Reason: {reason} -> Groups: {', '.join(groups)}")

print(f"\nPipeline finished.")
print(f"MAFFT alignments with modified headers (ID|Species|GeneSymbol) are in: {os.path.abspath(mafft_output_dir)}")
print(f"RAxML results (based on modified alignments) are in: {os.path.abspath(raxml_output_dir)}")

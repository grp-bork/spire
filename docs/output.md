# Output
The following outputs are expected for the respective steps in the pipeline

## Filtered Data

Human/Quality filtered Fastqs:
* `<id>/<id>.filtered.pair.1.fq.gz`
* `<id>/<id>.filtered.pair.2.fq.gz`
* `<id>/<id>.filtered.singles.fq.gz`
* `read_count_after_qc.txt`

## MAGs:

Assembly (megahit)
* `assemblies/<id>-assembled.fa.gz`:

Gene Calling (prodigal)
* `prodigal/<id>.psa_megahit.prodigal.faa.gz`
* `prodigal/<id>.psa_megahit.prodigal.fna.gz`
* `prodigal/<id>.psa_megahit.prodigal.gff.gz`

Calculate Depths (jgi_summarize_bam_contig_depths)
* `bins/<id>_aligned_to_<id>.depths`

Binning (metabat2)
* `bins/bins/<id>.psa_megahit.psb_metabat2.{000001}.fa.gz`

Per-Bin Gene Calling (seqtk)
* `per_bin_genecalls/<id>_per_bin_genecalls.tar.gz`

Assembly Stats (assembly-stats)
* `assembly_stats/<id>.assembly_stats.tsv`

Assembly and Bin Mash Sketching (mash)
* `assembly_mash_sketching/<id>-assembled.fa.gz.msh`
* `bin_mash_sketching/<id>/<id>.psa_megahit.psb_metabat2.{000001.}fa.msh`

## Annotation:

rRNA Detection (barrnap)
* `rrna/<id>-barrnap.arc.gff`
* `rrna/<id>-barrnap.bac.gff`
* `rrna/<id>-barrnap.euk.gff`
* `rrna/<id>-barrnap.mito.gff`

ARG Annotation (abricate, rgi)
* `abricate/<id>.psa_megahit.abricate.megares.tsv`
* `rgiv6/<id>/*`

Virulence Factor Annotation (abricate)
* `abricate/<id>.psa_megahit.abricate.vfdb.tsv`

sORFs Detection (macrel)
* `macrel/<id>.smorfs.faa`

Genome Quality Assessment (checkm2, gunc)
* `gunc/<id>.GUNC.maxCSS_level.tsv`
* `gunc/<id>.GUNC.maxCSS_level_gunc5.tsv`
* `gunc/.GUNC.all_levels/<id>.psa_megahit.psb_metabat2.{000001}.progenomes_2.1.all_levels.tsv`
* `checkm2/<id>.checkm2.tsv`

Functional Annotation (eggnog-mapper)
* `eggnog_mapper/<id>.emapper.annotations.gz`
* `eggnog_mapper/<id>.emapper.hits.gz`
* `eggnog_mapper/<id>.emapper.seed_orthologs.gz`

Taxonomic Classification (gtdbtk)
* `gtdbtk/<id>/gtdbtk.log`
* `gtdbtk/<id>/gtdbtk.warnings.log`
* `gtdbtk/<id>/identify/gtdbtk.ar53.markers_summary.tsv`
* `gtdbtk/<id>/identify/gtdbtk.bac120.markers_summary.tsv`
* `gtdbtk/<id>/identify/gtdbtk.failed_genomes.tsv`
* `gtdbtk/<id>/identify/gtdbtk.translation_table_summary.tsv`
* `gtdbtk/<id>/align/gtdbtk.bac120.filtered.tsv`
* `gtdbtk/<id>/align/gtdbtk.bac120.msa.fasta.gz`
* `gtdbtk/<id>/align/gtdbtk.bac120.user_msa.fasta.gz`
* `gtdbtk/<id>/classify/gtdbtk.bac120.classify.tree.2.tree`
* `gtdbtk/<id>/classify/gtdbtk.bac120.summary.tsv`
* `gtdbtk/<id>/classify/gtdbtk.bac120.tree.mapping.tsv`
* `gtdbtk/<id>/classify/gtdbtk.backbone.bac120.classify.tree`

# CADD-SV Snakemake to Nextflow DSL2 Conversion Mapping

## Workflow Overview

The CADD-SV workflow annotates SVs with ~30 different genomic features across 3 regions:
1. **Main SV span** (the actual SV coordinates)
2. **100bp upstream flank** (100bp before the SV)
3. **100bp downstream flank** (100bp after the SV)

Then combines all features into a matrix and applies Random Forest scoring.

## Phase 1: BED Preparation Rules → Nextflow Processes

| Snakemake Rule | Nextflow Process | Input | Output | Chr Prefix |
|----------------|------------------|-------|--------|------------|
| prep_chr1 | CADDSV_PREP_BED | input VCF/BED | wchr.bed | with chr |
| prep_chr2 | CADDSV_STRIP_CHR | wchr.bed | nochr.bed | without chr |
| prep_merg1 | CADDSV_MERGE_BED | wchr.bed, nochr.bed | merged beds | both |
| prep_chr_100bpup | CADDSV_FLANK_UP | wchr.bed, genome | 100bpup beds | both |
| prep_chr_100bpdown | CADDSV_FLANK_DOWN | wchr.bed, genome | 100bpdown beds | both |

## Phase 2: Annotation Rules → Nextflow Processes

### Conservation & Constraint Scores
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| cadd_PC_phylop | CADDSV_ANNO_CADD_PC | PhastCons/CADD_PC_PhyloP_scores.bed.gz | map max,sum |
| cadd2 | CADDSV_ANNO_CADD2 | CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz | count |
| gerp | CADDSV_ANNO_GERP | gerp/gerp_score2_hg38_MAM_90q.bed.gz | map max |
| gerp2 | CADDSV_ANNO_GERP2 | gerp/gerp_score2_hg38_MAM_90q.bed.gz | count |
| LINSIGHT | CADDSV_ANNO_LINSIGHT | linsight/LINSIGHT_hg38_sort.bed.gz | sum |
| ccr | CADDSV_ANNO_CCR | ccr/ccrs.all.bed.gz | map max |
| MPC | CADDSV_ANNO_MPC | MPC/transcript_constraints_hg38liftover.bg.gz | map mean |
| ultraconserved | CADDSV_ANNO_ULTRACONS | ultraconserved/ultraconserved_hg38_muliz120M_sort.bed | coverage |

### Gene Model Annotations
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| gene_model_tmp | CADDSV_GENE_TMP | ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz | tabix + gtf2bed |
| gene_model | CADDSV_GENE_MODEL | gm_tmp.bed | coverage per feature type |
| gene_model_dist | CADDSV_GENE_DIST | ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed | closest -d |
| gene_names | CADDSV_GENE_NAMES | ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed | map distinct |
| pli | CADDSV_ANNO_PLI | gnomad/pli_exac.csv | R script |

### Regulatory Features
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| CTCF | CADDSV_ANNO_CTCF | CTCF/Wangetal_hg38.bed | coverage |
| GC | CADDSV_ANNO_GC | GC/gc5Base.bedGraph.gz | map mean |
| chromHMM_MAX | CADDSV_ANNO_CHROMHMM | chromhmm/chromHMM_GRCh38.bg.gz | map max (25 cols) |
| encode (x13) | CADDSV_ANNO_ENCODE | encode/{mark}/*_merged_90quant.bed.gz | map max,sum |
| encode2 | CADDSV_MERGE_ENCODE | all encode outputs | paste |
| RemapTF | CADDSV_ANNO_REMAP | ReMap/ReMap2_overlapTF_hg38.bg.gz | map mean |
| Fantom5_counts | CADDSV_ANNO_FANTOM5 | fantom5/F5.hg38.enhancers_sort.bed | coverage counts |
| HI | CADDSV_ANNO_HI | DDD_HI/hg38_HI_Predictions_version3_sort.bed | map max |
| deepc | CADDSV_ANNO_DEEPC | deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz | max value |
| EP | CADDSV_ANNO_EP | enhancer-promoter-links/sorted_encode.bed | closest + R |

### 3D Genome / HiC Features
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| fire (x5 cell lines) | CADDSV_ANNO_FIRE | FIRE/fire_{cellline}.bed | map max,min |
| fire2 | CADDSV_MERGE_FIRE | all fire outputs | paste |
| HIC_hESC | CADDSV_ANNO_HIC_HESC | hic/hESC/combined/sorted.total.combined.domain | closest + R |
| HIC_encode (x4) | CADDSV_ANNO_HIC_ENCODE | Encode-HIC/{cell}/sorted_{tad}.bed.gz | closest + R |
| mergetad | CADDSV_MERGE_HIC | all HIC outputs | paste |
| microsynteny | CADDSV_ANNO_SYNTENY | synteny/microsynteny.bed | closest + R |
| genomegitar1 (x19) | CADDSV_ANNO_GENOMEGITAR | genomegitar/{id}/DI_sort.bed | map max,min |
| genomegitar2 | CADDSV_MERGE_GENOMEGITAR | all genomegitar outputs | paste |
| genomegitar3 | CADDSV_GENOMEGITAR_MINMAX | DI.bed | awk min/max |

## Phase 3: Matrix Assembly Rules → Nextflow Processes

| Snakemake Rule | Nextflow Process | Description |
|----------------|------------------|-------------|
| complete_script | CADDSV_ASSEMBLE_MATRIX | Combine all annotations for main span |
| complete_script_100bpup | CADDSV_ASSEMBLE_MATRIX_UP | Combine all annotations for upstream flank |
| complete_script_100bpdown | CADDSV_ASSEMBLE_MATRIX_DOWN | Combine all annotations for downstream flank |

## Phase 4: Scoring Rules → Nextflow Processes

| Snakemake Rule | Nextflow Process | Description |
|----------------|------------------|-------------|
| scoring | CADDSV_SCORE | R script with RF models |
| sort | CADDSV_SORT_OUTPUT | Sort and add header |

## Conda Environments

| Environment File | Purpose | Key Tools |
|------------------|---------|-----------|
| envs/SV.yml | Main annotation | bedtools, tabix, gtf2bed, R |
| envs/prepBED.yml | BED preparation | bedtools |

## Required Annotation Files (from annotations/ directory)

Total: ~50+ annotation files across subdirectories
EOFcat > rule_mapping.md << 'EOF'
# CADD-SV Snakemake to Nextflow DSL2 Conversion Mapping

## Workflow Overview

The CADD-SV workflow annotates SVs with ~30 different genomic features across 3 regions:
1. **Main SV span** (the actual SV coordinates)
2. **100bp upstream flank** (100bp before the SV)
3. **100bp downstream flank** (100bp after the SV)

Then combines all features into a matrix and applies Random Forest scoring.

## Phase 1: BED Preparation Rules → Nextflow Processes

| Snakemake Rule | Nextflow Process | Input | Output | Chr Prefix |
|----------------|------------------|-------|--------|------------|
| prep_chr1 | CADDSV_PREP_BED | input VCF/BED | wchr.bed | with chr |
| prep_chr2 | CADDSV_STRIP_CHR | wchr.bed | nochr.bed | without chr |
| prep_merg1 | CADDSV_MERGE_BED | wchr.bed, nochr.bed | merged beds | both |
| prep_chr_100bpup | CADDSV_FLANK_UP | wchr.bed, genome | 100bpup beds | both |
| prep_chr_100bpdown | CADDSV_FLANK_DOWN | wchr.bed, genome | 100bpdown beds | both |

## Phase 2: Annotation Rules → Nextflow Processes

### Conservation & Constraint Scores
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| cadd_PC_phylop | CADDSV_ANNO_CADD_PC | PhastCons/CADD_PC_PhyloP_scores.bed.gz | map max,sum |
| cadd2 | CADDSV_ANNO_CADD2 | CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz | count |
| gerp | CADDSV_ANNO_GERP | gerp/gerp_score2_hg38_MAM_90q.bed.gz | map max |
| gerp2 | CADDSV_ANNO_GERP2 | gerp/gerp_score2_hg38_MAM_90q.bed.gz | count |
| LINSIGHT | CADDSV_ANNO_LINSIGHT | linsight/LINSIGHT_hg38_sort.bed.gz | sum |
| ccr | CADDSV_ANNO_CCR | ccr/ccrs.all.bed.gz | map max |
| MPC | CADDSV_ANNO_MPC | MPC/transcript_constraints_hg38liftover.bg.gz | map mean |
| ultraconserved | CADDSV_ANNO_ULTRACONS | ultraconserved/ultraconserved_hg38_muliz120M_sort.bed | coverage |

### Gene Model Annotations
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| gene_model_tmp | CADDSV_GENE_TMP | ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz | tabix + gtf2bed |
| gene_model | CADDSV_GENE_MODEL | gm_tmp.bed | coverage per feature type |
| gene_model_dist | CADDSV_GENE_DIST | ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed | closest -d |
| gene_names | CADDSV_GENE_NAMES | ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed | map distinct |
| pli | CADDSV_ANNO_PLI | gnomad/pli_exac.csv | R script |

### Regulatory Features
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| CTCF | CADDSV_ANNO_CTCF | CTCF/Wangetal_hg38.bed | coverage |
| GC | CADDSV_ANNO_GC | GC/gc5Base.bedGraph.gz | map mean |
| chromHMM_MAX | CADDSV_ANNO_CHROMHMM | chromhmm/chromHMM_GRCh38.bg.gz | map max (25 cols) |
| encode (x13) | CADDSV_ANNO_ENCODE | encode/{mark}/*_merged_90quant.bed.gz | map max,sum |
| encode2 | CADDSV_MERGE_ENCODE | all encode outputs | paste |
| RemapTF | CADDSV_ANNO_REMAP | ReMap/ReMap2_overlapTF_hg38.bg.gz | map mean |
| Fantom5_counts | CADDSV_ANNO_FANTOM5 | fantom5/F5.hg38.enhancers_sort.bed | coverage counts |
| HI | CADDSV_ANNO_HI | DDD_HI/hg38_HI_Predictions_version3_sort.bed | map max |
| deepc | CADDSV_ANNO_DEEPC | deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz | max value |
| EP | CADDSV_ANNO_EP | enhancer-promoter-links/sorted_encode.bed | closest + R |

### 3D Genome / HiC Features
| Snakemake Rule | Nextflow Process | Annotation File | Operation |
|----------------|------------------|-----------------|-----------|
| fire (x5 cell lines) | CADDSV_ANNO_FIRE | FIRE/fire_{cellline}.bed | map max,min |
| fire2 | CADDSV_MERGE_FIRE | all fire outputs | paste |
| HIC_hESC | CADDSV_ANNO_HIC_HESC | hic/hESC/combined/sorted.total.combined.domain | closest + R |
| HIC_encode (x4) | CADDSV_ANNO_HIC_ENCODE | Encode-HIC/{cell}/sorted_{tad}.bed.gz | closest + R |
| mergetad | CADDSV_MERGE_HIC | all HIC outputs | paste |
| microsynteny | CADDSV_ANNO_SYNTENY | synteny/microsynteny.bed | closest + R |
| genomegitar1 (x19) | CADDSV_ANNO_GENOMEGITAR | genomegitar/{id}/DI_sort.bed | map max,min |
| genomegitar2 | CADDSV_MERGE_GENOMEGITAR | all genomegitar outputs | paste |
| genomegitar3 | CADDSV_GENOMEGITAR_MINMAX | DI.bed | awk min/max |

## Phase 3: Matrix Assembly Rules → Nextflow Processes

| Snakemake Rule | Nextflow Process | Description |
|----------------|------------------|-------------|
| complete_script | CADDSV_ASSEMBLE_MATRIX | Combine all annotations for main span |
| complete_script_100bpup | CADDSV_ASSEMBLE_MATRIX_UP | Combine all annotations for upstream flank |
| complete_script_100bpdown | CADDSV_ASSEMBLE_MATRIX_DOWN | Combine all annotations for downstream flank |

## Phase 4: Scoring Rules → Nextflow Processes

| Snakemake Rule | Nextflow Process | Description |
|----------------|------------------|-------------|
| scoring | CADDSV_SCORE | R script with RF models |
| sort | CADDSV_SORT_OUTPUT | Sort and add header |

## Conda Environments

| Environment File | Purpose | Key Tools |
|------------------|---------|-----------|
| envs/SV.yml | Main annotation | bedtools, tabix, gtf2bed, R |
| envs/prepBED.yml | BED preparation | bedtools |

## Required Annotation Files (from annotations/ directory)

Total: ~50+ annotation files across subdirectories

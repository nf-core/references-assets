# Properly dealing with vcf files

## Assets with vcf files

These are the assets that have vcf files.

- assets/igenomes/Bos_taurus/Ensembl/UMD3.1.yml
- assets/igenomes/Homo_sapiens/Ensembl/GRCh37.yml
- assets/igenomes/Homo_sapiens/GATK/b37.yml
- assets/igenomes/Homo_sapiens/GATK/GRCh37decoy.yml
- assets/igenomes/Homo_sapiens/GATK/GRCh38.yml
- assets/igenomes/Homo_sapiens/GATK/hg19.yml
- assets/igenomes/Homo_sapiens/GATK/hg38.yml
- assets/igenomes/Homo_sapiens/UCSC/hg19.yml
- assets/igenomes/Mus_musculus/Ensembl/GRCm38.yml
- assets/igenomes/Rattus_norvegicus/Ensembl/Rnor_5.0.yml
- assets/igenomes/Rattus_norvegicus/Ensembl/Rnor_6.0.yml
- assets/igenomes/Sus_scrofa/Ensembl/Sscrofa10.2.yml

Currently VCFs are only stored in these assets as a single vcf entry for the same genome key.

My idea, listing all the VCFs and categorize them into group, that are ready to use for the pipelines.

## VCF files

I'm categorizing the vcf files into the following groups:

### dbsnp

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.dbsnp.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_144.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/dbsnp_138.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.dbsnp.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/dbsnp_138.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/dbsnp_144.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/dbsnp_146.hg38.vcf.gz

## known_indels

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.known_indels.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_phase1.indels.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz
s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz

## known_snps

s3://ngi-igenomes/igenomes/Bos_taurus/Ensembl/UMD3.1/Annotation/Variation/Bos_taurus.vcf
s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Variation/Homo_sapiens.vcf
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_omni2.5.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase3_v4_20130502.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_omni2.5.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_phase3_v4_20130502.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/1000G_omni2.5.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp142.vcf
s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Variation/Mus_musculus.vcf
s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
s3://ngi-igenomes/igenomes/Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Variation/Rattus_norvegicus.vcf
s3://ngi-igenomes/igenomes/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Variation/Rattus_norvegicus.vcf
s3://ngi-igenomes/igenomes/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Variation/Sus_scrofa.vcf

## germline_resource

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/af-only-gnomad.raw.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh37.PASS.AC.AF.only.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz

## hapmap

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/hapmap_3.3.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/hapmap_3.3_grch38_pop_stratified_af.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/hapmap_3.3.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/hapmap_3.3_b37_pop_stratified_af.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/hapmap_3.3.hg19.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/hapmap_3.3.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz

## pon

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GermlineResource/dummy_PON.gnomAD.GRCh37.WGS.AF.GT.01.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/dummy_PON.gnomAD.GRCh38.WGS.AF.GT.01.vcf.gz

## Other

No idea yet on how to categorize these files.

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/CEUTrio.HiSeq.WGS.b37.NA12878.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/CEUTrio.HiSeq.WGS.b37.bestPractices.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/NA12878.knowledgebase.snapshot.20131119.b37.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.variantEvalGoldStandard.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/NISTIntegratedCalls.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.variantEvalGoldStandard.vcf.gz
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/NISTIntegratedCalls.hg38.vcf.gz

## Dropping these files

s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf

They are vcf, but there is already a vcf.gz version of them.

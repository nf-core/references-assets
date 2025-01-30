#!/usr/bin/env bash

# Download or regenerate manifest file

# Command line flags
while getopts ":adr" opt; do
    case $opt in
        a)
            NEW_MANIFEST=true
            ;;
        d)
            DOWNLOAD_MANIFEST=true
            ;;
        r)
            REGENERATE_MANIFEST=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1;
            ;;
    esac
done

if [[ ${DOWNLOAD_MANIFEST} ]]; then
    rm -f manifest.txt
    echo "Downloading manifest"
    wget https://raw.githubusercontent.com/ewels/AWS-iGenomes/refs/heads/master/ngi-igenomes_file_manifest.txt -O manifest.txt
fi

if [[ ${REGENERATE_MANIFEST} ]]; then
    rm -f manifest.txt
    echo "Regenerating manifest"
    # cf https://github.com/ewels/AWS-iGenomes/pull/22
    aws s3 --no-sign-request ls --recursive s3://ngi-igenomes/igenomes/ | cut -d "/" -f 2- > tmp
    for i in `cat tmp`; do
        if [[ ! $i =~ /$ ]]; then
            echo s3://ngi-igenomes/igenomes/$i >> tmp_manifest
        fi
    done
    mv tmp_manifest manifest.txt

    rm tmp
fi

if [[ ${NEW_MANIFEST} ]]; then
    rm -f manifest.txt
    echo "Regenerating NEW manifest"
    # cf https://github.com/ewels/AWS-iGenomes/pull/22

    aws s3 --profile igenomes ls --recursive s3://nf-core-references-scratch/genomes/ | grep -v "pipeline_info" | cut -d "/" -f 2- > tmp
    for i in `cat tmp`; do
        if [[ ! $i =~ /$ ]]; then
            echo s3://nf-core-references-scratch/genomes/$i >> tmp_manifest
        fi
    done
    mv tmp_manifest manifest.txt

    rm tmp
fi

total_files=$(wc -l manifest.txt | cut -d " " -f 1)

echo "Number of files in manifest: $total_files"

cp manifest.txt leftover_manifest.txt

# Remove existing assets
rm -rf igenomes/

# Generate base info in species/genome/build.yml

# All source fasta.fai
# ALL fai are coming from a fasta file of the same name
# Hence I use it to generate fasta + fai (and catch with that the fasta that are not following the gemome.fa name scheme)

cat manifest.txt | grep "\.fai" | grep -v "Bowtie2Index" | grep -v "fai\.gz" > tmp_fai.txt

echo "Populating assets for fasta and fai"

for i in `cat tmp_fai.txt`;
do
    # Create a list of fasta files
    echo "${i::-4}" >> tmp_fasta.txt

    # Grab metadata
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    mkdir -p igenomes/${genome}

    echo "- genome: \"${build}\"" > igenomes/${genome}/${build_name}
    echo "  fasta: \"${i::-4}\"" >> igenomes/${genome}/${build_name}
    echo "  source: \"${genome}\"" >> igenomes/${genome}/${build_name}
    echo "  species: \"${species}\"" >> igenomes/${genome}/${build_name}
    echo "  fasta_fai: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source README
cat manifest.txt | grep "README" | grep -v "Archives" | grep -v "beagle" | grep -v "plink" | grep -v "PhiX\/Illumina\/RTA\/Annotation\/README\.txt" > tmp_readme.txt

echo "Populating assets for README"

for i in `cat tmp_readme.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  readme: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source gtf (removing the onces coming from gencode)
cat manifest.txt | grep "\.gtf" | grep -v "gtf\." | grep -v "STARIndex" | grep -v "Genes\.gencode" > tmp_gtf.txt

echo "Populating assets for GTF"

for i in `cat tmp_gtf.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  gtf: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source fasta.dict
cat manifest.txt | grep "\.dict" | grep -v "dict\.gz" | grep -v "dict\.old" > tmp_dict.txt

echo "Populating assets for fasta.dict"

for i in `cat tmp_dict.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  fasta_dict: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source genes.bed
cat manifest.txt | grep "genes\.bed" > tmp_bed.txt

echo "Populating assets for genes.bed"

for i in `cat tmp_bed.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  genes_bed: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source BowtieIndex
cat manifest.txt | grep "BowtieIndex" | grep -v "MDSBowtieIndex" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bowtie.txt

echo "Populating assets for BowtieIndex"

for i in `cat tmp_bowtie.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  bowtie1_index: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source Bowtie2Index
cat manifest.txt | grep "Bowtie2Index" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bowtie2.txt

echo "Populating assets for Bowtie2Index"

for i in `cat tmp_bowtie2.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  bowtie2_index: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source BWAIndex (we have version0.6.0, version0.5.x, and no version specified)
cat manifest.txt | grep "BWAIndex" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bwaindex.txt

echo "Populating assets for BWAIndex"

for i in `cat tmp_bwaindex.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    # Remove existing bwamem1_index if present
    # So that we only keep the latest version
    sed -i '\|bwamem1_index|d' igenomes/${genome}/${build_name}
    echo "  bwamem1_index: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source BWAmem2mem
cat manifest.txt | grep "BWAmem2Index" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bwamem2mem.txt

echo "Populating assets for BWAmem2Index"

for i in `cat tmp_bwamem2mem.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  bwamem2_index: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source Dragmap
cat manifest.txt | grep "dragmap" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_dragmap.txt

echo "Populating assets for DragmapHashtable"

for i in `cat tmp_dragmap.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  dragmap_hashtable: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source BismarkIndex
cat manifest.txt | grep "BismarkIndex\/genome\.fa" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bismark.txt

echo "Populating assets for BismarkIndex"

for i in `cat tmp_bismark.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  bismark_index: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source star Index
cat manifest.txt | grep "STARIndex" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_star.txt

echo "Populating assets for STARIndex"

for i in `cat tmp_star.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  star_index: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source Chromosomes fasta
cat manifest.txt | grep "Chromosomes" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_chromosomes.txt

echo "Populating assets for Chromosomes fasta"

for i in `cat tmp_chromosomes.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  chromosomes_fasta: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source AbundantSequences fasta
cat manifest.txt | grep "AbundantSequences" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_abundantsequences.txt

echo "Populating assets for AbundantSequences fasta"

for i in `cat tmp_abundantsequences.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  abundantsequences_fasta: \"${i}/\"" >> igenomes/${genome}/${build_name}
done

# All source refFlat (removing the ones coming from gencode)
cat manifest.txt | grep "refFlat\.txt" | grep -v "\.gz\.bak" | grep -v "Genes\.gencode" > tmp_refflat.txt

echo "Populating assets for refFlat"

for i in `cat tmp_refflat.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  genes_refflat: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source refgene
cat manifest.txt | grep "refGene\.txt" | grep -v "Archives" | grep -v "\.gz\.bak" > tmp_refgene.txt

echo "Populating assets for refgene"

for i in `cat tmp_refgene.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  genes_refgene: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source ChromInfo.txt
cat manifest.txt | grep "ChromInfo\.txt" | grep -v "Archives" > tmp_chrominfo.txt

echo "Populating assets for ChromInfo.txt"

for i in `cat tmp_chrominfo.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  chrom_info: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source GenomeSize.xml (removing the old ones)
cat manifest.txt | grep "GenomeSize\.xml" | grep -v "\.old" > tmp_GenomeSize.txt

echo "Populating assets for GenomeSize.xml"

for i in `cat tmp_GenomeSize.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  genome_size_xml: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source hairpin.fa
cat manifest.txt | grep "SmallRNA\/hairpin\.fa" > tmp_hairpin.txt

echo "Populating assets for hairpin.fa"

for i in `cat tmp_hairpin.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  hairpin_fasta: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source mature.fa
cat manifest.txt | grep "SmallRNA\/mature\.fa" > tmp_mature.txt

echo "Populating assets for mature.fa"

for i in `cat tmp_mature.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    echo "  mature_fasta: \"${i}\"" >> igenomes/${genome}/${build_name}
done

# All source vcf
cat manifest.txt | grep "\.vcf" | grep -v "\.idx" | grep -v "\.tbi" | grep -v "\.md5" > tmp_vcf.txt

echo "Populating assets for vcf"

a_dbsnp=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.dbsnp.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_144.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/dbsnp_138.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.dbsnp.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/dbsnp_138.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/dbsnp_144.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/dbsnp_146.hg38.vcf.gz")

a_known_indels=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.known_indels.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_phase1.indels.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz"
    "s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz")

a_known_snps=("s3://ngi-igenomes/igenomes/Bos_taurus/Ensembl/UMD3.1/Annotation/Variation/Bos_taurus.vcf"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Variation/Homo_sapiens.vcf"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_omni2.5.b37.vcf.gz"    
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.b37.vcf.gz"    
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase3_v4_20130502.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_omni2.5.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/1000G_phase3_v4_20130502.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/1000G_omni2.5.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp142.vcf"
    "s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Variation/Mus_musculus.vcf"
    "s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
    "s3://ngi-igenomes/igenomes/Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Variation/Rattus_norvegicus.vcf"
    "s3://ngi-igenomes/igenomes/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Variation/Rattus_norvegicus.vcf"
    "s3://ngi-igenomes/igenomes/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Variation/Sus_scrofa.vcf")

a_germline_resource=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/af-only-gnomad.raw.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh37.PASS.AC.AF.only.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz")

a_hapmap=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/hapmap_3.3.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/hapmap_3.3_grch38_pop_stratified_af.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/hapmap_3.3.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/hapmap_3.3_b37_pop_stratified_af.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/hapmap_3.3.hg19.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/hapmap_3.3.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz")

a_pon=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GermlineResource/dummy_PON.gnomAD.GRCh37.WGS.AF.GT.01.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/dummy_PON.gnomAD.GRCh38.WGS.AF.GT.01.vcf.gz")

a_other=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/CEUTrio.HiSeq.WGS.b37.NA12878.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/CEUTrio.HiSeq.WGS.b37.bestPractices.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/b37/NA12878.knowledgebase.snapshot.20131119.b37.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.variantEvalGoldStandard.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/NISTIntegratedCalls.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/Homo_sapiens_assembly38.variantEvalGoldStandard.vcf.gz"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/beta/NISTIntegratedCalls.hg38.vcf.gz")

a_skip=("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf"
    "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf")

for i in `cat tmp_vcf.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    source_vcf=$(echo $i | cut -d "/" -f 9)
    build_name="${build}.yml"

    if [[ $build =~ ^build ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $build =~ ^EF ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^Enterobacteriophage_lambda ]]; then
        build_name="${species}_${build}.yml"
    elif [[ $species =~ ^PhiX ]]; then
        build_name="${species}_${build}.yml"
    fi

    vcf_name=$(basename $i)
    if [ -z "${source_vcf}" ]; then
        source_vcf="unknown"
    elif [ "${source_vcf}" == "GermlineResource" ]; then
        source_vcf="GATKBundle"
    elif [ "${source_vcf}" == "Variation" ]; then
        source_vcf="unknown"
    elif [ "${source_vcf}" == "${vcf_name}" ]; then
        source_vcf="unknown"
    fi

    if [[ ${a_skip[@]} =~ $i ]]
    then
        echo "${i}"
        echo "skipping"
        continue
    elif [[ ${a_dbsnp[@]} =~ $i ]]
    then
        vcf_category="dbsnp"
    elif [[ ${a_known_indels[@]} =~ $i ]]
    then
        vcf_category="known_indels"
    elif [[ ${a_known_snps[@]} =~ $i ]]
    then
        vcf_category="known_snps"
    elif [[ ${a_germline_resource[@]} =~ $i ]]
    then
        vcf_category="germline_resource"
    elif [[ ${a_hapmap[@]} =~ $i ]]
    then
        vcf_category="hapmap"
    elif [[ ${a_pon[@]} =~ $i ]]
    then
        vcf_category="pon"
    elif [[ ${a_other[@]} =~ $i ]]
    then
        vcf_category="other"
    else
        echo "${i}"
        echo "not categorised"
        continue
    fi

    if ! [ -f igenomes/${genome}/${build_name} ]; then
        echo "- genome: \"${build}\"" > igenomes/${genome}/${build_name}
        echo "  source: \"${genome}\"" >> igenomes/${genome}/${build_name}
        echo "  species: \"${species}\"" >> igenomes/${genome}/${build_name}
        echo "  vcf:" >> igenomes/${genome}/${build_name}
    fi

    ! grep -q vcf igenomes/${genome}/${build_name} && echo "  vcf:" >> igenomes/${genome}/${build_name}

    echo "    ${vcf_category}:" >> igenomes/${genome}/${build_name}
    echo "      vcf: \"${i}\"" >> igenomes/${genome}/${build_name}

    grep -q ${i}.idx manifest.txt && echo "      vcf_idx: \"${i}.idx\"" >> igenomes/${genome}/${build_name}
    grep -q ${i}.tbi manifest.txt && echo "      vcf_tbi: \"${i}.tbi\"" >> igenomes/${genome}/${build_name}

    if [[ ${source_vcf} = "unknown" ]]; then
        continue
    fi

    echo "      vcf_source: \"${source_vcf}\"" >> igenomes/${genome}/${build_name}
done

echo "Fixing /GATK/GRCh37.yml name to GATK/GRCh37decoy.yml"

#  GATK/GRCh37.yml should actually be GATK/GRCh37decoy.yml
mv igenomes/GATK/GRCh37.yml igenomes/GATK/GRCh37decoy.yml

echo "Deleting assets from manifest"

# Delete all the lines that are used to generate the assets
for i in `cat tmp_*.txt`;
do
    sed -i "\|${i}|d" leftover_manifest.txt
done

leftover_files=$(wc -l leftover_manifest.txt | cut -d " " -f 1)
echo "Number of files in leftover manifest: $leftover_files"

let removed_files=$total_files-$leftover_files
echo "Number of files removed: $removed_files"

echo "Deleting tmp files"

rm -rf tmp_*

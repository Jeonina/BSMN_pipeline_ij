
#!/usr/bin/env bash

set -Ee -o pipefail

SECONDS=0

PIPEHOME=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
RESDIR=$PIPEHOME/resources

mkdir -p $RESDIR $RESDIR/hg19 $RESDIR/hg38

# Activate conda
eval "$(conda shell.bash hook)"
conda activate --no-stack bp

########################################
# b37 (GRCh37)
########################################
if [[ ! -f $RESDIR/hg19/human_g1k_v37_decoy.fasta.fai ]]; then
    echo "Downloading b37 reference..."

    cd $RESDIR/hg19

    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.dict.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta.amb
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta.ann
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta.bwt
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta.pac
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta.sa
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/dbsnp_138.b37.vcf.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_omni2.5.b37.vcf.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/hapmap_3.3.b37.vcf.gz
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz

    cd -

    gunzip -f $RESDIR/hg19/human_g1k_v37_decoy.*.gz

    echo "Indexing b37..."
    for V in $RESDIR/hg19/*.b37.vcf.gz; do
        gunzip -f $V
        bgzip -f -@ 4 ${V/.gz/}
        tabix -p vcf -f $V
    done

    samtools faidx $RESDIR/hg19/human_g1k_v37_decoy.fasta
    echo "b37 Done."
fi

########################################
# hg38 (FTP 사용)
########################################
if [[ ! -f $RESDIR/hg38/dbsnp_146.hg38.vcf.gz ]]; then
    echo "Downloading hg38 reference (FTP)..."

    cd $RESDIR/hg38

    # reference
    lftp -c 'set ftp:web-mode true; pget -n 10 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict'

    # index 생성
    bwa index Homo_sapiens_assembly38.fasta
    samtools faidx Homo_sapiens_assembly38.fasta

    # known sites
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz'

    # dbsnp
    lftp -c 'set ftp:web-mode true; pget -n 10 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz'

    cd -
fi

########################################
# hg19 UCSC
########################################
if [[ ! -f $RESDIR/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz ]]; then
    echo "Downloading hg19 UCSC..."

    cd $RESDIR/hg19

    wget -c https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.dict
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.fasta
    wget -c https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.fasta.fai

    lftp -c 'set ftp:web-mode true; pget -n 10 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 4 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz'
    lftp -c 'set ftp:web-mode true; pget -n 10 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz'

    echo "Indexing hg19..."
    for V in *.hg19*.vcf.gz; do
        gunzip -f $V
        bgzip -f -@ 4 ${V/.gz/}
        tabix -p vcf -f $V
    done

    cd -
fi

########################################
# liftOver
########################################
if [[ ! -f $RESDIR/hg19ToHg38.over.chain.gz ]]; then
    wget -c -P $RESDIR http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi

########################################
# Finish
########################################
conda deactivate

elapsed=$SECONDS
printf "\n>> Total $(($elapsed / 3600))h $(($elapsed % 3600 / 60))m $(($elapsed % 60))s\n"

#Modified from https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

git clone https://github.com/bulik/ldsc.git

cd /~/ldsc

conda env create --file environment.yml

source activate ldsc
#Once the above has completed, you can run:
./ldsc.py -h
./munge_sumstats.py -h

#munge trait 1
## add column names for each input
./munge_sumstats.py \
--sumstats /~/trait1.txt --out trait1 \
--N \
--snp \
--p  \
--signed-sumstats  \
--a1 \
--a2 \
--frq \
--chunksize 500000 \
--merge-alleles /~/w_hm3.snplist

#munge trait 2
## add column names for each input
./munge_sumstats.py \
--sumstats /~/trait2.txt --out trait2 \
--N \
--snp \
--p  \
--signed-sumstats  \
--a1 \
--a2 \
--frq \
--chunksize 500000 \
--merge-alleles /~/w_hm3.snplist

# LD Score Regression 
./ldsc.py \
--rg trait1.sumstats.gz,trait2.sumstats.gz \
--ref-ld-chr /~/eur_w_ld_chr/ \
--w-ld-chr /~/eur_w_ld_chr/ \
--out trait1_trait2

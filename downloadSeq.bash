mkdir -p saccharomyces_rnaSeq/{data,results/{map,filter,counts},scripts}

## latest version saccharomyces cerevisiae 
## https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/
cd saccharomyces_rnaSeq/data
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=19jQgKqXJkrSWjtHxbrt0KnYxB38FDwHJ' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=19jQgKqXJkrSWjtHxbrt0KnYxB38FDwHJ" -O ERR2929688_sample0.1.fastq.gz && rm -rf /tmp/cookies.txt

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1f7Yk6LU0mgVYyLgft1vxkubK4wUax_HX' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1f7Yk6LU0mgVYyLgft1vxkubK4wUax_HX" -O ERR2929687_sample0.1.fastq.gz && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1yfygGyea0zTFXN02m1vvhlp5rzDMNuER' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1yfygGyea0zTFXN02m1vvhlp5rzDMNuER" -O ERR2929686_sample0.1.fastq.gz && rm -rf /tmp/cookies.txt

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1Dv7ea7QJVZjygSXivpOfk6cPQfRqiNIc' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1Dv7ea7QJVZjygSXivpOfk6cPQfRqiNIc" -O ERR2929685_sample0.1.fastq.gz && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1yWJFe-QlrJeZX6deAlsiAhx7Skrc0l4H' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1yWJFe-QlrJeZX6deAlsiAhx7Skrc0l4H" -O ERR2929684_sample0.1.fastq.gz && rm -rf /tmp/cookies.txt

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1SwksI4CM5ffHSTBxyZ33ESI5rPiYtogF' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1SwksI4CM5ffHSTBxyZ33ESI5rPiYtogF" -O ERR2929683_sample0.1.fastq.gz && rm -rf /tmp/cookies.txt

cd ../results/counts/

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1rdSBRmKbVE1vVrS_0gKn_Jtjxm9XjSie' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1rdSBRmKbVE1vVrS_0gKn_Jtjxm9XjSie" -O final_counts.txt && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1ogX24QXctuEz_7qStqm2b2oYvqtLQ3ie' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1ogX24QXctuEz_7qStqm2b2oYvqtLQ3ie" -O final_counts.txt.summary && rm -rf /tmp/cookies.txt

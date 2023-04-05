### RNAseq analysis performed for the initial version of report
#### De novo transcriptome assembly (trinityrnaseq-v2.15.1)
        Trinity --seqType fq --max_memory 50G --SS_lib_type F --CPU 10 --no_normalize_reads --single E100000803_L01_36.fq.gz --output trinity_outQ61
        Trinity --seqType fq --max_memory 50G --SS_lib_type F --CPU 10 --no_normalize_reads --single E100000803_L01_40.fq.gz --output trinity_outQ64
        Trinity --seqType fq --max_memory 50G --SS_lib_type F --CPU 6 --no_normalize_reads --single E100000804_L01_16.fq.gz --output trinity_out
        Trinity --seqType fq --max_memory 50G --SS_lib_type F --CPU 6 --no_normalize_reads --single E100000804_L01_17.fq.gz --output trinity_outQ69
        Trinity --seqType fq --max_memory 50G --SS_lib_type F --CPU 6 --no_normalize_reads --single E100000804_L01_20.fq.gz --output trinity_outQ70

#### Transcripts BLAST search (ncbi-blast-2.13.0+) with in house database of MARKET species genomes

        makeblastdb -in MSExtraWGS.fasta -dbtype nucl ##### MSExtraWGS.fasta contains whole genome seqeunces in fasta format
        
        blastn -query trinity_outQ61.Trinity.fasta -db MSExtraWGS.fasta -outfmt '6 qseqid sseqid pident evalue score bitscore length qstart qend sstart send stitle' -max_target_seqs 2 -out trinity_outQ61.Trinity_blastn_MSTwoHits.txt -num_threads 30
        awk '$7>=100 {print}' trinity_outQ61.Trinity_blastn_MSTwoHits.txt >trinity_outQ61.Trinity_blastn_MSTwoHits100.txt
        
        blastn -query trinity_outQ64.Trinity.fasta -db MSExtraWGS.fasta -outfmt '6 qseqid sseqid pident evalue score bitscore length qstart qend sstart send stitle' -max_target_seqs 2 -out trinity_outQ64.Trinity_blastn_MSTwoHits.txt -num_threads 30
        awk '$7>=100 {print}' trinity_outQ64.Trinity_blastn_MSTwoHits.txt >trinity_outQ64.Trinity_blastn_MSTwoHits100.txt
        
        blastn -query trinity_outQ69.Trinity.fasta -db MSExtraWGS.fasta -outfmt '6 qseqid sseqid pident evalue score bitscore length qstart qend sstart send stitle' -max_target_seqs 2 -out trinity_outQ69.Trinity_blastn_MSTwoHits.txt -num_threads 30
        awk '$7>=100 {print}' trinity_outQ69.Trinity_blastn_MSTwoHits.txt >trinity_outQ69.Trinity_blastn_MSTwoHits100.txt
        
        blastn -query trinity_outQ70.Trinity.fasta -db MSExtraWGS.fasta -outfmt '6 qseqid sseqid pident evalue score bitscore length qstart qend sstart send stitle' -max_target_seqs 2 -out trinity_outQ70.Trinity_blastn_MSTwoHits.txt -num_threads 30
        awk '$7>=100 {print}' trinity_outQ70.Trinity_blastn_MSTwoHits.txt >trinity_outQ70.Trinity_blastn_MSTwoHits100.txt
        
        blastn -query trinity_out.Trinity.fasta -db MSExtraWGS.fasta -outfmt '6 qseqid sseqid pident evalue score bitscore length qstart qend sstart send stitle' -max_target_seqs 2 -out trinity_out.TrinityQ68_blastn_MSTwoHits.txt -num_threads 30
        awk '$7>=100 {print}' trinity_out.TrinityQ68_blastn_MSTwoHits.txt >trinity_out.TrinityQ68_blastn_MSTwoHits100.txt
        
#### Annotation (TransDecoder-TransDecoder-v5.7.0) Q61 transcripts aligned with Nyctereutes procyonoides genome

        cut -f1 trinity_outQ61.Trinity_blastn_MSTwoHits100.txt | sort -u >BestHitOrtholgs/outQ61BH
        perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' outQ61BH trinity_outQ61.Trinity.fasta >outQ61BH.fasta
        
        TransDecoder.LongOrfs -t outQ61BH.fasta --output_dir ./
        blastp -query longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
        hmmsearch --cpu 8 -E 1e-10 --domtblout pfam.domtblout Pfam-A.hmm longest_orfs.pep
        TransDecoder.Predict -t outQ61BH.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --output_dir ./
        
        cp outQ61BH.fasta.transdecoder.pep outQ61BH.fasta.transdecoder.fasta
        mkdir Ortholog_Mapping
        cp outQ61BH.fasta.transdecoder.fasta ../Ortholog_Mapping/
        
        cd Ortholog_Mapping
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz
        gunzip GCF_000001405.40_GRCh38.p14_protein.faa.gz
        
        mkdir OrthoMapQ61
        mv outQ61BH.fasta.transdecoder.fasta OrthoMapQ61/
        cp GCF_000001405.40_GRCh38.p14_protein.faa.gz OrthoMapQ61/
        cd OrthoMapQ61/
        sed -E '/^>/s/( +|\t).*//' GCF_000001405.40_GRCh38.p14_protein.faa >HumanProte.fasta
        rm -rf GCF_000001405.40_GRCh38.p14_protein.faa
        
        vim outQ61BH.fasta.transdecoder.fasta  #### clean the fasta file
        :%s/ TRINITY_.*/ /g
        :%s/\*/ /g
        :%s/\s\+$//

#### Ortholog mapping (OrthoFinder version_2.5.4)

       ./orthofinder -f /Full_PATH/OrthoMapQ61
        
        
        
    

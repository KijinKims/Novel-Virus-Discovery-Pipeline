#!/bin/sh

###############################################################
#							      #
#	Novel Virus Discovery Pipeline			      #
#							      #
#	This pipeline is for identifying chanses of 	      #
#	existing viruses in biological sample using	      #
#	NGS data.					      #
#							      #
#	Any theoretical questions can be solved with	      #
#	My presentation file. (17/12/11 presented)	      #
#							      #
#	And instruction for this is provided in enclosed      #
#	Manual file. Please refer to it.(Written in Korean)   #
#							      #
#	If you have any question, please email me.	      #
#	e-mail: skkujin@gmail.com			      #
#							      #
#	Ki-jin Kim 17/12/29				      #
#							      #
###############################################################

#Paths of Software/DB
HOSTDBPATH=/home/jwsong/Desktop/Pipeline/HOSTDB
TRIMMOMATICPATH=/home/jwsong/miniconda3/envs/pipeline/share/trimmomatic/trimmomatic.jar
CHECKPOINTPATH=/home/jwsong/Desktop/Pipeline
BLASTDBPATH=/home/jwsong/Desktop/Pipeline/BLASTDB
BACTERIADBPATH=/home/jwsong/Desktop/Pipeline/Bacteria_DB
ADAPTORPATH=/home/jwsong/miniconda3/envs/pipeline/share/trimmomatic/adapters
FAQCSPATH=/home/jwsong/miniconda3/envs/pipeline/bin/FaQCs
BOWTIEPATH=/home/jwsong/miniconda3/envs/pipeline/bin/bowtie2
SAMTOOLSPATH=/home/jwsong/miniconda3/envs/pipeline/bin/samtools
BEDTOOLSPATH=/home/jwsong/miniconda3/envs/pipeline/bin/bedtools
SPADESPATH=/home/jwsong/miniconda3/envs/pipeline/bin/spades.py
BLASTNPATH=/home/jwsong/miniconda3/envs/pipeline/bin/blastn
TBLASTXPATH=/home/jwsong/miniconda3/envs/pipeline/bin/tblastx

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

#User Input
echo "=====Please type name====="
read name

echo "=====Please type DB name.(Human/Mouse_Rat/Bat/Tick/Cro/Sorex/Penguin)====="
read DB

#Print start time
printf "\t*****Starting time: ${TIMESTAMP}*****\n"

#Clipping adaptor
printf "\n\n=====Trimming adaptor with trimmomatic=====\n\n"

Step=$(grep "CUTADAPT" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "CUTADAPT" ]
        then
                cd Adaptor_Clipping
		java -jar ${TRIMMOMATICPATH} PE -phred33 ../$1 ../$2 ${name}_1.fq unpaired.fq ${name}_2.fq unpaired.fq ILLUMINACLIP:${ADAPTORPATH}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		cd ..
fi

#Remove PhiX sequence
printf "\n\n=====PhiX removal with Bowtie2=====\n\n"

Step=$(grep "BOWTIE2" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BOWTIE2" ]
        then
                cd Host_removal
                ${BOWTIEPATH} -x ${HOSTDBPATH}/PhiX -1 ../Adaptor_Clipping/${name}_1.fq -2 ../Adaptor_Clipping/${name}_2.fq -S ${name}_mapped_and_unmapped.sam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi

fi

printf "\n\n=====Convert sam files into bam files=====\n\n"

Step=$(grep "SAMTOBAM" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMTOBAM" ]
        then
                ${SAMTOOLSPATH} view -bS ${name}_mapped_and_unmapped.sam > ${name}_mapped_and_unmapped.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\n=====Filter required unmapped reads=====\n\n"

Step=$(grep "SAMFLAGS" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMFLAGS" ]
        then
                ${SAMTOOLSPATH} view -b -f 12 -F 256 ${name}_mapped_and_unmapped.bam > ${name}_bothEndsUnmapped.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\n=====Sort of the bam files=====\n\n"

Step=$(grep "BAMSORT" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMSORT" ]
        then
                ${SAMTOOLSPATH} sort -n ${name}_bothEndsUnmapped.bam -o ${name}_bothEndsUnmapped_sorted.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\n=====Convert bam files into fastq files=====\n\n"

Step=$(grep "BAMTOFASTQ" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMTOFASTQ" ]
        then
                ${BEDTOOLSPATH} bamtofastq -i ${name}_bothEndsUnmapped_sorted.bam -fq ${name}_phiX_removed_1.fastq -fq2 ${name}_phiX_removed_2.fastq
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                rm ${name}_mapped_and_unmapped.sam
                rm ${name}_mapped_and_unmapped.bam
                rm ${name}_bothEndsUnmapped.bam
                rm ${name}_bothEndsUnmapped_sorted.bam
                cd ..
fi

#Remove host sequence
printf "\n\n=====Host removal with Bowtie2=====\n\n"

Step=$(grep "BOWTIE2" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BOWTIE2" ]
        then
                cd Host_removal
                ${BOWTIEPATH} -x ${HOSTDBPATH}/${DB} -1 ${name}_phiX_removed_1.fastq -2 ${name}_phiX_removed_2.fastq -S ${name}_mapped_and_unmapped.sam
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		rm ${name}_phiX_removed_1.fastq
		rm ${name}_phiX_removed_2.fastq
fi

printf "\n\n=====Convert sam files into bam files=====\n\n"

Step=$(grep "SAMTOBAM" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMTOBAM" ]
        then
		${SAMTOOLSPATH} view -bS ${name}_mapped_and_unmapped.sam > ${name}_mapped_and_unmapped.bam
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\n=====Filter required unmapped reads=====\n\n"

Step=$(grep "SAMFLAGS" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMFLAGS" ]
        then
		${SAMTOOLSPATH} view -b -f 12 -F 256 ${name}_mapped_and_unmapped.bam > ${name}_bothEndsUnmapped.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\n=====Sort of the bam files=====\n\n"

Step=$(grep "BAMSORT" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMSORT" ]
        then
		${SAMTOOLSPATH} sort -n ${name}_bothEndsUnmapped.bam -o ${name}_bothEndsUnmapped_sorted.bam
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\n=====Convert bam files into fastq files=====\n\n"

Step=$(grep "BAMTOFASTQ" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMTOFASTQ" ]
        then
		${BEDTOOLSPATH} bamtofastq -i ${name}_bothEndsUnmapped_sorted.bam -fq ${name}_host_removed_1.fastq -fq2 ${name}_host_removed_2.fastq
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		rm ${name}_mapped_and_unmapped.sam
		rm ${name}_mapped_and_unmapped.bam
		rm ${name}_bothEndsUnmapped.bam
		rm ${name}_bothEndsUnmapped_sorted.bam
		cd ..
fi

#Filter reads of ungood quality
printf "\n\n=====Filtering with FaQCs=====\n\n"

Step=$(grep "FILTER" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "FILTER" ]
        then
		cd Quality_filter
		${FAQCSPATH} -p ../Host_removal/${name}_host_removed_1.fastq ../Host_removal/${name}_host_removed_2.fastq -d ${name}
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		cd ..
fi

printf "\n\n=====De novo assembly with SPAdes=====\n\n"

Step=$(grep "DENOVO" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "DENOVO" ]
        then
		cd De_novo_Assembly
		${SPADESPATH} --pe1-1 ../Quality_filter/${name}/QC.1.trimmed.fastq --pe1-2 ../Quality_filter/${name}/QC.2.trimmed.fastq -o ${name}
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		cd ${name}
                cp contigs.fasta ../${name}.fas
                cd ../..
fi

#BLAST search
#Blastn(somewhat similar)
printf "\n\n=====Blast Search with Blastn=====\n\n"

Step=$(grep "BLASTN" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BLASTN" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_blastn.csv" 
                ${BLASTNPATH} -query ../De_novo_Assembly/${name}.fas -task blastn -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_blastn.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

#Megablast(high similiar)
printf "\n\n=====Blast Search with Megablast (high similarity)=====\n\n"

Step=$(grep "MEGABLAST" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "MEGABLAST" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_megablast.csv"
                ${BLASTNPATH} -query ../De_novo_Assembly/${name}.fas -task megablast -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_megablast.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

#Bacteria screening with megablast
printf "\n\n=====Bacteria Blast Search with Megablast (high similarity)=====\n\n"

Step=$(grep "BACMEGABLAST" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BACMEGABLAST" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_bacteria_megablast.csv"
                ${BLASTNPATH} -query ../De_novo_Assembly/${name}.fas -task megablast -db ${BACTERIADBPATH}/Bacteria -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_bacteria_megablast.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

printf "\n\n=====Blast Search with DCmegablast (considering discontinuity)=====\n\n"

#DCmegablast(ignore 3rd of triplet code)
Step=$(grep "DCMEGABLAST" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "DCMEGABLAST" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_dcmegablast.csv"
                ${BLASTNPATH} -query ../De_novo_Assembly/${name}.fas -task dc-megablast -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_dcmegablast.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

#Blastx(amino acid)
printf "\n\n=====Blast Search with Blastx (amino acid sequence)=====\n\n"

Step=$(grep "BLASTX" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BLASTX" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_blastx.csv"
                ${TBLASTXPATH} -query ../De_novo_Assembly/${name}.fas -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_blastx.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

printf "\t*****End time: ${TIMESTAMP}*****\n"

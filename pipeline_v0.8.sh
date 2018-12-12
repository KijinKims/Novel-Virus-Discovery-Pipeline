#!/bin/sh

HOSTDBPATH=/home/molecularvirology/Desktop/Kijin/Pipeline/HOSTDB
TRIMMOMATICPATH=/home/molecularvirology/Documents/Trimmomatic-0.36/trimmomatic-0.36.jar
RAYPATH=/home/molecularvirology/Documents/Ray-2.3.1/ray-build/Ray
CHECKPOINTPATH=/home/molecularvirology/Desktop/Kijin/Pipeline
BLASTDBPATH=/home/molecularvirology/Desktop/Kijin/Pipeline/BLASTDB
ADAPTORPATH=/home/molecularvirology/Documents/Trimmomatic-0.36/adapters
FAQCSPATH=/home/molecularvirology/Documents/FaQCs-master

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

echo "Please type name."
read name

echo "Please type DB name.(Human/Mouse_Rat/Bat/Tick/Cro/Sorex/Penguin)"
read DB

logfile=${name}_log

printf "\tStarting time: ${TIMESTAMP}\n"

printf "\n\nTrimming adaptor with trimmomatic\n\n"

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

printf "\n\nPhiX removal with Bowtie2\n\n"

Step=$(grep "BOWTIE2" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BOWTIE2" ]
        then
                cd Host_removal
                bowtie2 -x ${HOSTDBPATH}/PhiX -1 ../Adaptor_Clipping/${name}_1.fq -2 ../Adaptor_Clipping/${name}_2.fq -S ${name}_mapped_and_unmapped.sam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi

fi

printf "\n\nConvert sam files into bam files\n\n"

Step=$(grep "SAMTOBAM" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMTOBAM" ]
        then
                samtools view -bS ${name}_mapped_and_unmapped.sam > ${name}_mapped_and_unmapped.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\nFilter required unmapped reads\n\n"

Step=$(grep "SAMFLAGS" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMFLAGS" ]
        then
                samtools view -b -f 12 -F 256 ${name}_mapped_and_unmapped.bam > ${name}_bothEndsUnmapped.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\nSort of the bam files\n\n"

Step=$(grep "BAMSORT" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMSORT" ]
        then
                samtools sort -n ${name}_bothEndsUnmapped.bam ${name}_bothEndsUnmapped_sorted
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\nConvert bam files into fastq files\n\n"

Step=$(grep "BAMTOFASTQ" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMTOFASTQ" ]
        then
                bedtools bamtofastq -i ${name}_bothEndsUnmapped_sorted.bam -fq ${name}_phiX_removed_1.fastq -fq2 ${name}_phiX_removed_2.fastq
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


printf "\n\nHost removal with Bowtie2\n\n"

Step=$(grep "BOWTIE2" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BOWTIE2" ]
        then
                cd Host_removal
                bowtie2 -x ${HOSTDBPATH}/${DB} -1 ${name}_phiX_removed_1.fastq -2 ${name}_phiX_removed_2.fastq -S ${name}_mapped_and_unmapped.sam
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		rm ${name}_phiX_removed_1.fastq
		rm ${name}_phiX_removed_2.fastq
fi

printf "\n\nConvert sam files into bam files\n\n"

Step=$(grep "SAMTOBAM" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMTOBAM" ]
        then
		samtools view -bS ${name}_mapped_and_unmapped.sam > ${name}_mapped_and_unmapped.bam
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\nFilter required unmapped reads\n\n"

Step=$(grep "SAMFLAGS" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "SAMFLAGS" ]
        then
		samtools view -b -f 12 -F 256 ${name}_mapped_and_unmapped.bam > ${name}_bothEndsUnmapped.bam
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\nSort of the bam files\n\n"

Step=$(grep "BAMSORT" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMSORT" ]
        then
		samtools sort -n ${name}_bothEndsUnmapped.bam ${name}_bothEndsUnmapped_sorted
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
fi

printf "\n\nConvert bam files into fastq files\n\n"

Step=$(grep "BAMTOFASTQ" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BAMTOFASTQ" ]
        then
		bedtools bamtofastq -i ${name}_bothEndsUnmapped_sorted.bam -fq ${name}_host_removed_1.fastq -fq2 ${name}_host_removed_2.fastq
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

printf "\n\nFiltering with FaQCs\n\n"

Step=$(grep "FILTER" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "FILTER" ]
        then
		cd Quality_filter
		perl ${FAQCSPATH}/FaQCs.pl -p ../Host_removal/${name}_host_removed_1.fastq ../Host_removal/${name}_host_removed_2.fastq -d ${name}
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		cd ..
fi

printf "\n\nDe novo assembly with SPAdes\n\n"

Step=$(grep "DENOVO" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "DENOVO" ]
        then
		cd De_novo_Assembly
		spades.py --pe1-1 ../Quality_filter/${name}/QC.1.trimmed.fastq --pe1-2 ../Quality_filter/${name}/QC.2.trimmed.fastq -o ${name}
		if [ $? -ne 0 ]
                        then
                                exit 1
                fi
		cd ${name}
                cp contigs.fasta ../${name}.fas
                cd ../..
fi

printf "\n\nBlast Search with Blastn\n\n"

Step=$(grep "BLASTN" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BLASTN" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_blastn.csv" 
                blastn -query ../De_novo_Assembly/${name}.fas -task blastn -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_blastn.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

printf "\n\nBlast Search with Megablast (high similarity)\n\n"

Step=$(grep "MEGABLAST" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "MEGABLAST" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_megablast.csv"
                blastn -query ../De_novo_Assembly/${name}.fas -task megablast -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_megablast.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

printf "\n\nBlast Search with DCmegablast (considering discontinuity)\n\n"

Step=$(grep "DCMEGABLAST" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "DCMEGABLAST" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_dcmegablast.csv"
                blastn -query ../De_novo_Assembly/${name}.fas -task dc-megablast -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_dcmegablast.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi

printf "\n\nBlast Search with Blastx (amino acid sequence)\n\n"

Step=$(grep "BLASTX" ${CHECKPOINTPATH}/checkpoint.txt)
if [ "${Step}" != "BLASTX" ]
        then
                cd Screening
                printf " Contig No., Contig Len., Subject Accession, E-Value, Bitscore, Align Len., Subject Title, Status, Identity, Start, End\n" | tee -a "${name}_blastx.csv"
                tblastx -query ../De_novo_Assembly/${name}.fas -db ${BLASTDBPATH}/Virus -evalue 1.0e-5 -outfmt "10 qseqid qlen sacc evalue bitscore length stitle pident sstart send" -max_target_seqs 10 | tee -a "${name}_blastx.csv"
                if [ $? -ne 0 ]
                        then
                                exit 1
                fi
                cd ..
fi


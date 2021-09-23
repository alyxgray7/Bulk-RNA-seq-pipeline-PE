rule convertFQ2FA:
    input:
        r1 = "samples/star/{sample}_bam/Unmapped.out.mate1",
        r2 = "samples/star/{sample}_bam/Unmapped.out.mate2"
    output:
        "samples/star/{sample}_bam/{sample}_unmapped.fa"
    message:
        """--- Converting Unmapped FASTQ to FASTA."""
    shell:
        """
        cat {input.r1} | sed -n '1~4s/^@/>/p;2~4p' > {output} &&
        cat {input.r2} | sed -n '1~4s/^@/>/p;2~4p' >> {output}
        """


rule countUniqUnmappedSeqs:
    input:
        "samples/star/{sample}_bam/{sample}_unmapped.fa"
    output:
        "data/unmappedSeqs/{sample}_overRepseqCount.txt"
    message:
        """--- Counting unique unmapped sequences."""
    shell:
        """
        cat {input} | grep -v "^>" | sort | uniq -c | sort -nr > {output}
        """

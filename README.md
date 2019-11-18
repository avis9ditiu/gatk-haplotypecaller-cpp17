## Introduction
This Hardware-friendly version of GATK4 HaplotypeCaller.

## Getting started

	git clone https://github.com/avis9ditiu/gatk-haplotypecaller-cpp17.git
	cd gatk-haplotypecaller-cpp17
    mkdir build
    cd build
    cmake ..
    make
    ./gatk -I ../input/chrM.sam -R ../reference/chrM.fa -O ../output/chrM.vcf 

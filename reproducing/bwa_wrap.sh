#!/bin/bash

set -eu

USAGE="$(cat <<-END
Align reads using BWA.
    -b <path>,   --bwa <path>
        BWA executable [default: bwa].
    -i <prefix>, --input <prefix>
        Input prefix. There must be files <prefix>{1,2}.fq[.gz]
    -o <prefix>, --output <prefix>
        Output BAM/CRAM filename.
    -f <fasta>,  --fasta-ref <fasta>
        Fasta reference.
    -s <string>, --sample <string>
        Sample name [default: none].
    -t <int>,    --threads <int>.
        Number of threads [default: 4].
    -- [args]
        Arguments supplied to BWA MEM.
END
)"

bwa=bwa
bwa_args=""
threads=4
sample=""

while (( "$#" )); do
    case "$1" in
        -b|--bwa)
            bwa="$2"
            shift 2
            ;;
        -i|--input)
            in_prefix="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            fasta="$2"
            shift 2
            ;;
        -s|--sample)
            sample="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -h|--help)
            echo "${USAGE}"
            exit 0
            ;;
        --)
            shift
            bwa_args="$@"
            break
            ;;
        *)
            echo "Error: Unexpected argument $1" >&2
            exit 1
            ;;
    esac
done

compress=false
if [[ -f ${in_prefix}1.fq.gz ]]; then
    echo "Decompressing reads"
    pigz -d -p ${threads} ${in_prefix}{1,2}.fq.gz
    compress=true
fi

echo "Aligning reads"
if [ -z ${sample} ]; then
    bwa mem ${fasta} -t ${threads} $bwa_args ${in_prefix}{1,2}.fq > ${output}.unsort.sam
else
    bwa mem ${fasta} -R "@RG\tID:${sample}\tSM:${sample}" -t ${threads} $bwa_args \
        ${in_prefix}{1,2}.fq > ${output}.unsort.sam
fi

echo "Sorting alignments"
samtools sort -@ ${threads} --reference ${fasta} -o ${output} ${output}.unsort.sam

echo "Indexing alignments"
samtools index -@ ${threads} ${output}
rm ${output}.unsort.sam

if [[ ${compress} = true ]]; then
    echo "Compressing reads back"
    pigz -p ${threads} ${in_prefix}{1,2}.fq
fi

echo "Success"

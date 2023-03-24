#!/bin/bash

set -e

wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Evaluate diploid genotypes.
    -b <file>,  --benchmark <file>
        Benchmark BED file.
    -p <file>,  --parascopy <file>
        Parascopy variants[_pooled].bed.gz file.
    -f <file>,  --fasta-ref <file>
        Fasta reference.
    -r <file>,  --regions <file>
        BED file with all duplication regions.
    -R <file>,  --ref-cn <file>
        BED file with the reference copy numbers.
    -m <str>,   --mode <str>
        Benchmarking mode.
    -o <dir>,   --output <dir>
        Output directory.
    -n <int> <int>,   --cn-bounds <int> <int>
        Minimal and maximal reference aggr. copy number values.
        See get_comp_regions for default.
        Does not affect "all" mode.
END
)"

command="$0 $*"
bounds_arg=""

while (( "$#" )); do
    case "$1" in
        -b|--benchmark)
            benchmark="$2"
            shift 2
            ;;
        -p|--parascopy)
            parascopy="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            genome="$2"
            shift 2
            ;;
        -r|--regions)
            regions="$2"
            shift 2
            ;;
        -R|--ref-cn)
            ref_cn="$2"
            shift 2
            ;;
        -m|--mode)
            mode="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -n|--cn-bounds)
            bounds_arg="--cn-bounds $2 $3"
            shift 3
            ;;
        -h|--help)
            echo "${USAGE}"
            exit 0
            ;;
        *)
            >&2 echo "Error: Unexpected argument $1"
            exit 1
            ;;
    esac
done

if [[ -z ${regions} ]] || [[ ! -f ${regions} ]]; then
    >&2 echo "Error: Regions BED file is not provided or does not exist!"
    exit 1
elif [[ -z ${ref_cn} ]] || [[ ! -f ${ref_cn} ]]; then
    >&2 echo "Error: Reference copy number BED file is not provided or does not exist!"
    exit 1
fi

set -u

if [[ "$(basename ${parascopy})" = *pooled* ]]; then
    pooled_arg="--pooled"
else
    pooled_arg=""
fi

if [[ ${mode} = "unset" ]]; then
    >&2 echo "Error: mode is not provided!"
    exit 1
elif [[ ${mode} == non-* ]]; then
    mode=${mode#non-}
    bedtools_cmd=subtract
else
    bedtools_cmd=intersect
fi

if [[ ${mode} = all ]]; then
    ln -s "$(readlink -f ${regions})" ${output}/parascopy.bed
else
    if [[ ${mode} = anycn ]]; then
        filtering="--qual 0 --filters any --cn any"
    elif [[ ${mode} = nonref ]]; then
        filtering="--qual 0 --filters any --cn non-ref"
    elif [[ ${mode} = unfilt_refcn ]]; then
        filtering="--qual 0 --filters any --cn only-ref"
    elif [[ ${mode} = ref_pscn ]]; then
        filtering="--qual 0 --filters any --cn pscn-ref"
    elif [[ ${mode} = refcn ]]; then
        filtering="--qual 0 --filters only-pass --cn only-ref"
    elif [[ ${mode} = confid ]] || [[ ${mode} = highq_pass ]]; then
        filtering="--qual 20 --filters only-pass --cn only-ref"
    elif [[ ${mode} = highq_filter ]]; then
        filtering="--qual 20 --filters not-pass --cn only-ref"
    elif [[ ${mode} = lowq_pass ]]; then
        filtering="--qual -20 --filters only-pass --cn only-ref"
    elif [[ ${mode} = lowq_filter ]]; then
        filtering="--qual -20 --filters not-pass --cn only-ref"
    else
        >&2 echo "Error: unknown mode ${mode}"
        exit 1
    fi
    ${wdir}/get_comp_regions.py -i ${parascopy} -f ${genome} -o ${output}/parascopy.bed -r ${ref_cn} \
        $filtering $pooled_arg $bounds_arg
fi

bedtools intersect -a ${benchmark} -b ${regions} | \
    bedtools ${bedtools_cmd} -a - -b ${output}/parascopy.bed | \
    bedtools sort -i - -g ${genome}.fai | \
    bedtools merge -i - > ${output}/comparison.bed

#!/bin/bash

set -e

wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Evaluate diploid genotypes.
    -b <pref>,   --benchmark <prefix>
        Benchmark prefix. Must contain files <prefix>.vcf.gz and <prefix>.bed.
    -p <dir>,   --parascopy <dir>
        Parascopy directory.
    -s <str>,   --suffix <str>
        Parascopy suffix [default: ""].
        If not empty, output file will be named parascopy_<suffix>.
    -f <file>,  --fasta-ref <file>
        Fasta reference. Must contain <file>.sdf index.
    -c <dir>
        Path to directory with freebayes.vcf.gz and gatk.vcf.gz files [default: "calls"].
    -o <dir>,   --output <dir>
        Output directory. May contain comparison.bed file.
    -m <str>,   --mode <str>
        Select benchmarking regions based on the mode:
        - all:   use all input benchmarking regions,
        - refcn:  use benchmarking regions where Parascopy predicts psCN = 2 with filter = PASS,
        - unfilt_refcn: use benchmarking regions where Parascopy predicts psCN = 2 irrespective of any filters,
        - confid: use benchmarking regions where Parascopy predicts psCN = 2 with quality over 20,
        If the mode starts with "non-", use benchmarking regions minus "refcn", "unfilt_refcn" or "confid" regions.
    -r <file>,  --regions <file>
        Limit the analysis to the regions (BED file).
    -C <file>,  --ref-cn <file>
        Input BED file with reference copy number values.
        For format, see "parascopy call --assume-cn".
    -n <int> <int>,   --cn-bounds <int> <int>
        Limit analysis between two reference copy numbers.
        See get_comp_regions.py for default.
END
)"

calls="calls"
suffix=""
mode="unset"
command="$0 $*"

in_regions=""
ref_cn=""
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
        -s|--suffix)
            suffix="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            genome="$2"
            shift 2
            ;;
        -c|--calls)
            calls="$2"
            shift 2
            ;;
        -m|--mode)
            mode="$2"
            shift 2
            ;;
        -r|--regions)
            in_regions="$2"
            shift 2
            ;;
        -C|--ref-cn)
            ref_cn="$2"
            shift 2
            ;;
        -n|--cn-bounds)
            bounds_arg="--cn-bounds $2 $3"
            shift 3
            ;;
        -o|--output)
            output="$2"
            shift 2
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

if [[ -z ${benchmark} ]] || [[ ! -f ${benchmark}.vcf.gz ]] || [[ ! -f ${benchmark}.bed ]]; then
    >&2 echo "Error: Benchmark BED or VCF files (-b, --benchmark) are not provided OR do not exist!"
    exit 1
elif [[ -z ${parascopy} ]] || [[ ! -d ${parascopy} ]]; then
    >&2 echo "Error: Parascopy directory (-p, --parascopy) is not provided OR does not exist!"
    exit 1
elif [[ -z ${genome} ]] || [[ ! -f ${genome} ]]; then
    >&2 echo "Error: Fasta reference (-f, --fasta-ref) is not provided or does not exist!"
    exit 1
elif [[ ! -d ${genome}.sdf ]]; then
    >&2 echo "Error: Fasta reference (-f, --fasta-ref) does not have an .sdf index!"
    exit 1
# elif [[ ! -f ${calls}/freebayes.vcf.gz ]] || [[ ! -f ${calls}/gatk.vcf.gz ]]; then
#     >&2 echo "Error: Variant call files ${calls}/freebayes.vcf.gz OR ${calls}/gatk.vcf.gz do not exist!"
#     exit 1
elif [[ -z ${output} ]]; then
    >&2 echo "Error: Output directory (-o, --output) is not provided!"
    exit 1
fi

if [[ -z ${suffix} ]]; then
    par_name="parascopy"
else
    par_name="parascopy_${suffix}"
fi

set -u

mkdir -p ${output}
cmp_regions="${output}/comparison.bed"
filter_regions="${output}/filter.bed"
log="${output}/eval.log"
echo "Command: $command" > ${log}

if [[ ! -f ${cmp_regions} ]]; then
    >&2 echo "** Generating variant calling regions."
    ${wdir}/extract_regions.sh -b ${benchmark}.bed -f ${genome} -p ${parascopy}/variants.bed.gz \
        -m ${mode} -o ${output} -r ${in_regions} -R ${ref_cn} $bounds_arg &>> ${log}
    bedtools slop -i ${cmp_regions} -b 100 -g ${genome}.fai | \
        bedtools merge -i - -d 10 > ${filter_regions}
else
    >&2 echo "!! Variant calling regions already present."
fi

if ! (grep -vq '^#' ${cmp_regions}); then
    >&2 echo "!! Empty comparison regions, stopping !!"
    exit 0
fi

if [[ ! -f ${output}/${par_name}.vcf.gz.tbi ]]; then
    rm -f ${output}/${par_name}.vcf.gz
    rtg vcffilter -i ${parascopy}/variants.vcf.gz -o ${output}/${par_name}.vcf.gz \
        --include-bed=${filter_regions} -j ${wdir}/ploidy_gq0.js &>> ${log}
else
    >&2 echo "!! Skipping output Parascopy VCF file."
fi

if [[ ! -f ${output}/eval-${par_name}/done ]]; then
    >&2 echo "** Running Parascopy vcfeval."
    rm -rf ${output}/eval-${par_name}
    rtg vcfeval -b ${benchmark}.vcf.gz -c ${output}/${par_name}.vcf.gz \
        -t ${genome}.sdf -e ${cmp_regions} -o ${output}/eval-${par_name} &>> ${log}
else
    >&2 echo "!! Skipping Parascopy vcfeval."
fi

for f in ${calls}/*.vcf.gz; do
    caller="$(basename ${f} .vcf.gz)"
    if [[ -f ${output}/eval-${caller}/done ]]; then
        >&2 echo "!! Skipping ${caller}."
    else
        >&2 echo "** Running vcfeval for ${caller}."
        rm -rf ${output}/eval-${caller}
        rtg vcfeval -b ${benchmark}.vcf.gz -c ${calls}/${caller}.vcf.gz \
            -t ${genome}.sdf -e ${cmp_regions} -o ${output}/eval-${caller} &>> ${log}
    fi
done

printf "\n"
${wdir}/write_summary.py -t 20 -T ${output}/eval-*

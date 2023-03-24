#!/bin/bash

set -e

wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Simulate reads using ART Illumina.
    -f <prefix>,   --fasta-ref <prefix>
        Prefix to the reference fasta file.
        Files <prefix>1.fa and <prefix>2.fa must exist.
    -o <dir>,   --output <dir>
        Output directory.
    -@ <int>,   --threads <int>
        Number of threads.
    -s <int>,   --seed <int>
        Seeds (<seed> + 1) to (<seed> + 2 * <threads>) will be used.
        If not set, selects a random seed.
    -c <int>,   --coverage <int>
        Read coverage.
    -R,   --skip-renaming
        Do not combine and rename reads.
    -- <args>
        Arguments passed to art_illumina.
        Must not contain -c, -f.
END
)"

command="$0 $*"
combine_reads=true

while (( "$#" )); do
    case "$1" in
        -f|--fasta-ref)
            fasta_prefix="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -@|--threads)
            threads="$2"
            shift 2
            ;;
        -s|--seed)
            seed="$2"
            shift 2
            ;;
        -c|--coverage)
            coverage="$2"
            shift 2
            ;;
        -R|--skip-renaming)
            combine_reads=false
            shift
            ;;
        --)
            shift
            args="$*"
            break
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

if [[ -z ${fasta_prefix} ]] || [[ ! -f ${fasta_prefix}1.fa ]] || [[ ! -f ${fasta_prefix}2.fa ]]; then
    >&2 echo "Error: file ${fasta_prefix}1.fa or ${fasta_prefix}2.fa does not exist!"
    exit 1
elif [[ -z ${output} ]]; then
    >&2 echo "Error: Output directory is not set!"
    exit 1
elif [[ -z ${threads} ]] || [[ ${threads} -lt 1 ]]; then
    >&2 echo "Error: Number of threads is not set or too low!"
    exit 1
elif [[ -z ${coverage} ]] || [[ ${coverage} -lt 1 ]]; then
    >&2 echo "Error: Coverage is not set or too low!"
    exit 1
elif [[ -z ${args} ]]; then
    >&2 echo "Error: ART arguments are not set!"
    exit 1
fi

if [[ -z ${seed} ]]; then
    seed=${RANDOM}
fi

set -u

mkdir -p ${output}
log=${output}/simul.log

echo "Command: $command" | tee ${log}
export cov1="$(bc <<< "scale=4; ${coverage}/${threads}/2")"
echo "    Running ${threads} processes, each with coverage ${cov1}" | tee -a ${log}
echo "    Using random seeds $((seed + 1)) to $((seed + 2 * threads))" | tee -a ${log}

export fasta_prefix
export seed
export threads
export output
export args

simulate() {
    set -u

    i="$1"
    if [[ ${i} -le ${threads} ]]; then
        hap=1
        j=${i}
    else
        hap=2
        j=$((i - threads))
    fi

    out_prefix="$(printf "%s/%03d.hap%d." ${output} ${j} ${hap})"
    (set -x;
        art_illumina -i ${fasta_prefix}${hap}.fa --rndSeed $((seed + i)) -o ${out_prefix} \
            -f ${cov1} $args &> ${out_prefix}log
    )
}

export -f simulate

echo "Starting simulation" | tee -a ${log}
seq 1 $((threads * 2)) | xargs -P${threads} -i sh -c 'simulate {}' |& tee -a ${log}

if [[ ${combine_reads} = true ]]; then
    echo "Combine and rename reads" | tee -a ${log}
    seq 1 2 | xargs -P2 -i sh -c "${wdir}/rename_fq.awk ${output}/*.{}.fq > reads{}.fq"
else
    echo "Do not combine and rename reads" | tee -a ${log}
fi
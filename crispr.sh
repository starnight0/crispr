#!/bin/bash

function usage
{
    echo "Usage: crispr.sh ";
    echo "  -i [FASTA FILE]";
    echo "  -L [Target Length]";
    echo "  -o [OUT PREFIX]";
    echo "  -p [PAM seq]";
    echo "  -bt [bowtie index]";
    echo "  --dry : only print but not run.";
}

dry=0;
ref="";
while [ "$1" != "" ]; do
    case $1 in 
        -i | --in-fasta )  shift
            fasta=$1
            ;;
        -p | --pam )       shift
            pam=$1
            ;;
        --pam-pos ) shift
            pam_pos=$1
            ;;
        -L | --target-length ) shift
            tlen=$1
            ;;
        -o | --out-prefix ) shift
            out_prefix=$1
            ;;
        --bt | --bowtie-index ) shift
            index=$1
            ;;
        --ref ) shift
            ref=$1
            ;;
        --dry ) dry=1
            ;;
        -h | --help )  usage; exit 1
            ;;
    esac
    shift
done

if [ $fasta = "" ]; then
    echo "in_fasta_file is needed, set use '-i'"
    exit 1
fi
if [ $pam = "" ]; then
    echo "pam sequence is needed, set use '-p'"
    exit 1
fi
if [ $tlen = "" ]; then
    echo "target length (include PAM) is needed, set use '-L'"
    exit 1
fi
if [ $out_prefix = "" ]; then
    echo "out_prefix is needed, set use '-o'"
    exit 1
fi
if [ $index = "" ]; then
    echo "bowtie index is needed, set use '--bt'"
    exit 1
fi
if [ ! $pam_pos ]; then
    pam_pos=3
    echo "PAM positon is set as '3', you can reset it use '--pam-pos'"
fi
if [ ! $ref]; then
    ref=$fasta
fi
echo "in_fasta_file = $fasta"
echo "ref           = $ref"
echo "pam           = $pam"
echo "target_length = $tlen"
echo "out_prefix    = $out_prefix"
echo "bowtie_index  = $index"
echo "pam_pos       = $pam_pos"

bin_home=/media/f/my_program3/src/crispr
if [ $dry -eq 1 ]; then
    echo "
    $bin_home/search_pam --in-fasta $fasta --pam $pam --target-length $tlen --out-prefix $out_prefix. --pam-pos $pam_pos
    $bin_home/sort_and_count.sh $out_prefix.bed $out_prefix $index
    $bin_home/filt_pam --in-fasta $ref --pam $pam --pam-pos $pam_pos --in-bam $out_prefix.bam --out-prefix $out_prefix.
    $bin_home/sort_bam_by_name.sh $out_prefix.filtered.bam
    $bin_home/target_stat --in-bam $out_prefix.filtered.sn.bam --out-prefix $out_prefix.
    ";
else
    step=1;
    echo "run"

    echo "STEP $step:"
    step=`expr $step + 1`;
    echo "$bin_home/search_pam --in-fasta $fasta --pam $pam --target-length $tlen --out-prefix $out_prefix. --pam-pos $pam_pos"
    $bin_home/search_pam --in-fasta $fasta --pam $pam --target-length $tlen --out-prefix $out_prefix. --pam-pos $pam_pos

    echo "STEP $step:"
    step=`expr $step + 1`;
    echo "$bin_home/sort_and_count.sh $out_prefix.bed $out_prefix $index"
    $bin_home/sort_and_count.sh $out_prefix.bed $out_prefix $index

    echo "STEP $step:"
    step=`expr $step + 1`;
    echo "$bin_home/filt_pam --in-fasta $ref --pam $pam --pam-pos $pam_pos --in-bam $out_prefix.bam --out-prefix $out_prefix."
    $bin_home/filt_pam --in-fasta $ref --pam $pam --pam-pos $pam_pos --in-bam $out_prefix.bam --out-prefix $out_prefix.

    echo "STEP $step:"
    step=`expr $step + 1`;
    echo "$bin_home/sort_bam_by_name.sh $out_prefix.filtered.bam"
    $bin_home/sort_bam_by_name.sh $out_prefix.filtered.bam

    echo "STEP $step:"
    step=`expr $step + 1`;
    echo "$bin_home/target_stat --in-bam $out_prefix.filtered.sn.bam --out-prefix $out_prefix."
    $bin_home/target_stat --in-bam $out_prefix.filtered.sn.bam --out-prefix $out_prefix.
fi


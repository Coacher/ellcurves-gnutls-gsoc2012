#!/bin/bash

help() {
cat <<EOF
--------
$0: do benchmarking and generate reports

    Usage:
    benchmark [ -h ] [ -p PATCHDIR ] [ -s SOURCEDIR ] [ -o OUTDIR ]

    -h -- Prints this help message.

    -p -- Path to the dir with patches.

    -s -- Path to the dir with sources.

    -o -- Path to the dir for reports.

EOF
exit 0
}

PATCHDIR="$PWD"
SOURCEDIR="$PWD"
OUTDIR="$PWD"

# number of successive tests for each case
TESTRUNS="3"

# array of tested bit lengths
# for each number in this array there must be
# a file named "number_bits.patch"
BITS=(
"224"
"256"
"384"
"521" )

# array of tested features
# for each element in this array there must be
# a file named "element.patch"
PATCHES=(
"default"
"wmnaf"
)

die() {
    echo "die: $*" 1>&2
    exit 1
}

prepare_source() {
    WANT_AUTOMAKE="1.11" make autoreconf &> /dev/null || die "autoreconf failed"
}

configure_source() {
    echo -e -n "\tConfiguring sources...\t"
    WANT_AUTOMAKE="1.11" ./configure --disable-gtk-doc --disable-gcc-warnings &> /dev/null || die "configure failed"
    echo "done."
}

make_source() {
    echo -e -n "\tBuilding sources...\t"
    make &> /dev/null || die "make failed"
    echo "done."
}

clean_source() {
    make distclean &> /dev/null || die "make failed"
}

set_bit_length() {
    echo -e -n "Setting bit length to $1...\t"
    patch -s -Np1 -i "${PATCHDIR}/$1_bits.patch" || die "Failed to set bit length to $1"
    echo "done."
}

revert_bit_length() {
    patch -s -R -Np1 -i "${PATCHDIR}/$1_bits.patch" || die "Failed to revert bit length from $1"
}

enable_feature() {
    echo -e -n "Testing feature $1...\n"
    patch -s -Np1 -i "${PATCHDIR}/$1.patch" || die "Failed to apply patch for $1"
}

disable_feature() {
    patch -s -R -Np1 -i "${PATCHDIR}/$1.patch" || die "Failed to apply patch for $1"
}

make_header() {
    echo "# $1. ${TESTRUNS} successive run(s). ECC $2 bits." >> "${OUTDIR}/$1_results" || die "Failed to make header"
}

test_feature() {
    echo -e -n "\tRunning benchmark...\t"
    local i
    for (( i = 0; i < ${TESTRUNS}; ++i )); do
        ${SOURCEDIR}/src/gnutls-cli --benchmark-tls &>> "${OUTDIR}/$1_results" || die "Failed to run test"
    done
    echo "done."
}

flushmem() {
    if [ "${UID}" -eq 0 ]; then
        echo 3 > /proc/sys/vm/drop_caches
    else
        return 0
    fi
}

# main function
do_testing() {
    echo -e "\nBenchmark started on `date "+%d.%m.%Y"` at `date "+%T"`\n"

    cd ${SOURCEDIR} || die "Failed to step into ${SOURCEDIR}"
    prepare_source

    local n m
    for (( n = 0; n < ${#PATCHES[@]}; ++n )); do
        enable_feature ${PATCHES[n]}
        configure_source
        make_source

        # 192 is a default bit length
        make_header ${PATCHES[n]} 192

        flushmem
        test_feature ${PATCHES[n]}

        for (( m = 0; m < ${#BITS[@]}; ++m )); do
            set_bit_length ${BITS[m]}
            make_source

            make_header ${PATCHES[n]} ${BITS[m]}

            flushmem
            test_feature ${PATCHES[n]}

            revert_bit_length ${BITS[m]}
        done

        echo -e "\n"
    done
    echo -e "\nBenchmark finished on `date "+%d.%m.%Y"` at `date "+%T"`\n"
}

undo_testing() {
    echo -e "Cleaning up.\n"
    local n
    for (( n = ${#PATCHES[@]} - 1; n >= 0; --n )); do
        disable_feature ${PATCHES[n]}
    done
    clean_source
}

if [ $# -ne 0 ]; then
    while getopts ":hp:s:o:" Option
        do
            case $Option in
                h     ) help;;
                p     ) PATCHDIR="$OPTARG";;
                s     ) SOURCEDIR="$OPTARG";;
                o     ) OUTDIR="$OPTARG";;
                *     ) help;;
            esac
        done
    do_testing
    undo_testing
else
    do_testing
    undo_testing
fi

exit 0

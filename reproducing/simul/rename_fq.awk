#!/bin/awk -f

FNR == 1 {
    printf "[%s] Analyzing file %s\n", strftime("%H:%M:%S"), FILENAME > "/dev/tty"
    n = split(FILENAME, filename_split, "/")
    basename = filename_split[n]
    m = split(basename, basename_split, ".")
    replacement = ":"basename_split[1]"."basename_split[2]"\\1"
}
{
    if (NR % 4 == 1) {
        $0 = gensub(/((\/[12])?)$/, replacement, 1)
    }
    print
}

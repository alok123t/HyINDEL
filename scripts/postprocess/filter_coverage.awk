{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= 30) {print $0, v[NR] }}

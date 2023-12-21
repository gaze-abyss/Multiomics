export PATH=/path/SparCC3:$PATH
cd /path/

# step1
SparCC.py ./t1.txt -i 100 -c cor_sp_t1.txt > log
SparCC.py ./t2.txt -i 100 -c cor_sp_t2.txt >> log

# step2
MakeBootstraps.py t1.txt -n 1000 -t ./boot_#.txt -p boot1/ >> log
MakeBootstraps.py t2.txt -n 1000 -t ./boot_#.txt -p boot2/ >> log

# step3
for n in {0..999}; do echo $n && SparCC.py boot1/boot_${n}.txt -i 100 --cor_file=boot1/bootstrap_cor_${n}.txt >> log; done
for n in {0..999}; do echo $n && SparCC.py boot2/boot_${n}.txt -i 100 --cor_file=boot2/bootstrap_cor_${n}.txt >> log; done

# step4
PseudoPvals.py cor_sp_t1.txt boot1/bootstrap_cor_#.txt 1000 -o boot1/pvals.two_sided.txt -t two_sided >> log
PseudoPvals.py cor_sp_t2.txt boot2/bootstrap_cor_#.txt 1000 -o boot2/pvals.two_sided.txt -t two_sided >> log

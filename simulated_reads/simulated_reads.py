import os

ref_dir = "./references"
out_dir = "./simulated_reads"
os.makedirs(out_dir, exist_ok=True)

read_length = 150
read_pairs = 500
min_insert = 250
max_insert = 500
mean_quality = 33
seed = 42

for fname in os.listdir(ref_dir):
    if not fname.endswith(".fasta") and not fname.endswith(".fa"):
        continue

    base = os.path.splitext(fname)[0]
    ref_path = os.path.join(ref_dir, fname)
    out1 = os.path.join(out_dir, f"{base}_R1.fastq")
    out2 = os.path.join(out_dir, f"{base}_R2.fastq")

    cmd = (
        f"randomreads.sh ref={ref_path} out1={out1} out2={out2} "
        f"length={read_length} reads={read_pairs} paired=t "
        f"mininsert={min_insert} maxinsert={max_insert} "
        f"q={mean_quality} seed={seed} "
        f"adderrors"
    )

    print(cmd)

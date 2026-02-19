import argparse
import subprocess
import os


def run_macs3(bam_file, output_dir, genome, name_prefix, q):
    os.makedirs(output_dir, exist_ok=True)
    print("Running MACS3 peak calling")
    results = subprocess.run(
        ["macs3", "callpeak",
         "-t", bam_file,
         "-f", "BAM",
         "-g", genome,
         "-n", name_prefix,
         "--outdir", output_dir,
         "-q", str(q)],
        capture_output=True,
        text=True
    )
    print(f"MACS3 peak calling complete! Output in {output_dir}")
    if results.returncode != 0:
        print("Error:", results.stderr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MACS3 peak calling')
    parser.add_argument("-t", "--treatment", required=True, help="BAM file")
    parser.add_argument("-g", "--genome", default="hs", help="Genome Size")
    parser.add_argument("-n", "--name_prefix", default="p53", help="Name prefix")
    parser.add_argument("-q", "--qvalue", default="0.05", help="Q-Value (minimum FDR)")
    parser.add_argument("-o", "--outdir", default="macs3_output", help="Output directory")
    args = parser.parse_args()

    run_macs3(bam_file=args.treatment,
              output_dir=args.outdir,
              genome=args.genome,
              name_prefix=args.name_prefix,
              q=args.qvalue)
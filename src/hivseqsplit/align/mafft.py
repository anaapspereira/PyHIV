import subprocess
from pathlib import Path
from uuid import uuid4

from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord


# Align the test sequence with the reference sequence using MAFFT
def mafft_align(test_seq, ref_seq):
    temp_input_file = Path(f"{uuid4().hex}_input.fasta")
    temp_output_file = Path(f"{uuid4().hex}_output.fasta")
    try:
        with open(temp_input_file, 'w') as temp_file:
            SeqIO.write([SeqRecord(ref_seq, id="ref_seq"), SeqRecord(test_seq, id="test_seq")], temp_file, "fasta")

        mafft_cmd = ["mafft", "--localpair", "--maxiterate", "1000", temp_input_file]
        with open(temp_output_file, 'w') as output_file:
            subprocess.run(mafft_cmd, stdout=output_file, stderr=subprocess.PIPE, check=True)

        alignment = AlignIO.read(temp_output_file, "fasta")
        ref_aligned = None
        test_aligned = None
        for record in alignment:
            if record.id == "ref_seq":
                ref_aligned = str(record.seq)
            elif record.id == "test_seq":
                test_aligned = str(record.seq)

        if ref_aligned is None or test_aligned is None:
            print("Error: Reference sequence or test sequence not found in the alignment.")
            return None, None

        return test_aligned, ref_aligned
    except subprocess.CalledProcessError as e:
        print(f"Error in MAFFT alignment: {e}")
        return None, None
    finally:
        if temp_input_file.exists():
            temp_input_file.unlink()
        if temp_output_file.exists():
            temp_output_file.unlink()

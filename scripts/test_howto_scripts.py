"""Fixture tests for the docs/HOWTO.md helper scripts. Run: python3 -I scripts/test_howto_scripts.py (or python3 -m pytest scripts/)."""
import gzip, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import doc_blocks, make_families, one_protein_per_gene  # noqa: E402

# NCBI-style: gene g1 with two mRNAs (XP_2 is longer), gene g2 with one, a CDS hanging off its gene (g3, no mRNA row), a protein with an
# internal stop (g4), a CDS whose protein is missing from the FASTA (g5). Ensembl-style ids for gene e1 (two transcripts, equal length:
# the smaller id wins). XP_9 is in the FASTA only.
GFF = """##gff-version 3
c1\tRefSeq\tgene\t1\t900\t.\t+\t.\tID=gene-g1
c1\tRefSeq\tmRNA\t1\t900\t.\t+\t.\tID=rna-XM_1;Parent=gene-g1
c1\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=cds-XP_1;Parent=rna-XM_1;protein_id=XP_1
c1\tRefSeq\tmRNA\t1\t900\t.\t+\t.\tID=rna-XM_2;Parent=gene-g1
c1\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=cds-XP_2;Parent=rna-XM_2;protein_id=XP_2
c1\tRefSeq\tCDS\t60\t90\t.\t+\t0\tID=cds-XP_2;Parent=rna-XM_2;protein_id=XP_2
c1\tRefSeq\tgene\t1000\t1900\t.\t+\t.\tID=gene-g2
c1\tRefSeq\tmRNA\t1000\t1900\t.\t+\t.\tID=rna-XM_3;Parent=gene-g2
c1\tRefSeq\tCDS\t1000\t1030\t.\t+\t0\tID=cds-XP_3;Parent=rna-XM_3;protein_id=XP_3
c1\tRefSeq\tgene\t2000\t2900\t.\t+\t.\tID=gene-g3
c1\tRefSeq\tCDS\t2000\t2030\t.\t+\t0\tID=cds-XP_4;Parent=gene-g3;protein_id=XP_4
c1\tRefSeq\tmRNA\t3000\t3900\t.\t+\t.\tID=rna-XM_5;Parent=gene-g4
c1\tRefSeq\tCDS\t3000\t3030\t.\t+\t0\tID=cds-XP_5;Parent=rna-XM_5;protein_id=XP_5
c1\tRefSeq\tmRNA\t4000\t4900\t.\t+\t.\tID=rna-XM_6;Parent=gene-g5
c1\tRefSeq\tCDS\t4000\t4030\t.\t+\t0\tID=cds-XP_6;Parent=rna-XM_6;protein_id=XP_6
c2\tens\tmRNA\t1\t900\t.\t-\t.\tID=transcript:ET2;Parent=gene:e1
c2\tens\tCDS\t1\t30\t.\t-\t0\tID=CDS:EP2;Parent=transcript:ET2;protein_id=EP2
c2\tens\tmRNA\t1\t900\t.\t-\t.\tID=transcript:ET1;Parent=gene:e1
c2\tens\tCDS\t1\t30\t.\t-\t0\tID=CDS:EP1;Parent=transcript:ET1;protein_id=EP1
"""
FAA = ">XP_1 short isoform\nMAAA\n>XP_2 long isoform\nMAAAAA\nAA*\n>XP_3\nMCC*\n>XP_4\nMDD\n>XP_5\nMK*KK\n>XP_9 not in gff\nMEE\n>EP2\nMGG\n>EP1\nMHH\n"


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"): name = line[1:]; out[name] = ""
        else: out[name] += line
    return out


def test_one_protein_per_gene():
    with tempfile.TemporaryDirectory() as d:
        gff = os.path.join(d, "sp.gff.gz")
        with gzip.open(gff, "wt") as f: f.write(GFF)
        faa = os.path.join(d, "sp.protein.faa")
        with open(faa, "w") as f: f.write(FAA)
        out = os.path.join(d, "proteomes", "Sp.faa")
        one_protein_per_gene.main(["--gff", gff, "--proteins", faa, "--out", out])
        got = fasta(out)
        assert got == {"XP_2": "MAAAAAAA", "XP_3": "MCC", "XP_4": "MDD", "EP1": "MHH"}, got   # longest; '*' stripped; tie -> smaller id


def write_of(d, rows, hogs):
    os.makedirs(os.path.join(d, "Orthogroups")); os.makedirs(os.path.join(d, "Phylogenetic_Hierarchical_Orthogroups"))
    with open(os.path.join(d, "Orthogroups", "Orthogroups.tsv"), "w") as f:
        f.write("Orthogroup\tA\tB\tC\n" + "".join("\t".join(r) + "\n" for r in rows))
    with open(os.path.join(d, "Phylogenetic_Hierarchical_Orthogroups", "N1.tsv"), "w") as f:
        f.write("HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\n" + "".join("\t".join(r) + "\n" for r in hogs))


def test_make_families():
    with tempfile.TemporaryDirectory() as d:
        of = os.path.join(d, "of"); prot = os.path.join(d, "proteomes"); os.makedirs(prot)
        write_of(of, [("OG1", "a1", "b1", "c1"),          # strict
                      ("OG2", "a2", "b2", ""),            # partial (C absent)
                      ("OG3", "a3", "", ""),              # too small
                      ("OG4", "a4, a5", "b4", "c4")],     # multi-copy -> HOG
                 [("N1.HOG0", "OG4", "n0", "a4", "b4", "c4"), ("N1.HOG1", "OG4", "n0", "a5", "", "")])
        seqs = {"A": "a1 a2 a3 a4 a5", "B": "b1 b2 b4", "C": "c1 c4"}
        for lab, ids in seqs.items():
            ext = ".faa.gz" if lab == "C" else ".faa"              # both plain and gzip proteomes are found
            op = gzip.open if ext.endswith(".gz") else open
            with op(os.path.join(prot, lab + ext), "wt") as f:
                for i in ids.split():
                    f.write(f">{i}\n{'MK*K' if i == 'c4' else 'MKK*'}\n")
        out = os.path.join(d, "families")
        make_families.main(["--of", of, "--proteomes", prot, "--min", "2", "--out", out])
        ls = lambda t: sorted(os.listdir(os.path.join(out, t, "families")))
        assert ls("A_strict") == ["OG1.faa"]
        assert ls("B_partial") == ["OG2.faa"]
        assert ls("C_hog") == ["N1.HOG0.faa"]
        assert fasta(os.path.join(out, "C_hog", "families", "N1.HOG0.faa")) == {"a4": "MKK", "b4": "MKK"}   # c4: internal stop, dropped
        assert open(os.path.join(out, "A_strict", "taxon_map.tsv")).read() == "a1\tA\nb1\tB\nc1\tC\n"


def test_doc_blocks():
    doc = "# T\n\n```sh\nsetup\n```\n\n## 1. Go\n\n```sh\nX=1\necho $X\n```\n\n```text\nnot run\n```\n\n```sh\necho two\n```\n"
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as f: f.write(doc)
    try:
        b = doc_blocks.blocks(doc, "## 1.")
        assert b == [["X=1", "echo $X"], ["echo two"]]
        r = subprocess.run([sys.executable, "-I", os.path.join(HERE, "doc_blocks.py"), f.name, "--from", "## 1.", "--sub", "X=1", "X=2"],
                           capture_output=True, text=True, check=True)
        assert subprocess.run(["bash", "-c", r.stdout], capture_output=True, text=True).stdout == "2\ntwo\n"
        r = subprocess.run([sys.executable, "-I", os.path.join(HERE, "doc_blocks.py"), f.name, "--sub", "Y=1", "Y=2"], capture_output=True, text=True)
        assert r.returncode != 0 and "not found" in r.stderr
    finally:
        os.unlink(f.name)


if __name__ == "__main__":
    test_one_protein_per_gene(); test_make_families(); test_doc_blocks(); print("ok: 3 tests")

"""The Terminal tab's logic (runner.py): pipelines from tools.toml turned into bash, paths between
Windows and WSL, where `open` sends a file, and the QC summary read from the tools' output."""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import runner  # noqa: E402


def bash_syntax_ok(text):
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("no bash to check the syntax with")
    done = subprocess.run([bash, "-n", "-c", text], capture_output=True, text=True)
    return done.returncode == 0, done.stderr


def test_pipelines_build_valid_bash(tmp_path):
    reads = tmp_path / "my reads_R1.fq"
    files = {"reads1": str(reads), "reads2": "", "long": "", "ref": str(tmp_path / "ref.fa"),
             "seqs": str(tmp_path / "s.fa"), "out": str(tmp_path / "out.bam")}
    align = runner.pipelines()["align"]
    command, used = runner.build(align, runner.options_for(align, "Illumina paired-end"), files, 6, tmp_path / "run")
    assert used == ["fastp", "bwa", "samtools"]                        # only the tools this preset calls are checked
    assert command.startswith("( need fastp bwa samtools && cd ")
    assert "-w 6 -q 20 -l 50" in command and "bwa mem -t 6" in command  # nested placeholders filled
    assert "'" + runner.shell_path(reads) + "'" in command              # a path with a space is quoted
    assert '${R2:+-I "$R2"}' in command                                 # bash's own braces left alone
    nano = runner.options_for(align, "Nanopore")
    assert nano["cleanup"] == "filtlong" and nano["min_len"] == 50      # a preset's values, defaults for the rest
    assert runner.inputs_for(align, "Nanopore") == ["reads1", "reference"]
    # the diagram comes from tools.toml alone: choices named, tools with their subcommand, cleanup off shown as skipped
    assert runner.diagram(align, dict(nano, cleanup="none"), ["reads1", "reference"]) == [
        ("input", "Reads"), ("input", "Reference"), ("skipped", "Read cleanup: none"), ("step", "Aligner: minimap2_ont"),
        ("step", "samtools sort"), ("step", "samtools index"), ("step", "samtools flagstat"), ("output", "Open in openDNA")]
    asm = runner.pipelines()["assemble"]
    assert [t for k, t in runner.diagram(asm, runner.options_for(asm, "Hybrid (Illumina + Nanopore)"), ["long"]) if k == "step"] == [
        "Read cleanup: auto", "unicycler (Unicycler mode: normal)", "quast"]   # quast once: '-o quast' is not a second step
    # every preset of every pipeline: no placeholder left, and bash accepts it
    for key, pipeline in runner.pipelines().items():
        for preset in pipeline["presets"]:
            command, _ = runner.build(pipeline, runner.options_for(pipeline, preset), files, 2, tmp_path)
            assert not re.search(r"(?<!\$)\{\w+\}", command), (key, preset, command)
            ok, err = bash_syntax_ok(command)
            assert ok, (key, preset, err)
    ok, err = bash_syntax_ok(runner.session_setup(tmp_path) + runner.wrap("echo 'it''s' | wc -c") + runner.stop_children(42, 43)
                             + "\n" + runner.install(["mafft"]))
    assert ok, err
    assert "install" in runner.install(["mafft", "clustalo"]) and "create -y -n opendna-tools" in runner.install(["mafft"])
    assert runner.default_output(align, "pUC 19", str(reads), None) == str(tmp_path / "pUC_19_my reads.bam")
    assert runner.run_folder(tmp_path / "x.bam").name == "x_run"


def test_open_paths_and_summary(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "use_wsl", lambda: True)
    assert runner.shell_path("D:\\data\\x.fq").endswith("/mnt/d/data/x.fq")
    assert runner.windows_path("/mnt/d/data/x.fq") == "D:\\data\\x.fq"
    assert runner.windows_path("\\\\wsl$\\Ubuntu\\home\\a.bam") == "\\\\wsl$\\Ubuntu\\home\\a.bam"   # wslpath -w output kept
    monkeypatch.undo()
    aligned, plain = tmp_path / "a.fa", tmp_path / "p.fa"
    aligned.write_text(">a\nAC-GT\n>b\nACCGT\n")
    plain.write_text(">a\nACGT\n>b\nACCGT\n")
    assert [runner.open_kind(p) for p in (tmp_path / "x.bam", tmp_path / "v.vcf.gz", aligned, plain, tmp_path / "r.gb",
                                           tmp_path / "fastp.html")] == ["track", "track", "msa", "document", "document", "browser"]
    assert runner.open_kind(plain, msa=True) == "msa"
    assert not runner.is_alignment(">a\nAC-GT\n")                      # one record is not an alignment
    (tmp_path / "fastp.json").write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 1000},
                                                                 "after_filtering": {"total_reads": 920, "q30_rate": 0.941}}}))
    (tmp_path / "flagstat.txt").write_text("1000 + 0 in total (QC-passed reads + QC-failed reads)\n"
                                           "987 + 0 mapped (98.70% : N/A)\n")
    (tmp_path / "quast").mkdir()
    (tmp_path / "quast" / "report.tsv").write_text("Assembly\tassembly\n# contigs\t3\nTotal length\t4641652\nN50\t4600000\n")
    (tmp_path / "quast" / "report.html").write_text("<html></html>")
    (tmp_path / "quast" / "icarus_viewers").mkdir()
    (tmp_path / "quast" / "icarus_viewers" / "contig_size_viewer.html").write_text("")   # not a report
    assert runner.summary(tmp_path) == ["fastp: 1,000 reads in, 920 kept (92%), Q30 94%",
                                        "samtools flagstat: 987 reads mapped (98.70%)",
                                        "QUAST: 3 contigs, 4,641,652 bp, N50 4,600,000 bp"]
    assert runner.reports(tmp_path) == [str(tmp_path / "quast" / "report.html")]
    assert runner.summary(tmp_path / "quast") == []                   # nothing known there: no lines

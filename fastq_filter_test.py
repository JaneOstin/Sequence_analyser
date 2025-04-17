import pytest
import logging
import os

from seq_analyser import filter_fastq
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord


def create_test_fastq(file_path, records):
    with open(file_path, "w") as file:
        SeqIO.write(records, file, "fastq")


class TestFilterFastq:
    @pytest.fixture
    def tmp_fastq_file(self, tmp_path):
        records = [
            SeqRecord(Seq.Seq("ACGT"), id="seq1", description="record1", letter_annotations={"phred_quality": [40, 40, 40, 40]}),
            SeqRecord(Seq.Seq("GGCCGG"), id="seq2", description="record2", letter_annotations={"phred_quality": [10, 10, 10, 10, 10, 10]})
        ]
        file_path = tmp_path / "input.fastq"
        create_test_fastq(str(file_path), records)
        return str(file_path)
    
    @pytest.fixture
    def empty_fastq_file(self, tmp_path):
        file_path = tmp_path / "empty.fastq"
        with open(str(file_path), "w") as f:
            f.write("")
        return str(file_path)

    def test_filtering_correct(self, tmp_fastq_file):
        output_filename = "filtered_output.fastq"
        filter_fastq(tmp_fastq_file, output_fastq=output_filename, quality_threshold=30)
        filtered_dir = os.path.join(os.path.dirname(tmp_fastq_file), "filtered")
        output_path = os.path.join(filtered_dir, output_filename)
        assert os.path.exists(output_path), "No output fastq file was created"
        records = list(SeqIO.parse(output_path, "fastq"))
        assert len(records) == 1, "Filtration doesn't work as expected"
        assert records[0].id == "seq1", "Filtration doesn't work as expected"

    def test_file_logging(self, tmp_fastq_file):
        logging.getLogger().handlers.clear()
        base_dir = '/'.join(tmp_fastq_file.split("/")[:-1]) + '/filtered'
        os.makedirs(base_dir, exist_ok=True)
        filter_fastq(tmp_fastq_file, output_fastq='', quality_threshold=30)
        log_file = os.path.join(base_dir, "fastq-filter.log")
        assert os.path.exists(log_file), "No log file was created"

    def test_empty_input_file(self, empty_fastq_file):
        output_filename = "empty_output.fastq"
        filter_fastq(empty_fastq_file, output_fastq=output_filename)
        filtered_dir = os.path.join(os.path.dirname(empty_fastq_file), "filtered")
        output_path = os.path.join(filtered_dir, output_filename)
        assert os.path.exists(output_path), "No output file was created for empty input file"
        records = list(SeqIO.parse(output_path, "fastq"))
        assert len(records) == 0, "Output file should be empty"

    def test_gc_boundaries(self, tmp_fastq_file):
        output_filename = "filtered_output.fastq"
        filter_fastq(tmp_fastq_file, output_fastq=output_filename, gc_bounds=(50, 100))
        filtered_dir = os.path.join(os.path.dirname(tmp_fastq_file), "filtered")
        output_path = os.path.join(filtered_dir, output_filename)
        records = list(SeqIO.parse(output_path, "fastq"))
        assert len(records) == 2, "Both records should pass filtration"

    def test_length_boundaries(self, tmp_fastq_file):
        output_filename = "filtered_output.fastq"
        filter_fastq(tmp_fastq_file, output_fastq=output_filename, length_bounds=(4, 6))
        filtered_dir = os.path.join(os.path.dirname(tmp_fastq_file), "filtered")
        output_path = os.path.join(filtered_dir, output_filename)
        records = list(SeqIO.parse(output_path, "fastq"))
        assert len(records) == 2, "Records on boundaties should be included"

    def test_non_tuple_length_bounds(self, tmp_fastq_file):
        output_filename = "filtered_output.fastq"
        filter_fastq(tmp_fastq_file, output_fastq=output_filename, length_bounds=4)
        filtered_dir = os.path.join(os.path.dirname(tmp_fastq_file), "filtered")
        output_path = os.path.join(filtered_dir, output_filename)
        records = list(SeqIO.parse(output_path, "fastq"))
        assert len(records) == 1, "Record must pass filtration"

    def test_incorrect_quality_threshold_type(self, tmp_fastq_file):
        with pytest.raises(TypeError):
            filter_fastq(tmp_fastq_file, quality_threshold="high")

    def test_non_existent_input_file(self, tmp_path):
        non_existent_file = str(tmp_path / "non_existent.fastq")
        with pytest.raises(FileNotFoundError):
            filter_fastq(non_existent_file)

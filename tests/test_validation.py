import pytest

from cytetype.preprocessing.validation import _is_gene_id_like, _id_like_percentage


class TestIsGeneIdLike:
    @pytest.mark.parametrize(
        "value",
        [
            "ENSG00000000003",
            "ENSG00000000003.14",
            "ENSMUSG00000000001",
            "ensg00000000003",
        ],
    )
    def test_ensembl_ids(self, value: str) -> None:
        assert _is_gene_id_like(value) is True

    @pytest.mark.parametrize(
        "value",
        [
            "NM_001301",
            "NR_046018",
            "XM_011541",
            "XR_001737",
        ],
    )
    def test_refseq_ids(self, value: str) -> None:
        assert _is_gene_id_like(value) is True

    @pytest.mark.parametrize(
        "value",
        [
            "7157",
            "672",
            "11286",
            "0",
        ],
    )
    def test_integer_entrez_ids(self, value: str) -> None:
        assert _is_gene_id_like(value) is True

    @pytest.mark.parametrize(
        "value",
        [
            "7157.0",
            "672.0",
            "11286.0",
            "0.0",
        ],
    )
    def test_float_stringified_entrez_ids(self, value: str) -> None:
        assert _is_gene_id_like(value) is True

    @pytest.mark.parametrize(
        "value",
        [
            "AFFY_HG_U133A.207163_S_AT",
            "ILLUMINA_HUMANHT_12_V4.ILMN_1762337",
        ],
    )
    def test_long_dotted_ids(self, value: str) -> None:
        assert _is_gene_id_like(value) is True

    @pytest.mark.parametrize(
        "value",
        [
            "TSPAN6",
            "DPM1",
            "SCYL3",
            "TP53",
            "BRCA1",
            "CD8A",
            "MS4A1",
        ],
    )
    def test_gene_symbols_not_flagged(self, value: str) -> None:
        assert _is_gene_id_like(value) is False

    @pytest.mark.parametrize(
        "value",
        [
            "",
            "   ",
            "7157.5",
        ],
    )
    def test_edge_cases(self, value: str) -> None:
        assert _is_gene_id_like(value) is False


class TestIdLikePercentage:
    def test_all_gene_symbols(self) -> None:
        values = ["TSPAN6", "DPM1", "SCYL3", "TP53", "BRCA1"]
        assert _id_like_percentage(values) == 0.0

    def test_all_ensembl_ids(self) -> None:
        values = [f"ENSG{i:011d}" for i in range(20)]
        assert _id_like_percentage(values) == 100.0

    def test_all_integer_entrez_ids(self) -> None:
        values = ["7157", "672", "3845", "11286", "9952"]
        assert _id_like_percentage(values) == 100.0

    def test_all_float_entrez_ids(self) -> None:
        values = ["7157.0", "672.0", "3845.0", "11286.0", "9952.0"]
        assert _id_like_percentage(values) == 100.0

    def test_mixed_float_entrez_and_symbols(self) -> None:
        symbols = ["TSPAN6", "DPM1", "SCYL3"]
        entrez = ["7157.0", "672.0", "3845.0", "11286.0", "9952.0", "904.0", "405.0"]
        pct = _id_like_percentage(symbols + entrez)
        assert pct == 70.0

    def test_empty_list(self) -> None:
        assert _id_like_percentage([]) == 100.0

import json, tempfile, unittest
from pathlib import Path
from sf_csa.core import (
    DEFAULT_MECHANISM_FAMILIES,
    classify_hit,
    classify_title,
    sequence_sha,
    structural_category,
    verify_release,
)

class CoreTests(unittest.TestCase):
    def test_whole_architecture_requires_both_coverages(self):
        h={"alntmscore":"0.8","qcov":"0.9","tcov":"0.6"}
        self.assertEqual(structural_category(h,0.7,0.5),"domain_or_partial_match")

    def test_title_does_not_promote_function(self):
        q={"accession":"Q","mechanism_group":"tonb_heme_receptor"}
        h={"target":"X","theader":"toluene transporter outer membrane barrel"}
        label,*_=classify_hit(q,h,"whole_architecture_match")
        self.assertNotIn(label,{"exact_function_supported","probable_same_function"})

    def test_msp_stays_unresolved(self):
        q={"accession":"M","mechanism_group":"msp_contested"}
        h={"target":"X","theader":"major surface protein Msp"}
        label,*_=classify_hit(q,h,"whole_architecture_match")
        self.assertEqual(label,"unresolved_or_conflicted")


class MechanismFamilyTests(unittest.TestCase):
    """The family table moved from literals in classify_title into the manifest."""

    def test_refine_narrows_a_matched_family(self):
        # A TonB-dependent transporter is an importer; a heme one is a heme receptor.
        self.assertEqual(classify_title("TonB-dependent receptor"), "susc_raga_importer")
        self.assertEqual(classify_title("TonB-dependent hemoglobin receptor"), "tonb_heme_receptor")

    def test_defaults_and_explicit_families_agree(self):
        for title in ("BamA outer membrane protein assembly factor", "PorG type IX secretion",
                      "major surface protein Msp", "trimeric porin", "toluene transporter"):
            self.assertEqual(
                classify_title(title), classify_title(title, DEFAULT_MECHANISM_FAMILIES), title
            )

    def test_a_campaign_can_replace_the_table(self):
        families = [{"group": "flagellin", "pattern": r"flagell"}]
        self.assertEqual(classify_title("flagellin FliC", families), "flagellin")
        # And the periodontal families no longer apply.
        self.assertEqual(classify_title("BamA omp85", families), "unknown")


class VerifyReleaseTests(unittest.TestCase):
    """Expectations come from the database manifest, never from the release itself."""

    def _release(self, root: Path, query_count: int) -> Path:
        output = root / "release"
        (output / "targets").mkdir(parents=True)
        (output / "SF_CSA_RELEASE_MANIFEST.json").write_text(json.dumps({
            "query_count": query_count, "proteome_count": 63,
            "target_statuses": {"WP_014224802.1": "CONDITIONAL"}, "queries": [],
        }))
        return output

    def _databases(self, root: Path, expected_queries: int) -> Path:
        path = root / "database_manifest.json"
        path.write_text(json.dumps({"release_expectations": {
            "query_count": expected_queries, "proteome_count": 63,
            "target_statuses": {"WP_014224802.1": "CONDITIONAL"},
        }}))
        return path

    def test_a_release_matching_its_config_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.assertEqual(verify_release(self._release(root, 8), self._databases(root, 8)), [])

    def test_a_release_of_the_wrong_size_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            errors = verify_release(self._release(root, 7), self._databases(root, 8))
            self.assertTrue(any("expected 8 targets" in e for e in errors), errors)

    def test_without_a_config_the_shape_is_reported_unverified(self):
        """A release must not be allowed to supply its own expected shape."""
        with tempfile.TemporaryDirectory() as tmp:
            errors = verify_release(self._release(Path(tmp), 8))
            self.assertTrue(any("was not verified" in e for e in errors), errors)


if __name__=="__main__": unittest.main()

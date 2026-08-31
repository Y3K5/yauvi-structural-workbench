window.YAUVI_PUBLIC_SHOWCASE = {
  "baseline": {
    "baseline_id": "structural-workbench-offline-qualification-v2-2026-08-26",
    "scientific_boundary": "Passing software tests are not passing external scientific qualification benchmarks.",
    "selection": "not network and not adapter",
    "total_deselected": 6,
    "total_passed": 495
  },
  "cases": [
    {
      "analysis_type": "structure_qc",
      "case_id": "HUC-01",
      "evidence_files": [
        {
          "label": "STRUCTURE_EVIDENCE.json",
          "path": "evidence/HUC-01/STRUCTURE_EVIDENCE.json",
          "sha256": "dffaf668899d38d5d8fb5afb0e3a7c91b91c04fd7aed1386aa6518f72bdae254"
        },
        {
          "label": "RESIDUE_QUALITY.tsv",
          "path": "evidence/HUC-01/RESIDUE_QUALITY.tsv",
          "sha256": "6cba8ab043109e9914b626ca64d77308c82494591533db3d4353cfca39ff142a"
        }
      ],
      "human_benefits": [
        "Prevents residue-numbering mistakes before mutation or site analysis.",
        "Makes experimental and predicted confidence interpretable without mixing them.",
        "Provides a reproducible accept, hold, or investigate-first boundary."
      ],
      "human_label": "Can I safely interpret these coordinates?",
      "human_question": "Before mapping a variant or functional residue, are sequence identity, residue numbering, provenance, and validation bound to the exact model?",
      "input_sha256": {
        "qc/model.pdb": "a598a5203e772bb96e4d2005f1f5cd5cb6d4d753511fce064025a258f34ac4b0",
        "qc/provenance.json": "f549541b8ef835652b46e8b803d5815344ea006e77dd8c371c1fca076531b9b7",
        "qc/reference.fasta": "daf3cb9543fe62dfb2104efcf1ea54cd710eedfdbf1b5c3b1001bc6afd06cc26",
        "qc/validation.json": "32ecc8b222d8815d2e46d6c10b55361ef1684f7804f6412071ebcd3904bff205"
      },
      "measurements": [
        {
          "help": "Reference residues with mapped coordinates",
          "label": "Sequence coverage",
          "value": "100.0%"
        },
        {
          "help": "Identity within the accepted mapping",
          "label": "Sequence identity",
          "value": "100.0%"
        },
        {
          "help": "Imported, not recomputed by StructQC",
          "label": "Validation",
          "value": "imported"
        },
        {
          "help": "Detected coordinate discontinuities",
          "label": "Chain breaks",
          "value": "0"
        }
      ],
      "non_claim": "Native conformation, biological function, or experimental correctness beyond the imported validation record.",
      "observed_result": "Both synthetic reference residues mapped exactly; predicted provenance and imported validation remained explicit.",
      "test_state": "passed_synthetic_case",
      "tool": "StructQC"
    },
    {
      "analysis_type": "membrane_orientation",
      "case_id": "HUC-02",
      "evidence_files": [
        {
          "label": "MEMBRANE_ORIENTATION.json",
          "path": "evidence/HUC-02/MEMBRANE_ORIENTATION.json",
          "sha256": "9f7e22278fad595725fea60e7d1018c5128fafd793d7ed3fbce143d7a83bbb95"
        },
        {
          "label": "RESIDUE_ORIENTATION.tsv",
          "path": "evidence/HUC-02/RESIDUE_ORIENTATION.tsv",
          "sha256": "cde3e060b467c5cabd52a649e24975827a7f51f2c198dfb0a3fc9719e4dba0c2"
        }
      ],
      "human_benefits": [
        "Helps choose candidate extracellular loops for follow-up assays.",
        "Separates membrane-core residues from flanking domains for construct planning.",
        "Provides a common frame for comparing mutations or models."
      ],
      "human_label": "Which side of a membrane might a receptor expose?",
      "human_question": "Can a single-pass membrane protein be placed in a consistent coordinate frame and divided into membrane and sided residue sets?",
      "input_sha256": {
        "membrane/synthetic_tm_receptor.pdb": "00c5d35eb773846d18b2b80727650a1f4fa38b9a3c59f12a29916622be20ca1a",
        "membrane/synthetic_tm_topology.json": "10fd5fbef527426c374d5b7ca29533f72f216f2b55bc967ab145ffcf31862829"
      },
      "measurements": [
        {
          "help": "Geometry route selected by the tool",
          "label": "Structure label",
          "value": "tm_helix_experimental"
        },
        {
          "help": "Named computational method",
          "label": "Orientation method",
          "value": "tm_helix_axis_v2"
        },
        {
          "help": "Coordinate-bound residue annotations",
          "label": "Residues reviewed",
          "value": "57"
        },
        {
          "help": "Geometry-derived candidates; not intact-cell evidence",
          "label": "Modeled surface set",
          "value": "28"
        }
      ],
      "non_claim": "Native topology, intact-cell exposure, antibody accessibility, expression, receptor function, or Mark 1 alpha-helical qualification.",
      "observed_result": "The invented receptor was routed as tm_helix_experimental using the experimental tm_helix_axis_v2 path with checksum-bound synthetic spans.",
      "test_state": "passed_synthetic_case",
      "tool": "MembraneOrient"
    },
    {
      "analysis_type": "conformational_state",
      "case_id": "HUC-03",
      "evidence_files": [
        {
          "label": "STATE_ENSEMBLE.json",
          "path": "evidence/HUC-03/STATE_ENSEMBLE.json",
          "sha256": "aecbaa94b27aa5f121d77fcb509a267c7c78f935849a2a05e10f351144d561b5"
        },
        {
          "label": "FRAME_METRICS.tsv",
          "path": "evidence/HUC-03/FRAME_METRICS.tsv",
          "sha256": "41d8d79d763522575ee9658795a32cbcf9c64c1339c064e350b33f3cfc1e124b"
        }
      ],
      "human_benefits": [
        "Compares mutant, apo, ligand-bound, or ensemble structures consistently.",
        "Identifies ambiguous frames instead of forcing a state call.",
        "Supports selection of conformations for further simulation or experiments."
      ],
      "human_label": "Which experimental conformation does my model resemble?",
      "human_question": "When two bounded reference states are declared, which reference is geometrically closer and is the margin interpretable?",
      "input_sha256": {
        "portfolio/inactive_reference.pdb": "38b2c93aab091181db5a29cb23a60c925f5c749ce3c184c66ef7c65ccf16a9c0",
        "portfolio/provenance.json": "a828e7f096f984f5d09ba8f24e53cc140d4e6d2ee91319d03a01644463d01349",
        "portfolio/query.pdb": "5bbc5614de1bacaee2aa4fa595d66053826ee86536e0bc6eb7c9a183ac51d82d",
        "portfolio/reference.fasta": "642b170bab9800ff5cf1fdaab5dfd86ce8ba782b855a20ca15cbb62f8c4c4902",
        "portfolio/reference_set.json": "95f98e8d40655a8bf3e74d7ad77e38a7417aa3125df87454fcd0bdd7fc3cd1f2"
      },
      "measurements": [
        {
          "help": "Bounded structural label",
          "label": "Resemblance label",
          "value": "active_like"
        },
        {
          "help": "After declared C-alpha alignment",
          "label": "Best RMSD",
          "value": "0.000 Å"
        },
        {
          "help": "Distance separation from the alternative reference",
          "label": "Reference margin",
          "value": "1.672 Å"
        },
        {
          "help": "Unresolved frames would remain in the denominator",
          "label": "Interpretable frames",
          "value": "1/1"
        }
      ],
      "non_claim": "Biochemical activity, activation, inhibition, mechanism, efficacy, or a time-resolved transition pathway.",
      "observed_result": "The synthetic query was active_like; all 1 frame(s) were retained in the population accounting.",
      "test_state": "passed_synthetic_case",
      "tool": "StateAtlas"
    },
    {
      "analysis_type": "functional_site_state",
      "case_id": "HUC-04",
      "evidence_files": [
        {
          "label": "SITE_CONTEXT.json",
          "path": "evidence/HUC-04/SITE_CONTEXT.json",
          "sha256": "49650292449ea356a2ed22087cb44c7a096a77fa23c0125ef3dddbcc1b3f1464"
        },
        {
          "label": "SITE_RESIDUES.tsv",
          "path": "evidence/HUC-04/SITE_RESIDUES.tsv",
          "sha256": "6799b4dd64a427a8143bbbd0e073c91dffd6f3d9be4f5453fcf67174bc25793d"
        }
      ],
      "human_benefits": [
        "Reveals missing or mismapped catalytic and binding residues.",
        "Supports control-mutation and construct-boundary discussions.",
        "Keeps annotations, observed chemistry, and pocket predictions separate."
      ],
      "human_label": "Are the declared functional residues structurally present?",
      "human_question": "Do curated residues map exactly, retain role-compatible identities, and overlap a separately identified pocket?",
      "input_sha256": {
        "portfolio/annotations.json": "396bdb4fdb9c299db0fd71b65ec5413ab73a1eae2fb84ca60c65188f76eea82d",
        "portfolio/pockets.json": "13a802fe52338c9ccc01699939fb79c958a3576c1ca3d383961bc93e3f8701d9",
        "portfolio/provenance.json": "a828e7f096f984f5d09ba8f24e53cc140d4e6d2ee91319d03a01644463d01349",
        "portfolio/query.pdb": "5bbc5614de1bacaee2aa4fa595d66053826ee86536e0bc6eb7c9a183ac51d82d",
        "portfolio/reference.fasta": "642b170bab9800ff5cf1fdaab5dfd86ce8ba782b855a20ca15cbb62f8c4c4902"
      },
      "measurements": [
        {
          "help": "Declared residues found in the structure",
          "label": "Mapped sites",
          "value": "3"
        },
        {
          "help": "Identity fits the declared role vocabulary",
          "label": "Role compatible",
          "value": "3"
        },
        {
          "help": "Descriptive C-alpha geometry",
          "label": "Maximum separation",
          "value": "5.295 Å"
        },
        {
          "help": "Scores stay specific to each named tool",
          "label": "Pocket methods",
          "value": "fpocket"
        }
      ],
      "non_claim": "Observed catalysis, ligand affinity, inhibition, druggability, physiological function, or clinical relevance.",
      "observed_result": "3 of 3 synthetic site residues were role-compatible; pocket evidence remained method-specific.",
      "test_state": "passed_synthetic_case",
      "tool": "SiteContext"
    },
    {
      "analysis_type": "assembly_interface",
      "case_id": "HUC-05",
      "evidence_files": [
        {
          "label": "ASSEMBLY_CONTEXT.json",
          "path": "evidence/HUC-05/ASSEMBLY_CONTEXT.json",
          "sha256": "2139adacffaf600958061b95fd853f73ed7934dd9dfba28cef8458ed14419cf8"
        },
        {
          "label": "INTERFACES.tsv",
          "path": "evidence/HUC-05/INTERFACES.tsv",
          "sha256": "48f2f9e392848ad95d929d6af72f0ebb755f85408d3037541134c647ccd2d795"
        }
      ],
      "human_benefits": [
        "Identifies candidate interface-disrupting or interface-preserving mutations.",
        "Shows when a site or surface is occluded by an assembly partner.",
        "Separates monomer interpretation from oligomeric structural context."
      ],
      "human_label": "Which residues become part of an oligomer interface?",
      "human_question": "In a declared two-chain assembly, which subject residues contact a partner and how much surface becomes buried?",
      "input_sha256": {
        "portfolio/assembly.pdb": "f2ad4d9e50daf000069ce25310fc6585401308822821b492954e41b10d8b3fc6",
        "portfolio/provenance.json": "a828e7f096f984f5d09ba8f24e53cc140d4e6d2ee91319d03a01644463d01349",
        "portfolio/query.pdb": "5bbc5614de1bacaee2aa4fa595d66053826ee86536e0bc6eb7c9a183ac51d82d",
        "portfolio/reference.fasta": "642b170bab9800ff5cf1fdaab5dfd86ce8ba782b855a20ca15cbb62f8c4c4902"
      },
      "measurements": [
        {
          "help": "Expected and observed chain inventory agree",
          "label": "Assembly complete",
          "value": "true"
        },
        {
          "help": "Subject residues with a partner inside the cutoff",
          "label": "Contact residues",
          "value": "3"
        },
        {
          "help": "Method-specific buried SASA",
          "label": "Buried surface",
          "value": "62.106 Å²"
        },
        {
          "help": "FreeSASA is preferred when installed",
          "label": "SASA method",
          "value": "freesasa_lee_richards_default_single_thread"
        }
      ],
      "non_claim": "Native oligomer abundance, binding affinity, intact-cell accessibility, physiological interaction, or mechanism.",
      "observed_result": "3 synthetic subject residues contacted chain B and 62.106 Å² became buried.",
      "test_state": "passed_synthetic_case",
      "tool": "AssemblyContext"
    },
    {
      "analysis_type": "sf_csa",
      "case_id": "HUC-06",
      "evidence_files": [
        {
          "label": "SF_CSA_RELEASE_MANIFEST.json",
          "path": "evidence/HUC-06/release/SF_CSA_RELEASE_MANIFEST.json",
          "sha256": "a691d1364cce0b70e9a119549ba781f875a0f1bf4404ee87b929192e50de75a0"
        },
        {
          "label": "RELEASE_COMPARISON_MATRIX.tsv",
          "path": "evidence/HUC-06/release/RELEASE_COMPARISON_MATRIX.tsv",
          "sha256": "80f326307e5b9f8b8984afcacd0cb60a535a496dec155febac5bc57e94559cdb"
        },
        {
          "label": "targets · QRY_A · structure_hits.tsv",
          "path": "evidence/HUC-06/release/targets/QRY_A/structure_hits.tsv",
          "sha256": "e6eafe75bba0763b250f0546b190c6076479f80d970807f40626c7552602788e"
        },
        {
          "label": "targets · QRY_A · species_comparison.tsv",
          "path": "evidence/HUC-06/release/targets/QRY_A/species_comparison.tsv",
          "sha256": "5beccaa6df68fa793919e23bcec14b93783b1866cf2b3597f3ad14961ee8da88"
        },
        {
          "label": "targets · QRY_B · structure_hits.tsv",
          "path": "evidence/HUC-06/release/targets/QRY_B/structure_hits.tsv",
          "sha256": "04ba36169e12b9baf4486a12e456acc99d98faa8117f6d31ff76bb48672f33bd"
        },
        {
          "label": "targets · QRY_B · species_comparison.tsv",
          "path": "evidence/HUC-06/release/targets/QRY_B/species_comparison.tsv",
          "sha256": "f134549abd96dc4969018881105c164281039eb6dd9f778ee3b27d905b7d048e"
        },
        {
          "label": "CHECKSUMS.json",
          "path": "evidence/HUC-06/release/CHECKSUMS.json",
          "sha256": "eb83999d30d1c2ddb9e0dae01a50b9e0e7df3c2fd76d52589f2a5f2f2f2e348f"
        }
      ],
      "human_benefits": [
        "Shows where fold similarity and sequence homology agree or disagree.",
        "Prevents a strong structural match from silently becoming exact functional transfer.",
        "Preserves missing structures and unresolved relationships as visible evidence gaps."
      ],
      "human_label": "How does structure-based evidence differ from sequence-based evidence?",
      "human_question": "Can the same invented proteins be compared through separate structural and sequence legs without turning similarity into exact functional transfer?",
      "input_sha256": {
        "inputs/database_manifest.json": "d9760fa72a50abf7cb837e44ab2aebbad721842286d0742104b9ac6db23c907a",
        "inputs/queries/QRY_A.faa": "f7e596dd8deec8b294dc3b5fa2a67032138505467c44d4bb53b33e5f421ccff3",
        "inputs/queries/QRY_A.pdb": "11dae7081d173d7f06e79ac38cf990ce46dd8c8fb7460dea0eae2bd2e4ec246d",
        "inputs/queries/QRY_B.faa": "d0bd2cb71218974e243bad81a28b8e69656d4ec1122f46ce70c2c23d877c9096",
        "inputs/queries/QRY_B.pdb": "ed2eba7ea57ef77e13eee2dc644c7c5840d878c6029686993de267447e52a8b9",
        "inputs/query_manifest.json": "907364e1957383cf5abbfc7328b5fc9a85c524fbf04f5d4147755a6d2405c56d",
        "runtime-fixture/hits.json": "447fc07fe07dfd6aad9fa78ee2a7dbba30eeeb09d7c05a1d93e8308d0672aa24"
      },
      "known_findings": [
        "The title-trap protection currently acts during release verification, not direct classify_hit calls.",
        "The main pipeline does not currently feed computed sequence reciprocal-best-hit evidence back into structural classification.",
        "External mini-database runs with the installed Foldseek and DIAMOND binaries remain pending."
      ],
      "measurements": [
        {
          "help": "Canonical sf-csa verify returned exit code 0",
          "label": "Release audit",
          "value": "passed"
        },
        {
          "help": "Invented checksum-bound query structures",
          "label": "Queries",
          "value": "2"
        },
        {
          "help": "Canned Foldseek-shaped rows interpreted by the real pipeline",
          "label": "Structural rows",
          "value": "6"
        },
        {
          "help": "Canned DIAMOND-shaped rows retained in a separate table",
          "label": "Sequence rows",
          "value": "3"
        },
        {
          "help": "Structural similarity and sequence homology are not merged",
          "label": "Evidence legs",
          "value": "2 separate"
        },
        {
          "help": "No Foldseek or DIAMOND alignment was computed in this case",
          "label": "Alignment engines",
          "value": "stubbed"
        }
      ],
      "non_claim": "Real alignment performance, biological function, orthology, substrate transfer, activity, pathogenic importance, or external scientific qualification.",
      "observed_result": "2 invented queries completed a checksum-bound SF-CSA release; the release audit passed and structural classifications remained separate from sequence orthology candidates.",
      "runtime_disclosure": "The canonical pipeline, subprocess construction, TSV parsing, classification, checksums, and release audit ran. Foldseek and DIAMOND were deterministic test doubles and computed no alignments.",
      "test_state": "passed_stubbed_pipeline_case",
      "tool": "SF-CSA"
    }
  ],
  "citation": {
    "author": "Yuvraj Patel",
    "license": "Apache-2.0",
    "orcid": "https://orcid.org/0009-0002-2276-7336",
    "title": "YAUVI Structural Biology Platform — Mark 1",
    "version": "0.1.0.dev0"
  },
  "data_class": "synthetic_demonstrations_plus_public_qualification_summaries",
  "limitations": [
    "Invented coordinates validate software behavior only.",
    "Benefits are potential research uses, not findings from these synthetic cases.",
    "Passing these cases does not replace external scientific qualification benchmarks.",
    "The SF-CSA public case verifies orchestration and evidence boundaries with test doubles; it does not benchmark real alignments.",
    "Synthetic demonstrations do not externally qualify a workflow; independent public cases are displayed separately.",
    "Four passing public cases do not establish workflow-general accuracy, and the two partial cases remain release blockers.",
    "Qualification v2 freezes the expanded panels and remains blocked until source-locked public cases are adopted and executed.",
    "Four of six Qualification v2 panels have executed and passed. Executed panels passing is not scope qualification: two panels are unadopted, membrane covers only its beta_barrel stratum, and no scope has reproduced on an independent second machine."
  ],
  "non_claims": [
    "Structural resemblance is not biochemical activity.",
    "Membrane orientation is not native exposure.",
    "A mapped site is not observed catalysis.",
    "An interface is not binding affinity.",
    "Structural or sequence similarity is not exact functional transfer."
  ],
  "platform_identity": {
    "display_name": "YAUVI Structural Biology Platform — Mark 1",
    "distribution": "yauvi-structural-workbench",
    "distribution_scope": "Command-line distribution. The loopback browser workbench is excluded: its controller imports private control-plane modules outside the published boundary.",
    "edition": "Mark 1",
    "identity_policy": "The Mark 1 name identifies the integrated platform. Standalone package, CLI, module, evidence-contract, and scientific-method names remain unchanged.",
    "platform_id": "yauvi_structural_biology_platform_mark_1",
    "primary_name": "YAUVI Structural Biology Platform",
    "publication_authorized": false,
    "release_state": "pre_public_preparation",
    "schema_version": "1.0",
    "scientific_suite_name": "YAUVI Structural Workbench",
    "share_non_claim": "It is not a clinical tool, a biochemical activity assay, or a universal protein-scoring system.",
    "share_status": "Mark 1 is a pre-public scientific build. Historical named cases remain visible; the expanded Qualification v2 panels are frozen but not yet source-adopted or executed.",
    "share_summary": "A local-first structural bioinformatics platform that turns protein coordinate files into inspectable, checksum-bound evidence across six analysis workflows.",
    "short_name": "YAUVI SBP Mark 1",
    "start_command": "python -m pip install -e \".[dev]\" && structqc describe",
    "tagline": "From protein coordinates to inspectable structural evidence."
  },
  "product": "YAUVI Structural Biology Platform — Mark 1",
  "publication_roadmap": {
    "current_phase": "local_hardening",
    "edition": "Mark 1",
    "gates": [
      {
        "evidence": "Four passed and two partial",
        "gate_id": "six_scientific_cases",
        "label": "Six external scientific case gates",
        "state": "partial"
      },
      {
        "evidence": "Current local reviewer selection: 524 passed and 6 network or adapter tests deselected; cross-platform clean-install reproduction remains a separate gate",
        "gate_id": "offline_software_baseline",
        "label": "Complete recorded offline regression",
        "state": "passed"
      },
      {
        "evidence": "Not started or not verified",
        "gate_id": "public_history",
        "label": "More than six months of active public history",
        "state": "blocked"
      },
      {
        "evidence": "Not recorded",
        "gate_id": "independent_use",
        "label": "Independent installation and research use",
        "state": "blocked"
      },
      {
        "evidence": "Incomplete",
        "gate_id": "license_audit",
        "label": "License and third-party redistribution audit",
        "state": "blocked"
      },
      {
        "evidence": "The 1,095-word paper compiles with official Inara and the four-page draft was visually reviewed; final approvals, journal metadata, and exact tool records remain",
        "gate_id": "paper_and_disclosures",
        "label": "Paper, authorship, funding, conflicts, and AI disclosure",
        "state": "partial"
      },
      {
        "evidence": "Not granted by local preparation",
        "gate_id": "publication_approval",
        "label": "Explicit approval for the exact public repository and release",
        "state": "blocked"
      }
    ],
    "official_guidance": [
      {
        "label": "JOSS submission requirements",
        "url": "https://joss.readthedocs.io/en/latest/submitting.html"
      },
      {
        "label": "JOSS review criteria",
        "url": "https://joss.readthedocs.io/en/latest/review_criteria.html"
      },
      {
        "label": "JOSS review checklist",
        "url": "https://joss.readthedocs.io/en/latest/review_checklist.html"
      }
    ],
    "phases": [
      {
        "deliverables": [
          "Six task-first showcase narratives",
          "Separate synthetic and public qualification evidence",
          "Passing macOS and Linux installation matrix on Python 3.10-3.12",
          "Resolved or explicitly narrowed membrane-orientation and conformational-state claims",
          "Completed third-party license and redistribution audit",
          "Complete AI-use, authorship, funding, and conflict records"
        ],
        "duration": "2-4 weeks",
        "label": "Harden and explain",
        "phase_id": "local_hardening",
        "state": "in_progress"
      },
      {
        "deliverables": [
          "Freely cloneable public repository",
          "Public issue and contribution pathways",
          "Continuous integration and tagged alpha release",
          "Public documentation and static evidence showcase",
          "Changelog, support policy, governance, and security reporting"
        ],
        "duration": "Approval boundary",
        "label": "Open the project",
        "phase_id": "approved_public_launch",
        "state": "awaiting_explicit_approval"
      },
      {
        "deliverables": [
          "Active public development distributed across the period",
          "Tagged releases and documented changes",
          "Independent installation feedback",
          "Documented research use or workflow integration",
          "Public issues, discussions, or external contributions",
          "Software changes traceable to user and scientific feedback"
        ],
        "duration": "More than 6 months",
        "label": "Develop through use",
        "phase_id": "open_development",
        "state": "blocked_public_history_not_started"
      },
      {
        "deliverables": [
          "Stable tagged release",
          "JOSS paper with required sections and verified references",
          "Approved author, affiliation, ORCID, funding, conflict, and AI-use statements",
          "Submission through the JOSS process",
          "Public response to reviewer issues",
          "Final archived software release and DOI when requested"
        ],
        "duration": "After every gate passes",
        "label": "Submit and review",
        "phase_id": "joss_submission",
        "state": "blocked"
      }
    ],
    "policy": "Local preparation does not authorize repository publication, release, archival deposit, or JOSS submission.",
    "product": "YAUVI Structural Biology Platform — Mark 1",
    "schema_version": "1.0",
    "scientific_suite": "YAUVI Structural Workbench",
    "submission_rule": "All scientific, packaging, public-history, research-use, and approval gates must be satisfied before submission."
  },
  "qualification": {
    "cases": [
      {
        "analysis_type": "structure_qc",
        "biological_context": "A structure can proceed to residue-level review only when its identity and provenance are explicit.",
        "case_label": "1CRN plus AlphaFold P69905 v6",
        "failed_checks": [],
        "finding": "Exact 46-residue mapping; raw clashscore 0.0; mean predicted-model pLDDT 98.064; unknown provenance remained incomplete.",
        "independent_reference": "wwPDB validation and AlphaFold DB model-confidence records",
        "remaining_limit": "Coordinate and imported validation evidence do not establish native conformation or function.",
        "source_links": [
          {
            "label": "RCSB and wwPDB files",
            "url": "https://www.rcsb.org/docs/programmatic-access/file-download-services"
          },
          {
            "label": "AlphaFold Protein Structure Database",
            "url": "https://www.alphafold.ebi.ac.uk/"
          }
        ],
        "status": "passed"
      },
      {
        "analysis_type": "membrane_orientation",
        "biological_context": "The current local method is promising for the tested beta-barrels but is not reliable across the tested alpha-helical stratum.",
        "case_label": "Five beta barrels and three alpha-helical OPM structures",
        "failed_checks": [
          "alpha_helical_mean_normal_error",
          "alpha_helical_rotation_invariance"
        ],
        "finding": "Beta-barrel mean normal error 7.442° passed; alpha-helical mean error 31.609° and 1U19 rotation Jaccard 0.88 failed.",
        "independent_reference": "OPM-oriented coordinates; deposited membrane normal is Z and REMARK records bilayer half-thickness",
        "remaining_limit": "Agreement with OPM placement does not establish native-cell exposure or topology in a tested cell.",
        "source_links": [
          {
            "label": "OPM and PPM reference system",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3245162/"
          }
        ],
        "status": "partial"
      },
      {
        "analysis_type": "conformational_state",
        "biological_context": "A resolved resemblance can support state comparison; an unresolved result prevents an unjustified active or inactive claim.",
        "case_label": "ABL active and inactive holdouts",
        "failed_checks": [
          "held_out_8SSN_inactive_like"
        ],
        "finding": "2V7A was active_like at 0.990 Å; 8SSN safely remained unresolved at 6.239 Å.",
        "independent_reference": "KinCore ABL1 experimental-chain labels",
        "remaining_limit": "A state label is geometric resemblance to curated references, not kinase activity.",
        "source_links": [
          {
            "label": "KinCore ABL1 classifications",
            "url": "https://dunbrack.fccc.edu/kincore/GENE/ABL1"
          }
        ],
        "status": "partial"
      },
      {
        "analysis_type": "functional_site_state",
        "biological_context": "Exact residue mapping can support mutational-control planning while keeping annotation separate from observed chemistry.",
        "case_label": "M-CSA glutamate racemase, PDB 1B73",
        "failed_checks": [],
        "finding": "6 curated residues mapped exactly at positions 7, 8, 70, 147, 178, 180; pocket evidence remained missing.",
        "independent_reference": "M-CSA entry 1, glutamate racemase P56868 / PDB 1B73",
        "remaining_limit": "Exact mapping of M-CSA residues is not an observation of catalysis.",
        "source_links": [
          {
            "label": "M-CSA entry 1",
            "url": "https://www.ebi.ac.uk/thornton-srv/m-csa/entry/1/"
          }
        ],
        "status": "passed"
      },
      {
        "analysis_type": "assembly_interface",
        "biological_context": "The result identifies assembly-bound interfaces and burial that are absent from isolated-chain interpretation.",
        "case_label": "Hemoglobin biological assembly 1, PDB 4HHB",
        "failed_checks": [],
        "finding": "4 expected chains recovered; 1760.003 Å² subject surface buried using FreeSASA.",
        "independent_reference": "wwPDB biological assembly 1 for hemoglobin 4HHB",
        "remaining_limit": "A deposited assembly and computed interface are not binding affinity or intact-cell accessibility.",
        "source_links": [
          {
            "label": "RCSB structure 4HHB",
            "url": "https://www.rcsb.org/structure/4HHB"
          },
          {
            "label": "FreeSASA",
            "url": "https://freesasa.github.io/"
          }
        ],
        "status": "passed"
      },
      {
        "analysis_type": "sf_csa",
        "biological_context": "Exact, homologous, analogous, and unresolved relationships remain distinct instead of collapsing into a function score.",
        "case_label": "CATH exact, homolog, fold-analogy, and unrelated controls",
        "failed_checks": [],
        "finding": "3 structure hits and 2 sequence hits were kept separate; Foldseek 10.941cd33 and diamond version 2.1.11 ran locally.",
        "independent_reference": "CATH v4.3 classifications plus public PDB controls",
        "remaining_limit": "A fold or sequence hit cannot be promoted to exact functional transfer.",
        "source_links": [
          {
            "label": "CATH downloads",
            "url": "https://cathdb.info/download"
          },
          {
            "label": "Foldseek",
            "url": "https://github.com/steineggerlab/foldseek"
          },
          {
            "label": "DIAMOND",
            "url": "https://www.nature.com/articles/s41592-021-01101-x"
          }
        ],
        "status": "passed"
      }
    ],
    "collection_id": "yauvi-structural-public-qualification-v1",
    "files": [
      {
        "label": "Readable qualification report",
        "path": "qualification/QUALIFICATION_REPORT.md"
      },
      {
        "label": "Machine-readable qualification results",
        "path": "qualification/QUALIFICATION_RESULTS.json"
      },
      {
        "label": "Source checksum lock",
        "path": "qualification/SOURCE_LOCK.json"
      },
      {
        "label": "Source verification record",
        "path": "qualification/SOURCE_VERIFICATION.json"
      }
    ],
    "overall_state": "incomplete_or_failed",
    "qualification_rule": "All required checks for all six workflows must pass. Software tests are a separate gate.",
    "result_sha256": "06c7e84e10e45335d506554d0587ea4979ca8c0b113c4c4e2d0a9c015d3f8eed",
    "source_artifact_count": 27,
    "source_lock_passed": true,
    "source_verification_sha256": "48091cfb8fc1385d751d3226c503f96b7312f3cfcfbbae596724629dd561f3d1",
    "workflow_counts": {
      "blocked": 0,
      "failed": 0,
      "partial": 2,
      "passed": 4
    }
  },
  "qualification_v2": {
    "collection_id": "yauvi-structural-public-qualification-v2",
    "files": [
      {
        "label": "Qualification v2 panel specification",
        "path": "qualification-v2/PANEL_MANIFEST.json"
      },
      {
        "label": "Qualification v2 status",
        "path": "qualification-v2/QUALIFICATION_V2_STATUS.json"
      },
      {
        "label": "Printable v2 audit",
        "path": "qualification-v2/QUALIFICATION_REPORT.html"
      },
      {
        "label": "V2 stratum gaps",
        "path": "qualification-v2/STRATUM_STATUS.tsv"
      },
      {
        "label": "V2 execution summary",
        "path": "qualification-v2/EXECUTION_SUMMARY.json"
      }
    ],
    "missing_records": 50,
    "overall_state": "blocked_panel_incomplete",
    "panels": [
      {
        "execution_state": "failed",
        "missing_count": 16,
        "record_count": 16,
        "state": "blocked_panel_incomplete",
        "workflow": "membrane_orientation"
      },
      {
        "execution_state": "not_executed",
        "missing_count": 18,
        "record_count": 0,
        "state": "blocked_panel_incomplete",
        "workflow": "conformational_state"
      },
      {
        "execution_state": "passed",
        "missing_count": 0,
        "record_count": 16,
        "state": "ready_for_execution",
        "workflow": "structure_qc"
      },
      {
        "execution_state": "passed",
        "missing_count": 0,
        "record_count": 16,
        "state": "ready_for_execution",
        "workflow": "functional_site_state"
      },
      {
        "execution_state": "passed",
        "missing_count": 0,
        "record_count": 16,
        "state": "ready_for_execution",
        "workflow": "assembly_interface"
      },
      {
        "execution_state": "not_executed",
        "missing_count": 16,
        "record_count": 0,
        "state": "blocked_panel_incomplete",
        "workflow": "sf_csa"
      }
    ],
    "scientific_execution": {
      "all_release_blocking_scopes_qualified": false,
      "cases_passed": 53,
      "cases_required": 114,
      "collection_note": "Three panels are fully composed and one is half composed. StructQC 16/16 with 2 controls and 7/7 coverage; site-context 16/16 with 1 control and 8/8; assembly-context 16/16 with 6/6; membrane 14/16 of the adopted beta_barrel stratum (16 of the panel's 32) with 6/6 coverage. 1BXW fails accuracy at 6.197 degrees against the collection 2.1 bound of 1.0, and 1QD6 fails at 8.391 degrees max with 8.391 degrees of rotational drift once collection 2.2 raised the rotation count from 5 to 8. 1QD6 passed at five rotations only because that sample never drew the second basin; arm64 and ubuntu x64 both reported the identical 8.390785 degrees at five. See PANEL_MANIFEST.json threshold_revisions. This flag stays false because it reports the whole collection: membrane's alpha_helical stratum, ABL StateAtlas and sf-csa remain unadopted, and no scope has completed the independent second-machine gate.",
      "counts_are_single_machine": true,
      "every_executed_panel_passed": false,
      "panels": [
        {
          "cases_adopted": 16,
          "cases_passed": 16,
          "cases_required": 16,
          "controls_passed": 0,
          "controls_total": 0,
          "coverage_required": 6,
          "coverage_unwitnessable": [],
          "coverage_witnessed": 6,
          "stratum_scope": "all_four_assembly_strata",
          "stratum_state": "passed",
          "workflow": "assembly_interface"
        },
        {
          "cases_adopted": 16,
          "cases_passed": 5,
          "cases_required": 32,
          "controls_passed": 0,
          "controls_total": 0,
          "coverage_required": 6,
          "coverage_unwitnessable": [],
          "coverage_witnessed": 6,
          "stratum_scope": "beta_barrel",
          "stratum_state": "failed",
          "workflow": "membrane_orientation"
        },
        {
          "cases_adopted": 16,
          "cases_passed": 16,
          "cases_required": 16,
          "controls_passed": 1,
          "controls_total": 1,
          "coverage_required": 8,
          "coverage_unwitnessable": [
            "curated_residue_missing_coordinates"
          ],
          "coverage_witnessed": 8,
          "stratum_scope": "all_four_site_context_strata",
          "stratum_state": "passed",
          "workflow": "functional_site_state"
        },
        {
          "cases_adopted": 16,
          "cases_passed": 16,
          "cases_required": 16,
          "controls_passed": 2,
          "controls_total": 2,
          "coverage_required": 7,
          "coverage_unwitnessable": [
            "modified_residues"
          ],
          "coverage_witnessed": 7,
          "stratum_scope": "all_four_structqc_strata",
          "stratum_state": "passed",
          "workflow": "structure_qc"
        }
      ],
      "panels_executed": 4,
      "panels_total": 6,
      "recorded_on": [
        {
          "machine": "x86_64",
          "platform": "Darwin",
          "python": "3.12.7"
        }
      ],
      "scope_qualification_note": "Executed panels passing is not scope qualification. A Mark 1 scope is qualified only when its panel composes in full, every case and control passes, and the result reproduces independently on a second machine. This summary reports execution only. 2 of 6 panels are unadopted and no second-machine reproduction is recorded, so no scope is qualified.",
      "second_machine_reproduction": "not_recorded",
      "workflows_executed": [
        "assembly_interface",
        "functional_site_state",
        "membrane_orientation",
        "structure_qc"
      ],
      "workflows_not_executed": [
        "conformational_state",
        "sf_csa"
      ]
    },
    "scientific_execution_performed": true
  },
  "release": {
    "publication_authorized": false,
    "release_candidate": false,
    "release_state": "pre_public_preparation",
    "submission_eligible": false,
    "version_control": "public_git_repository"
  },
  "reviewer_files": [
    {
      "label": "Reviewer quickstart",
      "path": "reviewer/REVIEWER_QUICKSTART.md"
    },
    {
      "label": "Publication roadmap",
      "path": "reviewer/JOSS_PUBLICATION_ROADMAP.json"
    },
    {
      "label": "JOSS paper preview",
      "path": "reviewer/PAPER_PREVIEW.md"
    },
    {
      "label": "JOSS preparation checklist",
      "path": "reviewer/JOSS_CHECKLIST.md"
    }
  ],
  "schema_version": "1.0",
  "showcase_id": "yauvi-public-evidence-showcase",
  "workflows": [
    {
      "analysis_type": "structure_qc",
      "external_benchmark": "public_case_passed",
      "external_benchmark_detail": "v2_panel_fully_composed_all_four_strata_executed_and_passed_second_machine_gate_outstanding",
      "inputs": [
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "PDB or mmCIF coordinates",
          "required": true,
          "role": "structure",
          "source_ids": [
            "pdb",
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "completeness_unevaluated",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Reference sequence",
          "required": false,
          "role": "reference_fasta",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "provenance_unknown",
          "extensions": [
            ".json"
          ],
          "label": "Provenance declaration",
          "required": false,
          "role": "provenance",
          "source_ids": []
        },
        {
          "absence_effect": "domain_confidence_unevaluated",
          "extensions": [
            ".json"
          ],
          "label": "Predicted aligned error",
          "required": false,
          "role": "pae",
          "source_ids": [
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "scientifically_incomplete",
          "extensions": [
            ".xml",
            ".json"
          ],
          "label": "wwPDB or local validation report",
          "required": false,
          "role": "validation_report",
          "source_ids": [
            "wwpdb_validation"
          ]
        }
      ],
      "measures": "Coordinate identity, models, chains, numbering, missing backbone atoms, provenance, confidence encoding, reference-sequence mapping, PAE, and imported validation evidence.",
      "non_claim": "It does not establish the native biological conformation, protein function, or experimental activity.",
      "optional_runtimes": {
        "gemmi": "available",
        "mkdssp": "missing",
        "molprobity_or_phenix": "missing"
      },
      "public_question": "Can I trust these coordinates?",
      "question": "Are these coordinates identity-bound, provenance-declared, and suitable for interpretation?",
      "required_runtimes": {},
      "scientific_scopes": [
        {
          "benchmark_collection": "qualification-v2-structqc",
          "known_limitations": [
            "Qualification v2 public panel is not complete."
          ],
          "release_blocking": true,
          "required_evidence": [
            "exact coordinates",
            "explicit provenance",
            "method-appropriate validation evidence"
          ],
          "scientific_state": "prototype",
          "scope_id": "coordinate_provenance_and_validation",
          "supported_subject_class": "experimental and predicted protein coordinate models"
        }
      ],
      "showcase_note": "Executed synthetic evidence is linked below.",
      "showcase_state": "passed_synthetic_case",
      "software_state": "prototype",
      "title": "Structure QC"
    },
    {
      "analysis_type": "membrane_orientation",
      "external_benchmark": "partial_public_case",
      "external_benchmark_detail": "NON_BLOCKING_from_collection_2_4_research_only_beta_barrel_executes_at_5_of_16_under_the_2_3_accuracy_gate_eleven_cases_sit_2_5_to_17_4_deg_from_OPM_mark_1_makes_no_accuracy_claim_for_membrane_orientation",
      "inputs": [
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "PDB or mmCIF coordinates",
          "required": true,
          "role": "structure",
          "source_ids": [
            "pdb",
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "alpha_helical_placement_incomplete",
          "extensions": [
            ".json"
          ],
          "label": "Coordinate-bound transmembrane spans",
          "required": false,
          "role": "topology_evidence",
          "source_ids": [
            "opm_ppm",
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "completeness_unevaluated",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Reference sequence",
          "required": false,
          "role": "reference_fasta",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "provenance_unknown",
          "extensions": [
            ".json"
          ],
          "label": "Provenance declaration",
          "required": false,
          "role": "provenance",
          "source_ids": []
        },
        {
          "absence_effect": "scientifically_incomplete",
          "extensions": [
            ".xml",
            ".json"
          ],
          "label": "wwPDB or local validation report",
          "required": false,
          "role": "validation_report",
          "source_ids": [
            "wwpdb_validation"
          ]
        }
      ],
      "measures": "A context-declared membrane placement, orientation transform, membrane depth, and modeled residue accessibility.",
      "non_claim": "Modeled orientation is not direct evidence of intact-cell exposure or topology in the tested organism.",
      "optional_runtimes": {},
      "public_question": "How might this protein sit in a membrane?",
      "question": "How does this protein sit in its declared membrane or surface context?",
      "required_runtimes": {},
      "scientific_scopes": [
        {
          "benchmark_collection": "qualification-v2-membrane-beta-barrel",
          "known_limitations": [
            "Independent second-machine reproduction remains required."
          ],
          "release_blocking": true,
          "required_evidence": [
            "exact coordinates",
            "declared membrane context"
          ],
          "scientific_state": "conditionally_qualified",
          "scope_id": "beta_barrel",
          "supported_subject_class": "transmembrane beta-barrel proteins"
        },
        {
          "benchmark_collection": "qualification-v2-membrane-alpha-helical",
          "known_limitations": [
            "Experimental method; not part of the Mark 1 qualified scope."
          ],
          "release_blocking": false,
          "required_evidence": [
            "exact coordinates",
            "checksum-bound transmembrane spans",
            "external topology for sidedness"
          ],
          "scientific_state": "prototype",
          "scope_id": "alpha_helical",
          "supported_subject_class": "single-pass and multipass alpha-helical membrane proteins"
        }
      ],
      "showcase_note": "Executed synthetic evidence is linked below.",
      "showcase_state": "passed_synthetic_case",
      "software_state": "conditionally_qualified",
      "title": "Membrane orientation"
    },
    {
      "analysis_type": "conformational_state",
      "external_benchmark": "partial_public_case",
      "external_benchmark_detail": "abl_exact_mapping_v2_implemented_held_out_panel_unadopted",
      "inputs": [
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "Query structure or trajectory topology",
          "required": true,
          "role": "structure",
          "source_ids": [
            "pdb",
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "single_structure_only",
          "extensions": [
            ".xtc",
            ".dcd",
            ".trr"
          ],
          "label": "Optional trajectory",
          "required": false,
          "role": "trajectory",
          "source_ids": []
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "Experimental active-state references",
          "required": true,
          "role": "active_reference",
          "source_ids": [
            "pdb",
            "sifts"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "Experimental inactive-state references",
          "required": true,
          "role": "inactive_reference",
          "source_ids": [
            "pdb",
            "sifts"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".json"
          ],
          "label": "Exact SIFTS or residue-equivalence map",
          "required": true,
          "role": "alignment_map",
          "source_ids": [
            "sifts",
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "completeness_unevaluated",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Reference sequence",
          "required": false,
          "role": "reference_fasta",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "provenance_unknown",
          "extensions": [
            ".json"
          ],
          "label": "Query provenance declaration",
          "required": false,
          "role": "provenance",
          "source_ids": []
        },
        {
          "absence_effect": "scientifically_incomplete",
          "extensions": [
            ".xml",
            ".json"
          ],
          "label": "wwPDB or local validation report",
          "required": false,
          "role": "validation_report",
          "source_ids": [
            "wwpdb_validation"
          ]
        }
      ],
      "measures": "Sequence-mapped alignment, RMSD, RMSF, frame-to-reference distance, deterministic clusters, and interpretable-frame populations.",
      "non_claim": "Active-like or inactive-like describes structural resemblance, not biochemical activation, inhibition, or efficacy.",
      "optional_runtimes": {
        "mdanalysis": "available"
      },
      "public_question": "Which conformation does it resemble?",
      "question": "Which experimentally bounded conformations does this structure or ensemble resemble?",
      "required_runtimes": {},
      "scientific_scopes": [
        {
          "benchmark_collection": "qualification-v2-abl-state-atlas",
          "known_limitations": [
            "Qualification v2 held-out ABL gate has not yet passed."
          ],
          "release_blocking": true,
          "required_evidence": [
            "two-sided experimental references",
            "exact SIFTS or explicit residue map",
            "at least 90 percent mapped ABL1 242-495 coverage"
          ],
          "scientific_state": "prototype",
          "scope_id": "abl_family",
          "supported_subject_class": "ABL-family experimental coordinate structures and ensembles"
        },
        {
          "benchmark_collection": "none",
          "known_limitations": [
            "Outside the Mark 1 qualified scope."
          ],
          "release_blocking": false,
          "required_evidence": [
            "two-sided experimental references",
            "exact declared alignment map"
          ],
          "scientific_state": "prototype",
          "scope_id": "other_proteins",
          "supported_subject_class": "other proteins with user-curated two-sided references"
        }
      ],
      "showcase_note": "Executed synthetic evidence is linked below.",
      "showcase_state": "passed_synthetic_case",
      "software_state": "prototype",
      "title": "Conformational resemblance"
    },
    {
      "analysis_type": "functional_site_state",
      "external_benchmark": "public_case_passed",
      "external_benchmark_detail": "v2_panel_fully_composed_all_four_strata_executed_and_passed_second_machine_gate_outstanding",
      "inputs": [
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "PDB or mmCIF coordinates",
          "required": true,
          "role": "structure",
          "source_ids": [
            "pdb",
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "completeness_unevaluated",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Reference sequence",
          "required": true,
          "role": "reference_fasta",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".json",
            ".tsv",
            ".csv"
          ],
          "label": "Site annotations (JSON, TSV, or CSV)",
          "required": true,
          "role": "site_annotations",
          "source_ids": [
            "mcsa",
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "activity_annotation_leg_missing",
          "extensions": [
            ".tsv",
            ".csv"
          ],
          "label": "Optional UniProt feature export for ActState",
          "required": false,
          "role": "uniprot_annotations",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "cofactor_identity_unresolved",
          "extensions": [
            ".json"
          ],
          "label": "Exact CCD/ChEBI component map",
          "required": false,
          "role": "component_map",
          "source_ids": [
            "pdb_ccd",
            "chebi"
          ]
        },
        {
          "absence_effect": "pocket_evidence_not_run",
          "extensions": [
            ".json"
          ],
          "label": "Method-declared pocket result",
          "required": false,
          "role": "pocket_result",
          "source_ids": []
        },
        {
          "absence_effect": "provenance_unknown",
          "extensions": [
            ".json"
          ],
          "label": "Provenance declaration",
          "required": false,
          "role": "provenance",
          "source_ids": []
        },
        {
          "absence_effect": "scientifically_incomplete",
          "extensions": [
            ".xml",
            ".json"
          ],
          "label": "wwPDB or local validation report",
          "required": false,
          "role": "validation_report",
          "source_ids": [
            "wwpdb_validation"
          ]
        }
      ],
      "measures": "Exact residue mapping, role-specific site completeness, ligand/cofactor observations, geometry, and method-specific pocket evidence.",
      "non_claim": "A complete or plausible site is not observed catalysis, binding, inhibition, or physiological function.",
      "optional_runtimes": {
        "gemmi": "available"
      },
      "public_question": "Are important functional residues present?",
      "question": "Are declared functional residues mapped, chemically plausible, and observed in this coordinate state?",
      "required_runtimes": {},
      "scientific_scopes": [
        {
          "benchmark_collection": "qualification-v2-site-context",
          "known_limitations": [
            "Qualification v2 M-CSA panel is not complete."
          ],
          "release_blocking": true,
          "required_evidence": [
            "exact residue mapping",
            "role-specific annotations",
            "exact component identifiers when applicable"
          ],
          "scientific_state": "prototype",
          "scope_id": "curated_functional_site_mapping",
          "supported_subject_class": "proteins with curated residue-role annotations"
        }
      ],
      "showcase_note": "Executed synthetic evidence is linked below.",
      "showcase_state": "passed_synthetic_case",
      "software_state": "prototype",
      "title": "Functional-site evidence"
    },
    {
      "analysis_type": "assembly_interface",
      "external_benchmark": "public_case_passed",
      "external_benchmark_detail": "v2_panel_fully_composed_all_four_strata_executed_and_passed_second_machine_gate_outstanding",
      "inputs": [
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "Isolated subject coordinates",
          "required": true,
          "role": "structure",
          "source_ids": [
            "pdb",
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "Expanded biological assembly",
          "required": true,
          "role": "assembly",
          "source_ids": [
            "pdb"
          ]
        },
        {
          "absence_effect": "completeness_unevaluated",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Reference sequence",
          "required": false,
          "role": "reference_fasta",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "provenance_unknown",
          "extensions": [
            ".json"
          ],
          "label": "Provenance declaration",
          "required": false,
          "role": "provenance",
          "source_ids": []
        },
        {
          "absence_effect": "scientifically_incomplete",
          "extensions": [
            ".xml",
            ".json"
          ],
          "label": "wwPDB or local validation report",
          "required": false,
          "role": "validation_report",
          "source_ids": [
            "wwpdb_validation"
          ]
        }
      ],
      "measures": "Assembly identity, stoichiometry evidence, heavy-atom contacts, interface residues, SASA, burial, and lower-bound status.",
      "non_claim": "Assembly geometry does not establish native abundance, intact-cell accessibility, or measured binding.",
      "optional_runtimes": {
        "freesasa": "## FreeSASA 2.1.3 ##",
        "gemmi": "available"
      },
      "public_question": "Which residues form an assembly interface?",
      "question": "Which residues contact or become buried in a declared biological assembly?",
      "required_runtimes": {},
      "scientific_scopes": [
        {
          "benchmark_collection": "qualification-v2-assembly-context",
          "known_limitations": [
            "Qualification v2 operator and FreeSASA panel is not complete."
          ],
          "release_blocking": true,
          "required_evidence": [
            "exact assembly operators",
            "entity and copy identity",
            "pinned FreeSASA runtime for qualification"
          ],
          "scientific_state": "prototype",
          "scope_id": "deposited_biological_assembly",
          "supported_subject_class": "deposited protein biological assemblies"
        }
      ],
      "showcase_note": "Executed synthetic evidence is linked below.",
      "showcase_state": "passed_synthetic_case",
      "software_state": "prototype",
      "title": "Assembly and interfaces"
    },
    {
      "analysis_type": "sf_csa",
      "external_benchmark": "public_case_passed",
      "external_benchmark_detail": "qualification_v1_public_case_passed_real_engines_v2_four_family_panel_unadopted",
      "inputs": [
        {
          "absence_effect": "blocked",
          "extensions": [
            ".pdb",
            ".cif",
            ".mmcif"
          ],
          "label": "Query structure",
          "required": true,
          "role": "query_structure",
          "source_ids": [
            "pdb",
            "alphafold_db"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Query sequence FASTA",
          "required": true,
          "role": "query_fasta",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".fasta",
            ".fa",
            ".faa"
          ],
          "label": "Local comparison proteome",
          "required": true,
          "role": "source_proteome",
          "source_ids": [
            "uniprot_proteomes"
          ]
        },
        {
          "absence_effect": "blocked",
          "extensions": [
            ".json"
          ],
          "label": "Organism-appropriate interpretation tables",
          "required": true,
          "role": "interpretation_tables",
          "source_ids": [
            "mcsa"
          ]
        },
        {
          "absence_effect": "provenance_unknown",
          "extensions": [
            ".json"
          ],
          "label": "Query provenance declaration",
          "required": false,
          "role": "provenance",
          "source_ids": []
        },
        {
          "absence_effect": "scientifically_incomplete",
          "extensions": [
            ".xml",
            ".json"
          ],
          "label": "wwPDB or local validation report",
          "required": false,
          "role": "validation_report",
          "source_ids": [
            "wwpdb_validation"
          ]
        }
      ],
      "measures": "Foldseek structural hits, DIAMOND sequence hits, declared mechanism-group context, divergence, and closed-vocabulary relationship evidence.",
      "non_claim": "Similarity, homology, or a shared fold never becomes exact identity or automatic functional transfer.",
      "optional_runtimes": {},
      "public_question": "How does this protein relate to other proteins?",
      "question": "How do structural similarity and sequence homology compare without collapsing them into one claim?",
      "required_runtimes": {
        "diamond": "diamond version 2.1.11",
        "foldseek": "10.941cd33"
      },
      "scientific_scopes": [
        {
          "benchmark_collection": "qualification-v2-sf-csa",
          "known_limitations": [
            "Qualification v2 four-family panel is not complete."
          ],
          "release_blocking": true,
          "required_evidence": [
            "checksum-pinned Foldseek database",
            "declared proteome",
            "closed interpretation tables",
            "Foldseek and DIAMOND runtimes"
          ],
          "scientific_state": "prototype",
          "scope_id": "curated_structure_sequence_comparison",
          "supported_subject_class": "curator-declared protein family comparison panels"
        }
      ],
      "showcase_note": "Canonical pipeline and release audit passed through deterministic test doubles. Foldseek and DIAMOND computed no alignments; external binary benchmarks remain pending.",
      "showcase_state": "passed_stubbed_pipeline_case",
      "software_state": "prototype",
      "title": "SF-CSA comparison"
    }
  ]
};

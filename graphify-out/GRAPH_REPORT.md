# Graph Report - .  (2026-06-26)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 484 nodes · 956 edges · 17 communities (14 shown, 3 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 67 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]

## God Nodes (most connected - your core abstractions)
1. `Grna` - 77 edges
2. `GrnaType` - 40 edges
3. `DesignResult` - 39 edges
4. `ScoreResult` - 35 edges
5. `compute_all_scores()` - 32 edges
6. `OffTargetHit` - 29 edges
7. `Genome` - 28 edges
8. `design_guides()` - 20 edges
9. `find_off_targets_brute_force()` - 19 edges
10. `find_off_targets()` - 16 edges

## Surprising Connections (you probably didn't know these)
- `TestComputeSpecificity` --uses--> `DesignResult`  [INFERRED]
  tests/test_design.py → crispr_design/core/design.py
- `TestDesignGuides` --uses--> `DesignResult`  [INFERRED]
  tests/test_design.py → crispr_design/core/design.py
- `TestDesignResult` --uses--> `DesignResult`  [INFERRED]
  tests/test_design.py → crispr_design/core/design.py
- `TestFilterGuides` --uses--> `DesignResult`  [INFERRED]
  tests/test_design.py → crispr_design/core/design.py
- `TestComputeSpecificity` --uses--> `GrnaType`  [INFERRED]
  tests/test_design.py → crispr_design/core/grna.py

## Import Cycles
- None detected.

## Communities (17 total, 3 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.06
Nodes (24): Grna, Return reverse complement of the guide sequence., Calculate GC content of the guide sequence., GC content as percentage., Full gRNA sequence including PAM., Serialize to dictionary., CRISPR guide RNA with PAM detection and scoring.      Attributes:         sequen, Validate the gRNA sequence and PAM. (+16 more)

### Community 1 - "Community 1"
Cohesion: 0.06
Nodes (41): Any, DesignResult, Result of a gRNA design scan., Number of off-target sites found., batch(), design(), main(), CRISPR Design Studio CLI — main entry point. (+33 more)

### Community 2 - "Community 2"
Cohesion: 0.06
Nodes (23): compute_all_scores(), Result from scoring a gRNA or off-target site., Compute MIT-CF score (Cong et al. 2013).      The score is the sum of position-s, Compute CFD score (Hsu et al. 2013).      The CFD score accounts for position-sp, Compute Doench 2014 rule set 1 score.      Predicts on-target activity for SpCas, Compute Doench 2016 rule set 2 score.      Improved prediction of on-target acti, Compute all scoring algorithms for a gRNA.      Args:         sequence: 20bp gRN, score_cfd() (+15 more)

### Community 3 - "Community 3"
Cohesion: 0.06
Nodes (25): Tests for sequence utilities., TestExtractGrnaFromFasta, TestFormatSequence, TestGCContent, TestHammingDistance, TestParseFasta, TestReverseComplement, TestValidateDNASequence (+17 more)

### Community 4 - "Community 4"
Cohesion: 0.06
Nodes (42): batch_endpoint(), BatchInput, design_endpoint(), _get_genome_dict(), _get_test_genome(), GrnaInput, health(), list_genomes() (+34 more)

### Community 5 - "Community 5"
Cohesion: 0.07
Nodes (20): Upload a custom genome FASTA file., upload_genome(), create_test_genome(), Genome, Return reverse complement of a DNA sequence., Convert entire genome to dictionary (memory intensive!).          Returns:, Get genomic context around a position.          Args:             chrom: Chromos, Close genome file handles. (+12 more)

### Community 6 - "Community 6"
Cohesion: 0.08
Nodes (13): detect_pam_type(), find_pams_in_sequence(), Create a Grna from a full target sequence (guide + PAM).          The target_seq, Validate PAM sequence for the given gRNA type., Find all PAM sites in a DNA sequence and extract potential gRNAs.      Args:, Detect which PAM types are present in a sequence.      Args:         sequence: D, Tests for gRNA class and PAM detection., Create guide from full target sequence (guide + PAM). (+5 more)

### Community 7 - "Community 7"
Cohesion: 0.08
Nodes (7): Tests for FastAPI server., TestBatchEndpoint, TestDesignEndpoint, TestGenomesEndpoint, TestHealthEndpoint, TestScanEndpoint, TestScoreEndpoint

### Community 8 - "Community 8"
Cohesion: 0.13
Nodes (11): _compute_specificity(), design_guides(), filter_guides(), Compute a composite specificity score.      Combines on-target activity predicti, Filter design results based on quality criteria.      Args:         results: Des, Design and rank gRNAs from a DNA sequence.      Scans the sequence for all valid, Tests for guide design module., Design guides from a known sequence with GG PAMs. (+3 more)

### Community 9 - "Community 9"
Cohesion: 0.09
Nodes (8): Tests for CLI interface., Test designing from a raw sequence string., Test that missing sequence is handled., TestCLIBatch, TestCLIDesign, TestCLIHelp, TestCLIReport, TestCLIScan

### Community 10 - "Community 10"
Cohesion: 0.10
Nodes (19): design_result(), guide_list_file(), off_target_hit(), Shared test fixtures., Create a file with multiple guide sequences for batch testing., A known DNA sequence with GG PAMs., A 20bp guide sequence., Create a small test genome for testing. (+11 more)

## Knowledge Gaps
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Grna` connect `Community 0` to `Community 1`, `Community 2`, `Community 4`, `Community 6`, `Community 8`, `Community 10`?**
  _High betweenness centrality (0.237) - this node is a cross-community bridge._
- **Why does `parse_fasta()` connect `Community 3` to `Community 1`?**
  _High betweenness centrality (0.193) - this node is a cross-community bridge._
- **Why does `Genome` connect `Community 5` to `Community 1`, `Community 4`?**
  _High betweenness centrality (0.106) - this node is a cross-community bridge._
- **Are the 18 inferred relationships involving `Grna` (e.g. with `DesignResult` and `OffTargetHit`) actually correct?**
  _`Grna` has 18 INFERRED edges - model-reasoned connections that need verification._
- **Are the 18 inferred relationships involving `GrnaType` (e.g. with `DesignResult` and `OffTargetHit`) actually correct?**
  _`GrnaType` has 18 INFERRED edges - model-reasoned connections that need verification._
- **Are the 11 inferred relationships involving `DesignResult` (e.g. with `Grna` and `GrnaType`) actually correct?**
  _`DesignResult` has 11 INFERRED edges - model-reasoned connections that need verification._
- **Are the 14 inferred relationships involving `ScoreResult` (e.g. with `DesignResult` and `TestComputeSpecificity`) actually correct?**
  _`ScoreResult` has 14 INFERRED edges - model-reasoned connections that need verification._
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

def test_internal_rbs_checker_detects_rbs():
    checker = InternalRBSChecker()
    checker.initiate()

    seq = "ATGAAACCCAGGAGGTTTCCC"
    ok, bad = checker.run(seq)

    assert ok is False
    assert bad is not None

def test_internal_rbs_checker_no_rbs():
    checker = InternalRBSChecker()
    checker.initiate()

    seq = "ATGAAACCCAGGTTTCCC"
    ok, bad = checker.run(seq)

    assert ok is True
    assert bad is None

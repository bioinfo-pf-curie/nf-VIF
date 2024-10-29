process BUILD_TEST_BLAT_DB{
    container "phinguyen2000/blat:255336f"

    input:
    path(fasta)

    output:
    path("test_db.2bit")

    """
    faToTwoBit $fasta test_db.2bit
    """
}
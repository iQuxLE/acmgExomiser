package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.Acmg2015EvidenceAssigner;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.*;
import org.monarchinitiative.exomiser.core.model.ChromosomalRegionIndex;
import org.monarchinitiative.exomiser.core.model.TranscriptAnnotation;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.is;
import static org.junit.jupiter.api.Assertions.*;
import static org.monarchinitiative.exomiser.core.model.Pedigree.Individual.Sex.MALE;
import static org.monarchinitiative.exomiser.core.model.Pedigree.justProband;

class PS1PM5AssignerTest {

    private final AlleleProto.ClinVar clinVarPathogenicStarRating2 = AlleleProto.ClinVar.newBuilder()
            .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
            .setReviewStatus("criteria provided, multiple submitters, no conflicts")
            .build();

    private final AlleleProto.ClinVar clinVarPathogenicStarRating1 = AlleleProto.ClinVar.newBuilder()
            .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
            .setReviewStatus("criteria provided, single submitter")
            .build();

    private final AlleleProto.ClinVar clinVarBenignStarRating2 = AlleleProto.ClinVar.newBuilder()
            .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.BENIGN)
            .setReviewStatus("criteria provided, multiple submitters, no conflicts")
            .build();


    private MVStore mvStore(){
        return new MVStore.Builder().compress().open();
    }

    private VariantDataService initializeCustomVariantDataservice(AlleleProto.ClinVar cvInterpretation, AlleleProto.AlleleKey... alleleKeys){
        TestVariantDataService.Builder builder = TestVariantDataService.builder()
                .setGenomeAssembly(GenomeAssembly.HG19);
        for (AlleleProto.AlleleKey alleleKey : alleleKeys){
            builder.put(alleleKey,cvInterpretation);
        }
        return builder.build();
    }

    private VariantDataService initializeVariantDataservice(){
        TestVariantDataService.Builder builder = TestVariantDataService.builder()
                .setGenomeAssembly(GenomeAssembly.HG19);

        return builder.build();
    }

    private VariantDataService variantDataServiceForGeneStats(String geneSymbol, ClinVarGeneStats stats){
        TestVariantDataService.Builder builder = TestVariantDataService.builder()
                .setGenomeAssembly(GenomeAssembly.HG19);
        builder.put(geneSymbol, stats);
        return builder.build();
    }

    private VariantEvaluation buildVariantEvaluation(int chr, int pos, String ref, String alt, String hgvs, String cdna, String geneSymbol, VariantEffect variantEffect) {
        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(variantEffect)
                .hgvsProtein(hgvs)
                .hgvsCdna(cdna)
                .build();

        return TestFactory.variantBuilder(chr, pos, ref, alt)
                .geneSymbol(geneSymbol)
                .variantEffect(variantEffect)
                .annotations(List.of(transcriptAnnotation))
                .build();
    }

    private AlleleProto.AlleleKey parseAlleleKey(String variant){
        String[] vars = variant.split("-");
        return AlleleProto.AlleleKey.newBuilder()
                .setChr(Integer.parseInt(vars[0]))
                .setPosition(Integer.parseInt(vars[1]))
                .setRef(vars[2])
                .setAlt(vars[3])
                .build();
    }

    @Test
    void throwsExceptionWithMismatchedIds() {
        assertThrows(IllegalArgumentException.class, () -> new Acmg2015EvidenceAssigner("Zaphod", justProband("Ford", MALE), null, null));
    }

    private final JannovarVariantAnnotator jannovarAnnotator = new JannovarVariantAnnotator(TestFactory.getDefaultGenomeAssembly(), TestFactory
            .buildDefaultJannovarData(), ChromosomalRegionIndex.empty());


    @ParameterizedTest
    @CsvSource({
            "10-123276892-C-G, true, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/374820/ diff cdna, same AAchange
            "10-123276892-CG-AT, false, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/374820/ (no snp))
            "10-123276892-CC-GG, false, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/374820/ mock of upper to test vs ref allele length > 1 (no snp)
            "10-123276893-A-T, false, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/13267/ same cdna + same AAchange
            "10-123247514-C-A, false, false", // https://www.ncbi.nlm.nih.gov/clinvar/variation/1698211/ wrong codon
            "10-123247514-C-G, false, false", // https://www.ncbi.nlm.nih.gov/clinvar/variation/661397/ wrong codon
            "10-299372-G-C, false, false", // https://www.ncbi.nlm.nih.gov/clinvar/variation/689498/ wrong chromosome + wrong codon
            "11-123276893-A-T, false, false", // mocked Variant - wrong chromosome
            //           "10-123276893-AC-CT, false, true" // https://www.ncbi.nlm.nih.gov/clinvar/variation/2131381/ cDna c.1023-1024 // jannovar doesnt annotate to Missense but MNV

    })
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/13267/
    void testAssignPS1orPM5againstChr10Pos123276893AT(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating2, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123276893, "A", "T",
                "p.(Cys342Ser)", "c.1024T>A", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }

    @ParameterizedTest
    @CsvSource({
            "10-123247514-C-A, false, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/1698211/ no match for PS1 cause same cDNA and same proteinChange
            "10-123247514-C-G, true, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/661397/ match for PS1 cause different cDNA and same proteinChange
            "10-123247515-T-C, false, true" // mocked Variant - match PM5 cause different/new proteinChange - does not hit PS1 cause of same reason

    })

    void testAssignPS1orPM5againstChr10Pos123247514CA(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating2, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }


    @ParameterizedTest
    @CsvSource({
            // VariantEvaluation is not MISSENSE so no matches
            "10-123247514-C-A, false, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/1698211/
            "10-123247514-C-G, false, false",  // https://www.ncbi.nlm.nih.gov/clinvar/variation/661397/
            "10-123247515-T-C, false, false" // mocked Variant
    })

    void testAssignPS1orPM5VariantEvaluationIsNotMissense(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating2, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.START_LOST); //mocked

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }

    @ParameterizedTest
    @CsvSource({
            // VariantEvaluation would fit but clinVarMap holds mocked Variant with Benign interpretation only in this Range
            "10-123247515-T-C, false, false" // mocked Variant (here made not pathogenic)

    })

    void testAssignPS1orPM5NotPathogenic(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarBenignStarRating2, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }

    @ParameterizedTest
    @CsvSource({
            // VariantEvaluation would fit but clinVarMap holds mocked Variant with starRating1 only in this Range
            "10-123247515-T-C, false, false" // mocked Variant (here made with no fitting starRating)
    })

    void testAssignPS1orPM5NoStarRating1(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating1, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }

    @ParameterizedTest
    @CsvSource({
            "10-123276856-G-A, false, true", // https://www.ncbi.nlm.nih.gov/clinvar/variation/1066784/
            "10-123276856-G-C, false, false" // https://www.ncbi.nlm.nih.gov/clinvar/variation/13265/ no hit (same AAchange
    })

    void testAssignPS1orPM5againstChr10Pos123276856GC(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating2, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123276856, "G", "C",
                "p.(Ser354Cys)", "c.1061C>G", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }

    @ParameterizedTest
    @CsvSource({
            "10-123276892-C-A, false, true", // https://www.ncbi.nlm.nih.gov/clinvar/variation/374819/
            "10-123276892-C-G, false, true", // https://www.ncbi.nlm.nih.gov/clinvar/variation/374820/
            "10-123276891-G-C, false, true", // https://www.ncbi.nlm.nih.gov/clinvar/variation/13275/
            "10-123276893-A-T, false, true", // https://www.ncbi.nlm.nih.gov/clinvar/variation/13267/
            "10-123276892-C-T, false, false" // dont match - same position proteinChange (same variant)
    })


    void testAssignPS1orPM5againstChr10Pos123276892CT(String variantStr, boolean expectedPs1, boolean expectedPm5) {
        AlleleProto.AlleleKey variant = parseAlleleKey(variantStr);
        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating2, variant);
        PS1PM5Assigner ps1PM5Assigner = new PS1PM5Assigner(variantDataService, jannovarAnnotator);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123276892, "C", "T",
                "p.(Cys342Tyr)", "c.1025G>A", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        ps1PM5Assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }
}
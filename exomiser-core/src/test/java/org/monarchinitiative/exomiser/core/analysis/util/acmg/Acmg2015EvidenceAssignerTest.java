/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2021 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.monarchinitiative.exomiser.core.analysis.util.acmg;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.checkerframework.checker.units.qual.A;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.monarchinitiative.exomiser.core.genome.*;

import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.*;
import org.monarchinitiative.exomiser.core.model.Pedigree.Individual;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.*;
import org.monarchinitiative.exomiser.core.phenotype.ModelPhenotypeMatch;
import org.monarchinitiative.exomiser.core.prioritisers.model.Disease;
import org.monarchinitiative.exomiser.core.prioritisers.model.InheritanceMode;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import java.util.List;
import java.util.Map;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.*;
import static org.monarchinitiative.exomiser.core.model.Pedigree.Individual.Sex.FEMALE;
import static org.monarchinitiative.exomiser.core.model.Pedigree.Individual.Sex.MALE;
import static org.monarchinitiative.exomiser.core.model.Pedigree.justProband;

class Acmg2015EvidenceAssignerTest {
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123276893, "A", "T",
                                                                    "p.(Cys342Ser)", "c.1024T>A", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                                                                    "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.START_LOST); //mocked

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123276856, "G", "C",
                "p.(Ser354Cys)", "c.1061C>G", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
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
        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123276892, "C", "T",
                "p.(Cys342Tyr)", "c.1025G>A", "FGFR2", VariantEffect.MISSENSE_VARIANT);

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.contains(AcmgCriterion.PS1), is(expectedPs1));
        assertThat(builder.contains(AcmgCriterion.PM5), is(expectedPm5));
    }

    @Test
    void testAssignsPVS1() {

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        // requires variant to be on a transcript predicted to undergo NMD in a LoF-intolerant gene for full PVS1
        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.START_LOST)
                .rankType(TranscriptAnnotation.RankType.EXON)
                .rank(1)
                .rankTotal(5)
                .build();
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                .geneSymbol("PTEN")
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                .variantEffect(VariantEffect.START_LOST)
                .annotations(List.of(transcriptAnnotation))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(PVS1).build()));
    }

    @ParameterizedTest
    @CsvSource({
            "MALE, AUTOSOMAL_DOMINANT, AUTOSOMAL_DOMINANT, true",
            "FEMALE, AUTOSOMAL_DOMINANT, AUTOSOMAL_DOMINANT, true",
            "UNKNOWN, AUTOSOMAL_DOMINANT, AUTOSOMAL_DOMINANT, true",

            "MALE, AUTOSOMAL_RECESSIVE, AUTOSOMAL_RECESSIVE, true",
            "FEMALE, AUTOSOMAL_RECESSIVE, AUTOSOMAL_RECESSIVE, true",
            "UNKNOWN, AUTOSOMAL_RECESSIVE, AUTOSOMAL_RECESSIVE, true",

            "MALE, AUTOSOMAL_DOMINANT, AUTOSOMAL_RECESSIVE, false",
            "FEMALE, AUTOSOMAL_DOMINANT, AUTOSOMAL_RECESSIVE, false",
            "UNKNOWN, AUTOSOMAL_DOMINANT, AUTOSOMAL_RECESSIVE, false",

            "MALE, X_RECESSIVE, X_RECESSIVE, true",
            "FEMALE, X_RECESSIVE, X_RECESSIVE, true",
            "UNKNOWN, X_RECESSIVE, X_RECESSIVE, false",

            "MALE, X_DOMINANT, X_DOMINANT, true",
            "FEMALE, X_DOMINANT, X_DOMINANT, true",
            "UNKNOWN, X_DOMINANT, X_DOMINANT, true",
    })
    void testAssignsPVS1(Individual.Sex probandSex, InheritanceMode diseaseInheritanceMode, ModeOfInheritance modeOfInheritance, boolean expectPvs1) {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", probandSex), jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.START_LOST)
                .rankType(TranscriptAnnotation.RankType.EXON)
                .rank(1)
                .rankTotal(5)
                .build();
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                .geneSymbol("PTEN")
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                .variantEffect(VariantEffect.START_LOST)
                .annotations(List.of(transcriptAnnotation))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(diseaseInheritanceMode).diseaseType(Disease.DiseaseType.DISEASE).build();
        List<Disease> knownDiseases = diseaseInheritanceMode.isCompatibleWith(modeOfInheritance) ? List.of(cowdenSyndrome) : List.of();
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, modeOfInheritance, List.of(variantEvaluation), knownDiseases, List.of());
        AcmgEvidence expected = expectPvs1 ? AcmgEvidence.builder().add(PVS1).build() : AcmgEvidence.empty();
        assertThat(acmgEvidence, equalTo(expected));
    }

    @Test
    void testAssignsPS2() {
        Individual proband = Individual.builder().id("proband").motherId("mother").fatherId("father").sex(MALE).status(Individual.Status.AFFECTED).build();
        Individual mother = Individual.builder().id("mother").sex(FEMALE).status(Individual.Status.UNAFFECTED).build();
        Individual father = Individual.builder().id("father").sex(MALE).status(Individual.Status.UNAFFECTED).build();
        Pedigree pedigree = Pedigree.of(proband, mother, father);
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", pedigree, jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. has frequency data - will not trigger PM2
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                // n.b. missense variant - will not trigger PVS1
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.het()),
                        SampleData.of("mother", SampleGenotype.homRef()),
                        SampleData.of("father", SampleGenotype.homRef())
                ))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        // n.b. low phenotype score - will not trigger PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), compatibleDiseaseMatches);

        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PS2).build()));
    }

    @Test
    void testAssignsPS2_hasFamilyHistory() {
        Individual proband = Individual.builder().id("proband").motherId("mother").fatherId("father").sex(MALE).status(Individual.Status.AFFECTED).build();
        Individual mother = Individual.builder().id("mother").sex(FEMALE).status(Individual.Status.UNAFFECTED).build();
        Individual father = Individual.builder().id("father").sex(MALE).status(Individual.Status.UNAFFECTED).build();
        Pedigree pedigree = Pedigree.of(proband, mother, father);
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", pedigree, jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. has frequency data - will not trigger PM2
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                // n.b. missense variant - will not trigger PVS1
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.het()),
                        SampleData.of("mother", SampleGenotype.het()), // Unaffected mother has same genotype - can't be PS2
                        SampleData.of("father", SampleGenotype.homRef())
                ))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        // n.b. low phenotype score - will not trigger PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), compatibleDiseaseMatches);

        assertThat(acmgEvidence, equalTo(AcmgEvidence.empty()));
    }

    @Test
    void testAssignsPM2() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty(), jannovarAnnotator, initializeVariantDataservice());
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 12345, "A", "G")
                // n.b. missing frequency data - will trigger PM2
                .frequencyData(FrequencyData.of())
                .build();
        // Requires variant to be in gene associated with a disorder in order that any ACMG criteria can be applied
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_RECESSIVE).diseaseType(Disease.DiseaseType.DISEASE).build();

        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());

        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PM2).build()));
    }

    @Test
    void testVariantMustBeInGeneWithKnownDiseaseAssociationForAcmgCriteriaToBeAssigned() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty(), jannovarAnnotator, initializeVariantDataservice());
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(1, 12345, "A", "G")
                // n.b. missing frequency data - should trigger PM2
                .frequencyData(FrequencyData.of())
                .build();
        // Requires variant to be in gene associated with a disorder in order that any ACMG criteria can be applied
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(), List.of());

        assertThat(acmgEvidence, equalTo(AcmgEvidence.empty()));
    }

    @Test
    void testAssignsPM3() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", null, jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89000000, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. has frequency data - will not trigger PM2
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                // n.b. missense variant - will not trigger PVS1
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.parseGenotype("0|1"))
                ))
                .build();

        VariantEvaluation pathogenic = TestFactory.variantBuilder(10, 89624227, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. start loss variant - will trigger PVS1
                .variantEffect(VariantEffect.START_LOST)
                .pathogenicityData(PathogenicityData.of(ClinVarData.builder().primaryInterpretation(ClinVarData.ClinSig.PATHOGENIC).reviewStatus("reviewed by expert panel").build()))
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.parseGenotype("1|0"))
                ))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_RECESSIVE).diseaseType(Disease.DiseaseType.DISEASE).build();
        // n.b. low phenotype score - will not trigger PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_RECESSIVE, List.of(variantEvaluation, pathogenic), List.of(cowdenSyndrome), compatibleDiseaseMatches);

        AcmgEvidence expected = AcmgEvidence.builder().add(AcmgCriterion.PM3).build();
        assertThat(acmgEvidence, equalTo(expected));
    }

    @Test
    void testAssignsBP2_InCisWithPathAR() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty(), jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89000000, "A", "G")
                // n.b. has frequency data - will not trigger PM2
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.parseGenotype("0|1"))
                ))
                .build();

        VariantEvaluation pathogenic = TestFactory.variantBuilder(10, 89624227, "A", "G")
                .pathogenicityData(PathogenicityData.of(ClinVarData.builder().primaryInterpretation(ClinVarData.ClinSig.PATHOGENIC).reviewStatus("reviewed by expert panel").build()))
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.parseGenotype("0|1"))
                ))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_RECESSIVE).diseaseType(Disease.DiseaseType.DISEASE).build();
        // n.b. low phenotype score - will not trigger PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_RECESSIVE, List.of(variantEvaluation, pathogenic), List.of(cowdenSyndrome), compatibleDiseaseMatches);

        AcmgEvidence expected = AcmgEvidence.builder().add(AcmgCriterion.BP2).build();
        assertThat(acmgEvidence, equalTo(expected));
    }

    @Test
    void testAssignsBP2_InTransWithPathAD() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty(), jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89000000, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. has frequency data - will not trigger PM2
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                // n.b. missense variant - will not trigger PVS1
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.parseGenotype("0|1"))
                ))
                .build();

        VariantEvaluation pathogenic = TestFactory.variantBuilder(10, 89624227, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. start loss variant - will trigger PVS1
                .variantEffect(VariantEffect.START_LOST)
                .pathogenicityData(PathogenicityData.of(ClinVarData.builder().primaryInterpretation(ClinVarData.ClinSig.PATHOGENIC).reviewStatus("reviewed by expert panel").build()))
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.parseGenotype("1|0"))
                ))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        // n.b. low phenotype score - will not trigger PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation, pathogenic), List.of(cowdenSyndrome), compatibleDiseaseMatches);

        AcmgEvidence expected = AcmgEvidence.builder().add(AcmgCriterion.BP2).build();
        assertThat(acmgEvidence, equalTo(expected));
    }

    @Test
    void testAssignsPM4() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                .geneSymbol("MUC6")
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                .variantEffect(VariantEffect.STOP_LOST)
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PM4).build()));
    }

    @Test
    void testAssignsPM4_NotAssignedPM4WhenPVS1Present() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());

        TranscriptAnnotation startLostAnnotation = TranscriptAnnotation.builder()
                .geneSymbol("PTEN")
                .accession("ENST00000371953.7")
                .variantEffect(VariantEffect.START_LOST)
                .rankType(TranscriptAnnotation.RankType.EXON)
                .rank(1)
                .rankTotal(9)
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                // haploinsufficient gene
                .geneSymbol("PTEN")
                .annotations(List.of(startLostAnnotation)) // prevent PM4 as PVS1 already triggered
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                .variantEffect(VariantEffect.STOP_LOST)
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PVS1).build()));
    }

    @Nested
    class ComputationalEvidence {

        @Test
        void testAssignsPP3() {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(
                            PathogenicityScore.of(PathogenicitySource.REVEL, 1.0f),
                            PathogenicityScore.of(PathogenicitySource.MVP, 1.0f)
                    ))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PP3, Evidence.STRONG).build()));
        }

        @ParameterizedTest
        @CsvSource ({
                "MVP, 1.0f, , ",
                "REVEL, 1.0f, PP3, STRONG"
        })
        void testAssignsPP3_singleScoreIsInsufficientUnlessItsRevel(PathogenicitySource pathogenicitySource, float pathogenicityScore, AcmgCriterion acmgCriterion, Evidence evidence) {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(PathogenicityScore.of(pathogenicitySource, pathogenicityScore)))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());

            AcmgEvidence expected = acmgCriterion == null ? AcmgEvidence.empty() : AcmgEvidence.builder().add(acmgCriterion, evidence).build();
            assertThat(acmgEvidence, equalTo(expected));
        }

        @Test
        void testAssignsPP3_majorityMustBePath() {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(
                            PathogenicityScore.of(PathogenicitySource.POLYPHEN, 1.0f),
                            PathogenicityScore.of(PathogenicitySource.MVP, 1.0f),
                            PathogenicityScore.of(PathogenicitySource.M_CAP, 0.0f)
                    ))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PP3).build()));
        }

        @Test
        void testPP3andPM4_majorityMustBePathOrBenign() {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(
                            PathogenicityScore.of(PathogenicitySource.POLYPHEN, 1.0f),
                            PathogenicityScore.of(PathogenicitySource.MVP, 1.0f),
                            PathogenicityScore.of(PathogenicitySource.M_CAP, 0.0f),
                            PathogenicityScore.of(PathogenicitySource.MUTATION_TASTER, 0.0f)
                    ))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
            assertThat(acmgEvidence, equalTo(AcmgEvidence.empty()));
        }

        @Test
        void testAssignsBP4() {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(
                            PathogenicityScore.of(PathogenicitySource.POLYPHEN, 0.0f),
                            PathogenicityScore.of(PathogenicitySource.MVP, 0.0f)
                    ))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.BP4).build()));
        }

        @ParameterizedTest
        @CsvSource ({
            "CADD, 0.0f, , ",
            "REVEL, 0.0f, BP4, VERY_STRONG"
        })
        void testAssignsBP4_singleScoreIsInsufficientIfNotRevel(PathogenicitySource pathogenicitySource, float pathogenicityScore, AcmgCriterion acmgCriterion, Evidence evidence) {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(PathogenicityScore.of(pathogenicitySource, pathogenicityScore)))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());

            AcmgEvidence expected = acmgCriterion == null ? AcmgEvidence.empty() : AcmgEvidence.builder().add(acmgCriterion, evidence).build();
            assertThat(acmgEvidence, equalTo(expected));
        }

        void testAssignsBP4_majorityMustBeBenign() {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(
                            PathogenicityScore.of(PathogenicitySource.POLYPHEN, 0.0f),
                            PathogenicityScore.of(PathogenicitySource.MVP, 0.0f),
                            PathogenicityScore.of(PathogenicitySource.M_CAP, 1.0f)
                    ))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.BP4).build()));
        }

        @ParameterizedTest
        @CsvSource ({
                "0.932, PP3, STRONG",
                "0.773, PP3, MODERATE",
                "0.644f, PP3, SUPPORTING",

                "0.290f, BP4, SUPPORTING",
                "0.183f, BP4, MODERATE",
                "0.016f, BP4, STRONG",
                "0.003f, BP4, VERY_STRONG",
        })
        public void testRevelOverridesAllOtherScores(float revelScore, AcmgCriterion acmgCriterion, Evidence evidence) {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                    .geneSymbol("PTEN")
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                    .pathogenicityData(PathogenicityData.of(
                            // REVEL, if present, should be used as the sole means of classifying the PP3/BP4
                            PathogenicityScore.of(PathogenicitySource.REVEL, revelScore),
                            PathogenicityScore.of(PathogenicitySource.MVP, 0.0f),
                            PathogenicityScore.of(PathogenicitySource.POLYPHEN, 1.0f),
                            PathogenicityScore.of(PathogenicitySource.SIFT, 0.0f),
                            PathogenicityScore.of(PathogenicitySource.M_CAP, 1.0f)
                    ))
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .build();
            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(acmgCriterion, evidence).build()));
        }
    }

    // PP4
    @Test
    void testAssignsPP4() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                .geneSymbol("PTEN")
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f))) // prevent PM2 assignment
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        // High phenotype match triggers - PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.6, cowdenSyndrome, List.of()));

        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), compatibleDiseaseMatches);
        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.PP4).build()));
    }

    @Nested
    class ClinicalEvidence {

        @ParameterizedTest
        @CsvSource(
                delimiter = ';',
                value = {
                        "criteria provided, single submitter; SUPPORTING",
                        "criteria provided, multiple submitters, no conflicts; STRONG",
                        "reviewed by expert panel; STRONG",
                        "practice guideline; STRONG",
                })
        void testAssignsPP5(String reviewStatus, AcmgCriterion.Evidence evidence) {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty(), jannovarAnnotator, initializeVariantDataservice());
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89000000, "A", "G")
                    // n.b. PTEN is a haploinsufficient gene
                    .geneSymbol("PTEN")
                    // n.b. has frequency data - will not trigger PM2
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                    // n.b. missense variant - will not trigger PVS1
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .pathogenicityData(PathogenicityData.of(ClinVarData.builder().primaryInterpretation(ClinVarData.ClinSig.PATHOGENIC).reviewStatus(reviewStatus).build()))
                    .build();

            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            // n.b. low phenotype score - will not trigger PP4
            List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), compatibleDiseaseMatches);

            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(PP5, evidence).build()));
        }

        @ParameterizedTest
        @CsvSource(
                delimiter = ';',
                value = {
                        "criteria provided, single submitter; SUPPORTING",
                        "criteria provided, multiple submitters, no conflicts; STRONG",
                        "reviewed by expert panel; STRONG",
                        "practice guideline; STRONG",
                })
        void testAssignsBP6(String reviewStatus, AcmgCriterion.Evidence evidence) {
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty(), jannovarAnnotator, initializeVariantDataservice());
            // https://www.ncbi.nlm.nih.gov/clinvar/variation/127667/
            VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89622915, "A", "G")
                    // n.b. PTEN is a haploinsufficient gene
                    .geneSymbol("PTEN")
                    // n.b. has frequency data - will not trigger PM2
                    .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AFRICAN_INC_AFRICAN_AMERICAN, 1.42f)))
                    // n.b. missense variant - will not trigger PVS1
                    .variantEffect(VariantEffect.MISSENSE_VARIANT)
                    .pathogenicityData(PathogenicityData.of(ClinVarData.builder().primaryInterpretation(ClinVarData.ClinSig.BENIGN).reviewStatus(reviewStatus).build()))
                    .build();

            Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
            // n.b. low phenotype score - will not trigger PP4
            List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
            AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), compatibleDiseaseMatches);

            assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(BP6, evidence).build()));
        }
    }

    @Test
    void testAssignsBA1() {
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, initializeVariantDataservice());
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                .geneSymbol("PTEN")
                // high allele freq - triggers BA1 assignment
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 5.0f)))
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();

        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), List.of());
        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(AcmgCriterion.BA1).build()));
    }

    @Test
    void testAssignsBS4() {
        Individual proband = Individual.builder().id("proband").motherId("mother").fatherId("father").sex(MALE).status(Individual.Status.AFFECTED).build();
        Individual mother = Individual.builder().id("mother").sex(FEMALE).status(Individual.Status.AFFECTED).build();
        Individual father = Individual.builder().id("father").sex(MALE).status(Individual.Status.UNAFFECTED).build();
        Pedigree pedigree = Pedigree.of(proband, mother, father);
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", pedigree, jannovarAnnotator, initializeVariantDataservice());
        // https://www.ncbi.nlm.nih.gov/clinvar/variation/484600/ 3* PATHOGENIC variant  - reviewed by expert panel
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 89624227, "A", "G")
                // n.b. PTEN is a haploinsufficient gene
                .geneSymbol("PTEN")
                // n.b. has frequency data - will not trigger PM2
                .frequencyData(FrequencyData.of(Frequency.of(FrequencySource.EXAC_AMERICAN, 0.1f)))
                // n.b. missense variant - will not trigger PVS1
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .sampleGenotypes(SampleGenotypes.of(
                        SampleData.of("proband", SampleGenotype.het()),
                        SampleData.of("mother", SampleGenotype.homRef()), // Affected mother has different genotype - can't be PS2
                        SampleData.of("father", SampleGenotype.homRef())
                ))
                .build();
        Disease cowdenSyndrome = Disease.builder().diseaseId("OMIM:158350").diseaseName("COWDEN SYNDROME 1; CWS1").inheritanceMode(InheritanceMode.AUTOSOMAL_DOMINANT).diseaseType(Disease.DiseaseType.DISEASE).build();
        // n.b. low phenotype score - will not trigger PP4
        List<ModelPhenotypeMatch<Disease>> compatibleDiseaseMatches = List.of(ModelPhenotypeMatch.of(0.5, cowdenSyndrome, List.of()));
        AcmgEvidence acmgEvidence = instance.assignVariantAcmgEvidence(variantEvaluation, ModeOfInheritance.AUTOSOMAL_DOMINANT, List.of(variantEvaluation), List.of(cowdenSyndrome), compatibleDiseaseMatches);

        assertThat(acmgEvidence, equalTo(AcmgEvidence.builder().add(BS4).build()));
    }

//    @Test
//    void testRatio() {
////        VariantDataService variantDataService = initializeCustomVariantDataservice(clinVarPathogenicStarRating2, variant);
//        ClinVarGeneStats clinVarGeneStats = new ClinVarGeneStats("Test", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 12, ClinVarData.ClinSig.BENIGN, 2)));
////        MVStore mvStore = mvStore();
//        VariantDataService variantDataService = variantDataServiceForGeneStats("Test", clinVarGeneStats);
//
////        MVMap<String, ClinVarGeneStats> stringClinVarGeneStatsMVMap = MvStoreUtil.openGeneStatsMVMap(mvStore);
////        stringClinVarGeneStatsMVMap.put("Test", clinVarGeneStats);
////        variantDataService.getClinVarGeneStats("Test");
//        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);
//        VariantEvaluation variantEvaluation = buildVariantEvaluation(10, 123247514, "C", "A",
//                "p.(Lys659Asn)", "c.1977G>T", "FGFR2", VariantEffect.MISSENSE_VARIANT);
//
//        AcmgEvidence.Builder builder = AcmgEvidence.builder();
//
//        assigner.assignPP2(builder, variantEvaluation);
//        assertThat(builder.contains(AcmgCriterion.PP2), is(true));
//
//    }
}
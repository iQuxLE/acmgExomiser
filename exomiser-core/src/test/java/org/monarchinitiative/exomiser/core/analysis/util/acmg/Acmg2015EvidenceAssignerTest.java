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
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.monarchinitiative.exomiser.core.genome.*;

import org.monarchinitiative.exomiser.core.model.*;
import org.monarchinitiative.exomiser.core.model.Pedigree.Individual;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.monarchinitiative.exomiser.core.phenotype.ModelPhenotypeMatch;
import org.monarchinitiative.exomiser.core.prioritisers.model.Disease;
import org.monarchinitiative.exomiser.core.prioritisers.model.InheritanceMode;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.*;
import static org.monarchinitiative.exomiser.core.model.Pedigree.Individual.Sex.FEMALE;
import static org.monarchinitiative.exomiser.core.model.Pedigree.Individual.Sex.MALE;
import static org.monarchinitiative.exomiser.core.model.Pedigree.justProband;

class Acmg2015EvidenceAssignerTest {
    // https://www.ncbi.nlm.nih.gov/clinvar/variation/1698211/
    private AlleleProto.AlleleKey variant1aChr10Pos123247514 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123247514)
            .setRef("T")
            .setAlt("G")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/661397/
    private AlleleProto.AlleleKey variant1bChr10Pos123247514 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123247514)
            .setRef("C")
            .setAlt("G")
            .build();

    // mocked Variant
    private AlleleProto.AlleleKey variant2Chr10Pos123247515 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123247515)
            .setRef("T")
            .setAlt("C")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/374820/
    private AlleleProto.AlleleKey variant3Chr10Pos123276892 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276892)
            .setRef("C")
            .setAlt("G")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/13267/
    private AlleleProto.AlleleKey variant4Chr10Pos123276893 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276893)
            .setRef("A")
            .setAlt("T")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/689498/
    private AlleleProto.AlleleKey variant5Chr11Pos108124619 = AlleleProto.AlleleKey.newBuilder()
            .setChr(11)
            .setPosition(299372)
            .setRef("G")
            .setAlt("C")
            .build();

    // mockVariants

    private AlleleProto.AlleleKey variant6Chr11Pos123276893 = AlleleProto.AlleleKey.newBuilder()
            .setChr(11)
            .setPosition(123276893)
            .setRef("A") //T
            .setAlt("T") //A
            .build();


    //    https://www.ncbi.nlm.nih.gov/clinvar/variation/1317030/
    private final AlleleProto.AlleleKey variant7Chr10Pos123245027 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123245027)
            .setRef("T")
            .setAlt("C")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/13264/
    private final AlleleProto.AlleleKey variant8Chr10Pos123276899 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276899)
            .setRef("A")
            .setAlt("G")
            .build();

    //https://www.ncbi.nlm.nih.gov/clinvar/variation/1066784/
    private final AlleleProto.AlleleKey variant9Chr10Pos123276856 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276856)
            .setRef("G")
            .setAlt("A")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/2065795/
    private final AlleleProto.AlleleKey variant11Chr10Pos123245026 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123245026)
            .setRef("A")
            .setAlt("G")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/13265/
    private final AlleleProto.AlleleKey variant10Chr10Pos123276856 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276856)
            .setRef("G")
            .setAlt("C")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/374819/
    private final AlleleProto.AlleleKey variant12Chr10Pos123276892 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276892)
            .setRef("C")
            .setAlt("A")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/374820/
    private final AlleleProto.AlleleKey variant13Chr10Pos123276892 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276892)
            .setRef("C")
            .setAlt("G")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/13275/
    private final AlleleProto.AlleleKey variant14Chr10Pos123276893 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276891)
            .setRef("G")
            .setAlt("C")
            .build();

    // https://www.ncbi.nlm.nih.gov/clinvar/variation/13267/
    private final AlleleProto.AlleleKey variant15Chr10Pos123276891 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(123276893)
            .setRef("A")
            .setAlt("T")
            .build();

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

    private final AlleleProto.ClinVar clinVarUncertainSignificanceStarRating2 = AlleleProto.ClinVar.newBuilder()
            .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.UNCERTAIN_SIGNIFICANCE)
            .setReviewStatus("criteria provided, multiple submitters, no conflicts")
            .build();


    private final MVStore mvStore = new MVStore.Builder().compress().open();


    @Test
    void throwsExceptionWithMismatchedIds() {
        assertThrows(IllegalArgumentException.class, () -> new Acmg2015EvidenceAssigner("Zaphod", justProband("Ford", MALE), null, null));
    }

    private final JannovarVariantAnnotator jannovarAnnotator = new JannovarVariantAnnotator(TestFactory.getDefaultGenomeAssembly(), TestFactory
            .buildDefaultJannovarData(), ChromosomalRegionIndex.empty());

    @Test
    void testAssignPS1_SamePositionSameNucleotideSameProteinChangeDifferentCdna() {
        TestVariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //match
                .put(variant3Chr10Pos123276892, clinVarPathogenicStarRating2)
                // PS1 no match cause same cdna, PM5 no match cause same proteinChange
                .put(variant4Chr10Pos123276893, clinVarPathogenicStarRating2)
                // no match cause wrong codon
                .put(variant1bChr10Pos123247514, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Cys342Ser)")
                .hgvsCdna("c.1024T>A")
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123276893, "A", "T")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), not(equalTo(AcmgEvidence.empty())));
        assertThat(builder.contains(AcmgCriterion.PS1), is(true));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(instance.getProcessedVariantCount(), is(1));
    }

    @Test
    void testAssignPS1_mockedInputcDNAchange_TwoHitChr10Pos123276893() {
        TestVariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //match
                .put(variant3Chr10Pos123276892, clinVarPathogenicStarRating2)
                // PS1 match cause of mocked input cDNA change, PM5 no match cause proteinChange equals
                .put(variant4Chr10Pos123276893, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Cys342Ser)")
                .hgvsCdna("c.1024T>T") // this is mocked: there is no such variant but we test to hit 2 variants
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123276893, "A", "T")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), not(equalTo(AcmgEvidence.empty())));
        assertThat(builder.contains(AcmgCriterion.PS1), is(true));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(instance.getProcessedVariantCount(), is(2));
    }
    @Test
    void testAssignPS1andPM5_TwoHitsassignPM1andOneHitPM5Chr10Pos123247514() {
        TestVariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //match for PS1 cause different cDNA and same proteinChange
                .put(variant1aChr10Pos123247514, clinVarPathogenicStarRating2)
                // match for PS1 cause different cDNA and same proteinChange
                .put(variant1bChr10Pos123247514, clinVarPathogenicStarRating2)
                // match for PM5 (mocked variant) cause different/new proteinChange - does not hit PS1 cause of same reason
                .put(variant2Chr10Pos123247515, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Lys659Asn)")
                .hgvsCdna("c.1977G>T")
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123247514, "C", "A")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), not(equalTo(AcmgEvidence.empty())));
        assertThat(builder.contains(AcmgCriterion.PS1), is(true));
        assertThat(builder.contains(AcmgCriterion.PM5), is(true));
        assertThat(instance.getProcessedVariantCount(), is(3));
    }

    @Test
    void testAssignPS1orPM5_completeMismatchNoHitDifferentChromosomeAndDifferentCodon() {
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //no match different Codon
                .put(variant1bChr10Pos123247514, clinVarPathogenicStarRating2)
                // no match wrong Chromosome
                .put(variant5Chr11Pos108124619, clinVarPathogenicStarRating2)
                // no match wrong Chromosome but right position
                .put(variant6Chr11Pos123276893, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator,variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Cys342Ser)")
                .hgvsCdna("c.1024T>A")
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123276893, "A", "T")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(instance.getProcessedVariantCount(), is(0));
    }

    @Test
    void testAssignPS1_NoMissense() {

        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                .put(variant1bChr10Pos123247514, clinVarPathogenicStarRating2)
                .put(variant2Chr10Pos123247515, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.START_LOST) // this is a fake START_LOST, the real on this position is MISSENSE
                .hgvsProtein("p.(Lys659Asn)")
                .hgvsCdna("c.1977G>T")
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123247514, "C", "A")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.START_LOST) // this is a fake START_LOST, the real on this position is MISSENSE
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(instance.getProcessedVariantCount(), is(0));
    }

    @Test
    void testAssignPS1_mocked_NotPathogenicSameCodon() {
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //no match - mocked
                .put(variant2Chr10Pos123247515, clinVarBenignStarRating2)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Lys659Asn)")
                .hgvsCdna("c.1977G>T")
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123247514, "C", "A")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();


        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(instance.getProcessedVariantCount(), is(0));
    }

    @Test
    void testAssignPS1_mocked_starRating1() {
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //no match - mocked
                .put(variant2Chr10Pos123247515, clinVarPathogenicStarRating1)
                .build();

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Lys659Asn)")
                .hgvsCdna("c.1977G>T")
                .build();

        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123247514, "C", "A")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();


        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        instance.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(instance.getProcessedVariantCount(), is(0));
    }

    @Test
    void testAssignPM5_SamePositionSameNucleotidDifferentProteinChangeChr10Pos123276856(){
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                .put(variant9Chr10Pos123276856, clinVarPathogenicStarRating2)
                //no match - cause same protein change
                .put(variant10Chr10Pos123276856, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Ser354Cys)")
                .hgvsCdna("c.1061C>G")
                .build();

        // NM_000141.5(FGFR2):c.1061C>G (p.Ser354Cys)  https://www.ncbi.nlm.nih.gov/clinvar/variation/13265/
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123276856, "G", "C")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), not(equalTo(AcmgEvidence.empty())));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(true));
        assertThat(assigner.getProcessedVariantCount(), is(1));
    }

    @Test
    void testAssignPM5_SamePositionSameNucleotidDifferentProteinChangeChr10Pos123276892(){
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                // all match ( remember when going over this if we only need the codon so we actually need a real codon check for stuff like that ? ) or does it not matter cause swapping posis after changes?
                .put(variant12Chr10Pos123276892, clinVarPathogenicStarRating2)
                .put(variant13Chr10Pos123276892, clinVarPathogenicStarRating2)
                .put(variant14Chr10Pos123276893, clinVarPathogenicStarRating2)
                .put(variant15Chr10Pos123276891, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Cys342Tyr)")
                .hgvsCdna("c.1025G>A")
                .build();

        // NM_000141.5(FGFR2):c.1025G>A (p.Cys342Tyr) https://www.ncbi.nlm.nih.gov/clinvar/variation/13263/
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123276892, "C", "T")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), not(equalTo(AcmgEvidence.empty())));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(true));
        assertThat(assigner.getProcessedVariantCount(), is(4));
    }

    @Test
    void testAssignPM5_dontFindSamePositionSameProteinChange(){
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //no match
                .put(variant11Chr10Pos123245026, clinVarPathogenicStarRating2)
                .build();

        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Met693Thr)")
                .hgvsCdna("c.2078T>C")
                .build();

        // NM_000141.5(FGFR2):c.2078T>C (p.Met693Thr)  https://www.ncbi.nlm.nih.gov/clinvar/variation/2065795/
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123245026, "A", "G")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(assigner.getProcessedVariantCount(), is(0));
    }

    @Test
    void testAssignPM5_NotPathogenicButSameCodon(){
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                //no match - cause not pathogenic
                .put(variant7Chr10Pos123245027, clinVarUncertainSignificanceStarRating2)
                .build();

        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Met693Thr)")
                .hgvsCdna("c.2078T>C")
                .build();

        // NM_000141.5(FGFR2):c.2078T>C (p.Met693Thr)  https://www.ncbi.nlm.nih.gov/clinvar/variation/2065795/
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123245026, "T", "C")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(assigner.getProcessedVariantCount(), is(0));
    }

    @Test
    void testAssignPM5_LowStarRatingButSameCodon(){
        VariantDataService variantDataService = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                // no match - cause of starRating (mocked clinVarData)
                .put(variant7Chr10Pos123245027, clinVarPathogenicStarRating1)
                .build();

        Acmg2015EvidenceAssigner assigner = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE), jannovarAnnotator, variantDataService);

        TranscriptAnnotation transcriptAnnotation = TranscriptAnnotation.builder()
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .hgvsProtein("p.(Met693Thr)")
                .hgvsCdna("c.2078T>C")
                .build();

        // NM_000141.5(FGFR2):c.2078T>C (p.Met693Thr) https://www.ncbi.nlm.nih.gov/clinvar/variation/2065795/
        VariantEvaluation variantEvaluation = TestFactory.variantBuilder(10, 123245026, "A", "G")
                .geneSymbol("FGFR2")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .annotations(List.of(transcriptAnnotation))
                .build();

        AcmgEvidence.Builder builder = AcmgEvidence.builder();
        assigner.assignPS1orPM5(builder, variantEvaluation);
        assertThat(builder.build(), equalTo(AcmgEvidence.empty()));
        assertThat(builder.contains(AcmgCriterion.PS1), is(false));
        assertThat(builder.contains(AcmgCriterion.PM5), is(false));
        assertThat(assigner.getProcessedVariantCount(), is(0));
    }



    @Test
    void testAssignsPVS1() {

        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", probandSex));
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", pedigree);
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", pedigree);
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty());
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty());
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", null);
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty());
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty());
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));

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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty());
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
            Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", Pedigree.empty());
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", justProband("proband", MALE));
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
        Acmg2015EvidenceAssigner instance = new Acmg2015EvidenceAssigner("proband", pedigree);
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
}
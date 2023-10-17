package org.monarchinitiative.exomiser.core.genome.dao;

import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.Strand;

import java.nio.file.Path;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;

class ClinVarDaoMvStoreTest {

    private final  AlleleProto.AlleleKey positionStartMinus1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1229)
            .setAlt("A")
            .setRef("G")
            .build();

    private final  AlleleProto.AlleleKey positionEndPlus1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1231)
            .setAlt("C")
            .setRef("G")
            .build();

    // Variants directly at the boundaries:
    private final  AlleleProto.AlleleKey positionStartMinus2 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1228)
            .setAlt("A")
            .setRef("G")
            .build();

    private final  AlleleProto.AlleleKey positionEndPlus2 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1232)
            .setAlt("T")
            .setRef("A")
            .build();


    // Variants just outside the boundaries:
    private final  AlleleProto.AlleleKey positionStartMinus3 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1227)
            .setAlt("A")
            .setRef("G")
            .build();

    private final AlleleProto.AlleleKey positionEndPlus3 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1233)
            .setAlt("C")
            .setRef("T")
            .build();

    // Variants are exactly matching on position:
    private final  AlleleProto.AlleleKey positionExactlyMatches = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1230)
            .setAlt("T")
            .setRef("G")
            .build();

    // Variants are exactly matching on position + but ref and alt are different
    private final  AlleleProto.AlleleKey positionExactlyMatchesButDifferentRefAlt = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1230)
            .setAlt("A")
            .setRef("G")
            .build();

    // Variants are far out of boundaries:
    private final  AlleleProto.AlleleKey positionFarOutOfBoundariesMinus = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(7700)
            .setAlt("A")
            .setRef("G")
            .build();

    private final  AlleleProto.AlleleKey positionFarOutOfBoundariesPlus = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1)
            .setAlt("A")
            .setRef("G")
            .build();

    // two variants that are in range but Alt.length and Ref.length > 1

    private final  AlleleProto.AlleleKey positionInButLengthIncorrect = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1229)
            .setAlt("AA")
            .setRef("G")
            .build();

    private final  AlleleProto.AlleleKey positionInButLengthIncorrect1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1231)
            .setAlt("A")
            .setRef("GG")
            .build();




    private MVStore newMvStore() {
        // open the store (in-memory if fileName is null)
        return new MVStore.Builder()
                .compress()
                .open();
    }

    private GenomeAssembly genomeAssembly;
    private MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;
    private MVStore mvStore;

    @BeforeEach
    public void setUp(){
        genomeAssembly = GenomeAssembly.HG19;
        mvStore = newMvStore();
        clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);

    }

    @Test
    public void FiveInsideAndFourOutsideBoundaries(@TempDir Path tempDir){
        clinVarMap.put(positionStartMinus1, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        clinVarMap.put(positionEndPlus1, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.LIKELY_PATHOGENIC)
                .build());

        clinVarMap.put(positionStartMinus2, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.BENIGN)
                .build());

        clinVarMap.put(positionEndPlus2, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.BENIGN_OR_LIKELY_BENIGN)
                .build());

        clinVarMap.put(positionStartMinus3, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.CONFLICTING_PATHOGENICITY_INTERPRETATIONS)
                .build());

        clinVarMap.put(positionEndPlus3, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        clinVarMap.put(positionExactlyMatches, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        clinVarMap.put(positionExactlyMatchesButDifferentRefAlt, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        clinVarMap.put(positionFarOutOfBoundariesMinus, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        clinVarMap.put(positionFarOutOfBoundariesPlus, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2, 2));
        assertThat(result.isEmpty(), is(false));
        assertThat(result.size(), is(6));
    }

    @Test
    public void variantsNonOverlap(@TempDir Path tempDir){

        clinVarMap.put(positionStartMinus3, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        clinVarMap.put(positionEndPlus3, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        GenomeAssembly genomeAssembly = GenomeAssembly.HG19;
        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        ClinVarDaoMvStore clinVarDaoMvStore = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDaoMvStore.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(result.isEmpty(), is(true));
        assertThat(result.size(), is(0));

    }

    @Test
    public void variantsOverlapDirectlyAtBoundaries(@TempDir Path tempDir){
        MVStore mvStore = newMvStore();
        MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);

        clinVarMap.put(positionStartMinus2, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        clinVarMap.put(positionEndPlus2, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        GenomeAssembly genomeAssembly = GenomeAssembly.HG19;
        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(result.isEmpty(), is(false));
        assertThat(result.size(), is(2));
    }

    @Test
    public void variantsOverlapInsideBoundaries(@TempDir Path tempDir){


        clinVarMap.put(positionStartMinus1, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        clinVarMap.put(positionEndPlus1, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.NOT_PROVIDED)
                .build());

        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(result.isEmpty(), is(false));
        assertThat(result.size(), is(2));
    }

    @Test
    public void variantsExactlyMatchOnPosition(){

        clinVarMap.put(positionExactlyMatches, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation);
        assertThat(result.isEmpty(), is(false));
        assertThat(result.size(), is(1));
    }

    @Test
    public void variantsAreFarOutOfBoundaries(@TempDir Path tempDir){

        clinVarMap.put(positionFarOutOfBoundariesMinus, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        clinVarMap.put(positionFarOutOfBoundariesPlus, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230,1230), "T", "A")
                .build();

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(result.isEmpty(), is(true));
        assertThat(result.size(), is(0));

    }

    @Test
    public void noVariantsGettingAddedNoOverlap(){
        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(result.isEmpty(), is(true));
        assertThat(result.size(), is(0));

    }
    @Test
    public void setPositionInButLengthIncorrectTest(){
        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(genomeAssembly.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1230, 1230),"T", "A" )
                .build();

        clinVarMap.put(positionInButLengthIncorrect, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        clinVarMap.put(positionInButLengthIncorrect1, AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build());

        ClinVarDaoMvStore clinVarDao = new ClinVarDaoMvStore(mvStore, genomeAssembly);

        var result  = clinVarDao.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(result.isEmpty(), is(false));
        assertThat(result.size(), is(2));

    }

    @Test
    void getClinVarData() {
        try (MVStore mvStore = new MVStore.Builder().open()) {
            var clinvarMap = MvStoreUtil.openClinVarMVMap(mvStore);
            AlleleProto.AlleleKey alleleKey = AlleleProto.AlleleKey.newBuilder().setChr(1).setPosition(200).setRef("A").setAlt("T").build();
            AlleleProto.ClinVar clinVar = AlleleProto.ClinVar.newBuilder()
                    .setAlleleId("12345")
                    .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                    .setReviewStatus("criteria_provided,_multiple_submitters,_no_conflicts")
                    .build();
            clinvarMap.put(alleleKey, clinVar);
            GenomeAssembly genomeAssembly = GenomeAssembly.HG19;
            ClinVarDao instance = new ClinVarDaoMvStore(mvStore, genomeAssembly);
            Variant clinVarVariant = VariantEvaluation.builder()
                    .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(200, 200), "A", "T")
                    .build();

            assertThat(instance.getClinVarData(clinVarVariant), equalTo(AlleleProtoAdaptor.toClinVarData(clinVar)));

            Variant nonClinVarVariant = VariantEvaluation.builder()
                    .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(200, 200), "A", "A")
                    .build();

            assertThat(instance.getClinVarData(nonClinVarVariant), equalTo(ClinVarData.empty()));
        }
    }
}
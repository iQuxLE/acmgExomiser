package org.monarchinitiative.exomiser.core.genome;

import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.svart.*;

import java.util.Map;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.jupiter.api.Assertions.*;

class TestVariantDataServiceTest {

    private final AlleleProto.AlleleKey positionStartMinus1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(1229)
            .setAlt("A")
            .setRef("G")
            .build();

    private final AlleleProto.AlleleKey positionEndPlus1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(1231)
            .setAlt("C")
            .setRef("G")
            .build();

    private final AlleleProto.AlleleKey positionPlus2 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(1232)
            .setAlt("T")
            .setRef("G")
            .build();

    private final AlleleProto.AlleleKey positionMinus2 = AlleleProto.AlleleKey.newBuilder()
            .setChr(10)
            .setPosition(1228)
            .setAlt("C")
            .setRef("A")
            .build();


    private MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;
    private MVStore mvStore;


    private MVStore newMvStore() {
        return new MVStore.Builder().compress().open();
    }

//    private static ClinVarData getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(AlleleProto.ClinVar clinVar){
//        ClinVarData.ClinSig primaryInterpretation = AlleleProtoAdaptor.fromClinVarPrototoClinSig1(clinVar.getPrimaryInterpretation());
//        return ClinVarData.builder().primaryInterpretation(primaryInterpretation).build();
//
//    }

    @BeforeEach
    public void setUp() {
        mvStore = newMvStore();
        clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
    }

    private static GenomicVariant getGenomicVariantFromProtoAlleleKey(AlleleProto.AlleleKey alleleKey, GenomeAssembly genomeAssembly){
        return GenomicVariant.builder().variant(genomeAssembly.getContigById(alleleKey.getChr()), Strand.POSITIVE,
                Coordinates.oneBased(alleleKey.getPosition(),alleleKey.getPosition()), alleleKey.getRef(), alleleKey.getAlt()).build();
    }

    private static ClinVarData getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(AlleleProto.ClinVar clinVar){
        ClinVarData.ClinSig primaryInterpretation = toClinSig(clinVar.getPrimaryInterpretation());
        return ClinVarData.builder().primaryInterpretation(primaryInterpretation).build();

    }

    private static ClinVarData.ClinSig toClinSig(AlleleProto.ClinVar.ClinSig protoClinSig) {
        switch (protoClinSig) {
            case BENIGN:
                return ClinVarData.ClinSig.BENIGN;
            case BENIGN_OR_LIKELY_BENIGN:
                return ClinVarData.ClinSig.BENIGN_OR_LIKELY_BENIGN;
            case LIKELY_BENIGN:
                return ClinVarData.ClinSig.LIKELY_BENIGN;
            case UNCERTAIN_SIGNIFICANCE:
                return ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE;
            case LIKELY_PATHOGENIC:
                return ClinVarData.ClinSig.LIKELY_PATHOGENIC;
            case PATHOGENIC_OR_LIKELY_PATHOGENIC:
                return ClinVarData.ClinSig.PATHOGENIC_OR_LIKELY_PATHOGENIC;
            case PATHOGENIC:
                return ClinVarData.ClinSig.PATHOGENIC;
            case CONFLICTING_PATHOGENICITY_INTERPRETATIONS:
                return ClinVarData.ClinSig.CONFLICTING_PATHOGENICITY_INTERPRETATIONS;
            case AFFECTS:
                return ClinVarData.ClinSig.AFFECTS;
            case ASSOCIATION:
                return ClinVarData.ClinSig.ASSOCIATION;
            case DRUG_RESPONSE:
                return ClinVarData.ClinSig.DRUG_RESPONSE;
            case OTHER:
                return ClinVarData.ClinSig.OTHER;
            case PROTECTIVE:
                return ClinVarData.ClinSig.PROTECTIVE;
            case RISK_FACTOR:
                return ClinVarData.ClinSig.RISK_FACTOR;
            case NOT_PROVIDED:
            case UNRECOGNIZED:
            default:
                return ClinVarData.ClinSig.NOT_PROVIDED;
        }
    }

    @Test
    public void testClinVarDaoMvStoreViaVariantDataService() {
        AlleleProto.ClinVar clinVar1 = AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build();

        AlleleProto.ClinVar clinVar2 = AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.LIKELY_PATHOGENIC)
                .build();

        AlleleProto.ClinVar clinVar3 = AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.BENIGN)
                .build();

        AlleleProto.ClinVar clinVar4 = AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.LIKELY_BENIGN)
                .build();

        TestVariantDataService service = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                .put(positionEndPlus1, clinVar1)
                .put(positionStartMinus1, clinVar2)
                .put(positionPlus2, clinVar3)
                .put(positionMinus2, clinVar4)
                .build();

        //VariantEvaluation extends GenomicIntervall
        VariantEvaluation variantEvaluation = VariantEvaluation.builder().variant(GenomeAssembly.HG19.getContigById(10), Strand.POSITIVE, Coordinates.oneBased(1230, 1230), "T", "G").build();

        Map<GenomicVariant, ClinVarData> resultMap = service.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2,2));
        assertThat(resultMap.size(), is(4));

        // transforming keys to match datatypes from resultMap
        GenomicVariant keyVariant1 = getGenomicVariantFromProtoAlleleKey(positionEndPlus1, GenomeAssembly.HG19);
        GenomicVariant keyVariant2 = getGenomicVariantFromProtoAlleleKey(positionStartMinus1, GenomeAssembly.HG19);
        GenomicVariant keyVariant3 = getGenomicVariantFromProtoAlleleKey(positionPlus2, GenomeAssembly.HG19);
        GenomicVariant keyVariant4 = getGenomicVariantFromProtoAlleleKey(positionMinus2, GenomeAssembly.HG19);
        ClinVarData clinVarValue1 = getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(clinVar1);
        ClinVarData clinVarValue2 = getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(clinVar2);
        ClinVarData clinVarValue3 = getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(clinVar3);
        ClinVarData clinVarValue4 = getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(clinVar4);

        assertTrue(resultMap.containsKey(keyVariant1));
        assertTrue(resultMap.containsKey(keyVariant2));
        assertTrue(resultMap.containsKey(keyVariant3));
        assertTrue(resultMap.containsKey(keyVariant4));
        assertEquals(clinVarValue1, resultMap.get(keyVariant1));
        assertEquals(clinVarValue2, resultMap.get(keyVariant2));
        assertEquals(clinVarValue3, resultMap.get(keyVariant3));
        assertEquals(clinVarValue4, resultMap.get(keyVariant4));
    }
}


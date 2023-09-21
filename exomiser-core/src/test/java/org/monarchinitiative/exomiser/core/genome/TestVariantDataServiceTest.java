package org.monarchinitiative.exomiser.core.genome;

import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
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


    private MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;
    private MVStore mvStore;


    private MVStore newMvStore() {
        return new MVStore.Builder().compress().open();
    }

    @BeforeEach
    public void setUp() {
        mvStore = newMvStore();
        clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
    }

    @Test
    public void testClinVarDaoMvStoreViaVariantDataService() {
        AlleleProto.ClinVar clinVar1 = AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.PATHOGENIC)
                .build();

        AlleleProto.ClinVar clinVar2 = AlleleProto.ClinVar.newBuilder()
                .setPrimaryInterpretation(AlleleProto.ClinVar.ClinSig.LIKELY_PATHOGENIC)
                .build();

        TestVariantDataService service = TestVariantDataService.builder()
                .setMVStore(mvStore)
                .setGenomeAssembly(GenomeAssembly.HG19)
                .putAK(positionEndPlus1, clinVar1)
                .putAK(positionStartMinus1, clinVar2)
                .build();

        GenomicInterval genomicInterval = GenomicInterval.of(GenomeAssembly.HG19.getContigById(10), Strand.POSITIVE, Coordinates.oneBased(1230, 1230));

        Map<GenomicVariant, ClinVarData> resultMap = service.findClinVarDataOverlappingGenomicInterval(genomicInterval);
        assertThat(resultMap.size(), is(2));

        // transforming keys to match datatypes from resultMap
        GenomicVariant keyVariant1 = AlleleProtoAdaptor.getGenomicVariantFromProtoAlleleKey(positionEndPlus1, GenomeAssembly.HG19);
        GenomicVariant keyVariant2 = AlleleProtoAdaptor.getGenomicVariantFromProtoAlleleKey(positionStartMinus1, GenomeAssembly.HG19);
        ClinVarData clinVarValue1 = AlleleProtoAdaptor.getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(clinVar1);
        ClinVarData clinVarValue2 = AlleleProtoAdaptor.getClinVarDataPrimaryInterpretationFromAlleleProtoClinVar(clinVar2);

        assertTrue(resultMap.containsKey(keyVariant1));
        assertTrue(resultMap.containsKey(keyVariant2));
        assertEquals(clinVarValue1, resultMap.get(keyVariant1));
        assertEquals(clinVarValue2, resultMap.get(keyVariant2));
    }
}


package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import com.google.common.collect.ImmutableMap;
import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.TestFactory;
import org.monarchinitiative.exomiser.core.genome.dao.*;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.GenomicVariant;
import org.monarchinitiative.svart.Strand;

import java.util.Map;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.jupiter.api.Assertions.*;

class BS1BS2AssignerTest {

    private Variant buildVariant(int chr, int pos, String ref, String alt) {
        return TestFactory.variantBuilder(chr, pos, ref, alt).build();
    }

    AllelePropertiesDaoAdapter newInstanceWithData(Map<AlleleProto.AlleleKey, AlleleProto.AlleleProperties> value) {
        MVStore mvStore = MvAlleleStoreTestUtil.newMvStoreWithData(value);
        AllelePropertiesDaoMvStore allelePropertiesDao = new AllelePropertiesDaoMvStore(mvStore);
        return new AllelePropertiesDaoAdapter(allelePropertiesDao);
    }

    @Test
    void assign() {
        FrequencyData frequencyData = FrequencyData.of("rs1101");

        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .geneSymbol("FGFR2")
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(1,1), "T", "A")
                .frequencyData(frequencyData)
                .build();

        System.out.println(variantEvaluation);
        System.out.println(frequencyData);

        Variant variant = buildVariant(1, 12345, "A", "T");
        AlleleProto.AlleleKey key = AlleleProtoAdaptor.toAlleleKey(variant);
        AlleleProto.AlleleProperties properties = AlleleProto.AlleleProperties.newBuilder().setRsId("rs54321")
                .putProperties("KG", 0.04f)
                .putProperties("ESP_AA", 0.003f)
                .build();
        FrequencyDao instance = newInstanceWithData(ImmutableMap.of(key, properties));
        assertThat(instance.getFrequencyData(variant),
                equalTo(FrequencyData.of("rs54321",
                        Frequency.of(FrequencySource.THOUSAND_GENOMES, 0.04f),
                        Frequency.of(FrequencySource.ESP_AFRICAN_AMERICAN, 0.003f))));


    }
}
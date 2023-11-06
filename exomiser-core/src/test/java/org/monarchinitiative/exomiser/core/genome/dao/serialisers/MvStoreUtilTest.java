/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2018 Queen Mary University of London.
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

package org.monarchinitiative.exomiser.core.genome.dao.serialisers;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.h2.mvstore.type.ObjectDataType;
import org.h2.mvstore.type.StringDataType;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.exomiser.core.proto.AlleleProto.AlleleKey;
import org.monarchinitiative.exomiser.core.proto.AlleleProto.AlleleProperties;
import org.monarchinitiative.exomiser.core.proto.AlleleProto.ClinVar;

import java.util.Map;

import static org.hamcrest.CoreMatchers.*;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class MvStoreUtilTest {

    @Test
    public void openAlleleMVMap() {
        MVStore mvStore = new MVStore.Builder().open();
        MVMap<AlleleKey, AlleleProperties> map = MvStoreUtil.openAlleleMVMap(mvStore);
        assertThat(map.isEmpty(), is(true));
        assertThat(mvStore.hasMap("alleles"), is(true));
    }

    @Test
    public void alleleMapBuilder() {
        MVMap.Builder<AlleleKey, AlleleProperties> alleleMapBuilder = MvStoreUtil.alleleMapBuilder();
        assertThat(alleleMapBuilder.getKeyType(), equalTo(AlleleKeyDataType.INSTANCE));
        assertThat(alleleMapBuilder.getValueType(), equalTo(AllelePropertiesDataType.INSTANCE));
    }

    @Test
    public void openClinVarMVMap() {
        MVStore mvStore = new MVStore.Builder().open();
        MVMap<AlleleKey, ClinVar> map = MvStoreUtil.openClinVarMVMap(mvStore);
        assertThat(map.isEmpty(), is(true));
        assertThat(mvStore.hasMap("clinvar"), is(true));
    }

    @Test
    public void clinvarMapBuilder() {
        MVMap.Builder<AlleleKey, ClinVar> alleleMapBuilder = MvStoreUtil.clinVarMapBuilder();
        assertThat(alleleMapBuilder.getKeyType(), equalTo(AlleleKeyDataType.INSTANCE));
        assertThat(alleleMapBuilder.getValueType(), equalTo(ClinVarDataType.INSTANCE));
    }
    @Test
    public void clinVarGeneStatsBuilder() {
        MVMap.Builder<String, ClinVarGeneStats> geneStatsMapBuilder = MvStoreUtil.geneStatsMapBuilder();
        assertThat(geneStatsMapBuilder.getKeyType(), equalTo(StringDataType.INSTANCE));
        assertThat(geneStatsMapBuilder.getValueType(), instanceOf(ObjectDataType.class));
    }

    @Test
    public void openClinVarGeneStatsMVMap(){
        MVStore mvStore = new MVStore.Builder().open();
        MVMap<String, ClinVarGeneStats> map = MvStoreUtil.openGeneStatsMVMap(mvStore);
        map.put("TEST", new ClinVarGeneStats("TEST", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 1))));
        assertThat(map.isEmpty(), is(false));
        assertThat(mvStore.hasMap("genestats"), is(true));
        assertThat(map.get("TEST"), equalTo(Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 1))));
    }
}
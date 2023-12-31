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

package org.monarchinitiative.exomiser.core.genome.dao;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.monarchinitiative.exomiser.core.model.SvMetaType;
import org.monarchinitiative.svart.VariantType;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;

class SvMetaTypeTest {

    @ParameterizedTest
    @CsvSource({
            "INS, INS,           true",
            "INS, DUP,           true",
            "DUP, CNV_GAIN,      true",
            "DEL, DEL,           true",
            "DEL, CNV_LOSS,      true",
            "INS, CNV_GAIN,      true",
            "TRA, BND,           true",
            "INS, INS_ME,        false",
            "INS, DEL,           false",
            "INS, CNV_LOSS,      false",
            "INS, BND,           false",
            "CNV, INS,           true",
            "CNV, DUP,           true",
            "CNV, CNV_GAIN,      true",
            "CNV, DEL,           true",
            "DEL, CNV,           true",
            "CNV, CNV_LOSS,      true",
            "CNV, INV,           false",
            "CNV, BND,           false",
    })
    public void testEquivalentTypes(VariantType a, VariantType b, boolean expected) {
        assertThat(SvMetaType.isEquivalent(a, b), equalTo(expected));
    }
}
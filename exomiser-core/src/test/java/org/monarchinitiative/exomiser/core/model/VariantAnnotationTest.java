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

package org.monarchinitiative.exomiser.core.model;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.sameInstance;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class VariantAnnotationTest {

    @Test
    public void empty() {
        VariantAnnotation instance = VariantAnnotation.of("", "", VariantEffect.SEQUENCE_VARIANT, List.of());
        assertThat(instance, sameInstance(VariantAnnotation.empty()));
    }

    @Test
    public void notEmpty() {
        VariantAnnotation instance = VariantAnnotation.of("GENE1", "HGNC:12345", VariantEffect.SEQUENCE_VARIANT, List.of());
        assertThat(instance.getGeneSymbol(), equalTo("GENE1"));
        assertThat(instance.getGeneId(), equalTo("HGNC:12345"));
        assertThat(instance.getVariantEffect(), equalTo(VariantEffect.SEQUENCE_VARIANT));
        assertThat(instance.getTranscriptAnnotations(), equalTo(List.of()));
        assertThat(instance.hasTranscriptAnnotations(), equalTo(false));
    }

}
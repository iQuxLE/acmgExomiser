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

import com.fasterxml.jackson.annotation.JsonIgnore;
import de.charite.compbio.jannovar.annotation.VariantEffect;

import java.util.List;

/**
 *
 * @since 13.0.0
 */
public interface VariantAnnotations {

    public String getGeneSymbol();

    public String getGeneId();

    public VariantEffect getVariantEffect();

    public List<TranscriptAnnotation> getTranscriptAnnotations();

    public boolean hasTranscriptAnnotations();

    @JsonIgnore
    public default boolean isNonCodingVariant() {
        return VariantEffectUtility.isNonCodingVariant(getVariantEffect());
    }

    /**
     * @return true if the Variant is in a coding region, otherwise false
     * @since 12.0.0
     */
    @JsonIgnore
    public default boolean isCodingVariant() {
        return VariantEffectUtility.affectsCodingRegion(getVariantEffect());
    }

}

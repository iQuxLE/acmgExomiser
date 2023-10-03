package org.monarchinitiative.exomiser.core.genome.dao;

import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.svart.GenomicInterval;
import org.monarchinitiative.svart.GenomicVariant;

import java.util.Map;


/**
 * Interface for providing {@link ClinVarData} about a {@link Variant}.
 *
 * @since 14.0.0
 */
public interface ClinVarDao {

    ClinVarData getClinVarData(Variant variant);

    Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval);

}

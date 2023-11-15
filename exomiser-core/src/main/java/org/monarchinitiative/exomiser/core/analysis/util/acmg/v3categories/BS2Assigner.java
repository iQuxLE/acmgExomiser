package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.genome.dao.SvFrequencyDao;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;

import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BS2;

public class BS2Assigner {

    private static final int RECESSIVE_X_LINKED_THRESHOLD = 2;
    private static final int DOMINANT_THRESHOLD = 5;

    private final VariantDataService variantDataService;

    private final SvFrequencyDao svFrequencyDao;

    public BS2Assigner(VariantDataService variantDataService, SvFrequencyDao svFrequencyDao) {
        this.variantDataService = variantDataService;
        this.svFrequencyDao = svFrequencyDao;
    }

    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation, ModeOfInheritance modeOfInheritance) {
        // method NOT IMPLEMENTED yet (mauybe with alleleKey information ?
        // needs AlleleCount information (is in SvFrequency)
//        svFrequencyDao.getFrequencyData(variantEvaluation);
        int alleleCount = variantDataService.getAlleleCount(variantEvaluation.getVariant());

        switch (modeOfInheritance) {
            case AUTOSOMAL_RECESSIVE:
            case X_RECESSIVE:
                if (alleleCount > RECESSIVE_X_LINKED_THRESHOLD) {
                    acmgEvidenceBuilder.add(BS2, AcmgCriterion.Evidence.STRONG);
                }
                break;
            case AUTOSOMAL_DOMINANT:
                if (alleleCount > DOMINANT_THRESHOLD) {
                    acmgEvidenceBuilder.add(BS2, AcmgCriterion.Evidence.STRONG);
                }
                break;
            default:
                // For other inheritance patterns, or if the inheritance pattern is unknown ?
                // the BS2 criterion might not be applicable or might require more complex logic?
                break;
        }
    }

    public void testAssign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation, ModeOfInheritance modeOfInheritance){
        svFrequencyDao.getFrequencyData(variantEvaluation);
        // test this and see how we get to AN or AC for a given variant
    }
}

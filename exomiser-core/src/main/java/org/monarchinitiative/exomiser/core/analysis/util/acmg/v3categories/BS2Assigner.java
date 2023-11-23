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

    //BS2 - "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked
    // (hemizygous) disorder, with full penetrance expected at an early age" - needs het/hom counts from gnomAD
    // (these should be in the newer v14 variant store, but I need to double-check)\
    //Varsome:
    //We first determine the mode of inheritance of the gene, then compares the allele count (see allele frequency for quality checks) to the corresponding threshold:
    //
    //recessive or X-linked genes: allele count greater than 2,
    //dominant genes: allele count greater than 5.

    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation, ModeOfInheritance modeOfInheritance) {
        // method NOT IMPLEMENTED yet (mauybe with alleleKey information ?
        // needs AlleleCount information (is in SvFrequency)
//        svFrequencyDao.getFrequencyData(variantEvaluation);
//        int alleleCount = variantDataService.getAlleleCount(variantEvaluation.getVariant());
        int alleleCount  = 2;

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

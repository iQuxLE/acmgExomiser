package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;

public class PM2Assigner {

    private final VariantDataService variantDataService;

    // PM2 - "Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes
    // Project, or Exome Aggregation Consortium" - This one just need refining for the recessive case, again using het/hom counts from gnomAD.

    //Varsome:
    //PM2
    //Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project,
    // or Exome Aggregation Consortium. (Pathogenic, Moderate)
    //
    //We first established the gene's mode of inheritance.
    //
    //The rule will trigger if the allele frequency is not found in GnomAD, with valid coverage, or:
    //
    //For dominant genes (including X-Linked and AD/AR) we check that the allele count is less than 5.
    //For recessive genes (AR): the rule will trigger if the homozygous allele count is less than 2. Alternatively
    // the rule also checks whether the allele frequency is less than 0.0001.
    //For mitochondrial variants, rule PM2 will trigger if the allele frequency is below 2e-05 per ClinGen Guidelines.

    //In line with the SVI Recommendation for Absence/Rarity (PM2) - Version 1.0, rule PM2 always triggers with strength
    // supporting. Rule PM2 may be phased out altogether in future.

    public PM2Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }

//    public void assign(AcmgEvidence.Builder evidenceBuilder, VariantEvaluation variantEvaluation, ModeOfInheritance modeOfInheritance){
//        if (modeOfInheritance.isDominant()){
//
//        } else if (modeOfInheritance.isRecessive()) {
//
//        }
//        var x = modeOfInheritance.getAbbreviation()
//    }
}

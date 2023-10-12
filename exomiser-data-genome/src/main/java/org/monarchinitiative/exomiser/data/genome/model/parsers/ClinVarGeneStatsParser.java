package org.monarchinitiative.exomiser.data.genome.model.parsers;

import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStatistics;
import org.monarchinitiative.exomiser.data.genome.model.Allele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

import static org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData.ClinSig.*;


public class ClinVarGeneStatsParser {

    public class ClinVarAlleleParser extends VcfAlleleParser {

        private static final Logger logger = LoggerFactory.getLogger(org.monarchinitiative.exomiser.data.genome.model.parsers.ClinVarAlleleParser.class);

        @Override
        List<Allele> parseInfoField(List<Allele> alleles, String info) {
            ClinVarData clinVarData = parseClinVarData(info);
            for (Allele allele : alleles) {
                if (!clinVarData.isEmpty()) {
                    allele.setClinVarData(clinVarData);
                }
                logger.debug("{}", allele);
            }
            return alleles;
        }

        /**
         *
         * @param info
         * @return
         */
        private ClinVarGeneStatistics parseClinVarData(String info) {
//        ##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID"> - get to the web record using: https://www.ncbi.nlm.nih.gov/clinvar/?term=99222[alleleid]
//        ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance for this single variant">
//        ##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.">

            ClinVarGeneStatistics.Builder geneStatsBuilder = ClinVarGeneStatistics.builder();
            String [] fields = info.split(";");
            for (String field : fields) {
                String[] keyValue = field.split("=");
                String key = keyValue[0];
                String value = keyValue[1];
                switch (key) {
                    case "GENEINFO":
                        String geneSymbol = value.split(":")[0];
                        geneStatsBuilder.geneSymbol(geneSymbol);
                        break;
                    case "CLNSIG":
                        String[] clinsigs = value.split(",_");
                        ClinVarData.ClinSig primary = parseClinSig(clinsigs[0]);
                        Set<ClinVarData.ClinSig> secondary = parseSecondaryClinSig(clinsigs);
                        geneStatsBuilder.primaryInterpretation(primary);
                        geneStatsBuilder.secondaryInterpretations(secondary);
                        break;
                    case "CLNREVSTAT":
                        //CLNREVSTAT criteria_provided,_conflicting_interpretations, criteria_provided,_multiple_submitters,_no_conflicts, criteria_provided,_single_submitter, no_assertion_criteria_provided, no_assertion_provided, no_interpretation_for_the_single_variant, practice_guideline, reviewed_by_expert_panel
                        //CLNREVSTAT counts: criteria_provided,_conflicting_interpretations=12678, criteria_provided,_multiple_submitters,_no_conflicts=34967, criteria_provided,_single_submitter=197277, no_assertion_criteria_provided=34308, no_assertion_provided=10980, no_interpretation_for_the_single_variant=500, practice_guideline=23, reviewed_by_expert_panel=8786
                        geneStatsBuilder.reviewStatus(value);
                        break;
                    case "CLNSIGINCL":
                        Map<String, ClinVarData.ClinSig> includedAlleles = parseIncludedAlleles(value);
                        geneStatsBuilder.includedAlleles(includedAlleles);
                        break;
                    default:
                        break;
                }
            }
            return clinVarBuilder.build();
        }

        private Map<String, ClinVarData.ClinSig> parseIncludedAlleles(String value) {
            //15127:other|15128:other|15334:Pathogenic|
            Map<String, ClinVarData.ClinSig> includedAlleles = new HashMap<>();
            String[] incls = value.split("\\|");
            for (String inc : incls) {
                String[] fields = inc.split(":");
                if (fields.length == 2) {
                    includedAlleles.put(fields[0], parseClinSig(fields[1]));
                }
            }
            return includedAlleles;
        }

        private Set<ClinVarData.ClinSig> parseSecondaryClinSig(String[] clinsigs) {
            if (clinsigs.length > 1) {
                Set<ClinVarData.ClinSig> secondaryClinSigs = EnumSet.noneOf(ClinVarData.ClinSig.class);
                for (int i = 1; i < clinsigs.length; i++) {
                    secondaryClinSigs.add(parseClinSig(clinsigs[i]));
                }
                return secondaryClinSigs;
            }
            return Collections.emptySet();
        }

        private ClinVarData.ClinSig parseClinSig(String clinsig) {
            // Unique CLNSIG counts
            // Affects=100, Benign=23963, Benign/Likely_benign=10827, Conflicting_interpretations_of_pathogenicity=12784,
            // Likely_benign=52064, Likely_pathogenic=15127, Pathogenic=46803, Pathogenic/Likely_pathogenic=3278,
            // Uncertain_significance=120418, association=148, drug_response=290, not_provided=10980, other=1796, protective=30,
            // risk_factor=411
            switch (clinsig) {
                case "Uncertain_significance":
                    return UNCERTAIN_SIGNIFICANCE;
                case "Benign":
                    return BENIGN;
                case "Benign/Likely_benign":
                    return BENIGN_OR_LIKELY_BENIGN;
                case "Likely_benign":
                    return LIKELY_BENIGN;
                case "Conflicting_interpretations_of_pathogenicity":
                    return CONFLICTING_PATHOGENICITY_INTERPRETATIONS;
                case "Likely_pathogenic":
                    return LIKELY_PATHOGENIC;
                case "Pathogenic/Likely_pathogenic":
                    return PATHOGENIC_OR_LIKELY_PATHOGENIC;
                case "Pathogenic":
                    return PATHOGENIC;
                case "Affects":
                    return AFFECTS;
                case "association":
                    return ASSOCIATION;
                case "drug_response":
                    return DRUG_RESPONSE;
                case "other":
                    return OTHER;
                case "protective":
                    return PROTECTIVE;
                case "risk_factor":
                    return RISK_FACTOR;
                case "not_provided":
                default:
                    return NOT_PROVIDED;
            }
        }
    }

}

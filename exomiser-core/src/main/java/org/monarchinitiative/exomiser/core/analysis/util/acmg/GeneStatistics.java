package org.monarchinitiative.exomiser.core.analysis.util.acmg;

import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import java.util.HashMap;
import java.util.Map;

/*
    data class for geneStatistics
    = give geneName when open
    = get example:
    stop_gained {P=1,       frameshift: { P=1,
                 LP=2,                    LP=2,
                 B=3,                     B=3,
                 LB=4,                    LB=4,
                 ...}                       ....}
 */

public class GeneStatistics {

    private String geneName;
    private Map<String, Map<ClinVarData.ClinSig, Integer>> consequenceClinicalSignificanceCounts;

    private GeneStatistics(String geneName, Map<String, Map<ClinVarData.ClinSig, Integer>> consequenceClinicalSignificanceCounts){
        this.geneName = geneName;
        this.consequenceClinicalSignificanceCounts = consequenceClinicalSignificanceCounts;

    }

    public Map<String, Map<ClinVarData.ClinSig, Integer>> getConsequenceClinicalSignificanceCounts(){
        return consequenceClinicalSignificanceCounts;
    }

    public String getGeneName(){
        return geneName;
    }

    public void countConsequence(String consequence, ClinVarData.ClinSig clinicalSignificance) {
     consequenceClinicalSignificanceCounts.putIfAbsent(consequence, new HashMap<>());
     Map<ClinVarData.ClinSig, Integer> clinicalSignificanceCounts = consequenceClinicalSignificanceCounts.get(consequence);

    int count = clinicalSignificanceCounts.getOrDefault(clinicalSignificance, 0);
    clinicalSignificanceCounts.put(clinicalSignificance, count+1);

    //countConsequence(stop_gained

    }

}
